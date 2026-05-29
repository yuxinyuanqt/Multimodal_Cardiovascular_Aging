#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import numpy as np
import optuna
from catboost import CatBoostRegressor
import pandas as pd
from catboost import Pool
from sklearn.metrics import r2_score


# In[ ]:


def objective(trial):
    # Hyperparameter suggestions
    param = {
        'learning_rate': trial.suggest_float('learning_rate', 1e-5, 1e-1, log=True),
        'depth': trial.suggest_int('depth', 4, 12),
        'l2_leaf_reg': trial.suggest_int('l2_leaf_reg', 1, 20),
        'boosting_type': trial.suggest_categorical("boosting_type", ["Ordered", "Plain"]),
        'bootstrap_type': trial.suggest_categorical("bootstrap_type", ['Bayesian', 'Bernoulli','MVS']),
        'cat_features': [0],  # Specify 'Sex' index as categorical feature
        'loss_function': 'RMSE',
        'eval_metric': "MAE",
        'metric_period' : 5,
        'task_type': "GPU",
        'devices': '1',
        'silent': True,
        'random_seed': 97,
        'early_stopping_rounds' : 100,
        'use_best_model': True
    }
    
    if param["bootstrap_type"] == "Bayesian":
        param["bagging_temperature"] = trial.suggest_float("bagging_temperature", 0, 10)
    elif param["bootstrap_type"] == "Bernoulli":
        param["subsample"] = trial.suggest_float("subsample", 0.1, 1)
        
    
    
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    val_r2_scores = []
    fold_idx=0   
    for train_index, val_index in kf.split(X):
        
        fold_idx += 1
        X_train, X_val = X[train_index], X[val_index]
        y_train, y_val = y[train_index], y[val_index]
        
        scaler = StandardScaler()
        X_train_scaled = X_train.copy()
        X_val_scaled = X_val.copy()

        # Scale only the continuous columns using column indices
        X_train_scaled[:, continuous_columns_indices] = scaler.fit_transform(X_train[:, continuous_columns_indices])
        X_val_scaled[:, continuous_columns_indices] = scaler.transform(X_val[:, continuous_columns_indices])

        
        X_train_scaled_df = pd.DataFrame(X_train_scaled, columns=all_columns_names)
        X_val_scaled_df = pd.DataFrame(X_val_scaled, columns=all_columns_names)
        
        X_train_scaled_df['Sex'] = X_train_scaled_df['Sex'].astype('int')
        X_val_scaled_df['Sex'] = X_val_scaled_df['Sex'].astype('int')
        
        train_pool = Pool(X_train_scaled_df, y_train, cat_features=[sex_column_index])
        val_pool = Pool(X_val_scaled_df, y_val, cat_features=[sex_column_index])
        
              
        # Fit the model
        model = CatBoostRegressor(**param)
        model.fit(train_pool, eval_set=val_pool)
        
        # Predict and evaluate
        val_pred = model.predict(val_pool)
        
        if np.isclose(np.std(val_pred), 0) or np.any(np.isnan(val_pred)):
            val_r2 = 0
        else:
            val_r2 = r2_score(y_val, val_pred)
            
        val_r2_scores.append(val_r2)
        
        ## Optuna pruning: report intermediate score and stop trial early if it underperforms
        trial.report(np.mean(val_r2_scores), fold_idx)
        if trial.should_prune():
            raise optuna.TrialPruned()
            

    return np.mean(val_r2_scores)


# In[ ]:


if __name__ == '__main__':
    
    df_catboost = pd.read_csv("/path_to_the_file/df_nocvd_cohort_train_test.csv")
    
    train = df_catboost[df_catboost['Center'] == 0].copy() # Cheadle Center 

    X = train.drop(columns=['eid','age_imaging_derived','Center']).values
    y = train['age_imaging_derived'].values
    
    # Identify columns for model training:
    sex_column_index = 0  # Categorical column index for 'sex'
    continuous_columns_indices = list(range(1, X.shape[1])) # Indices of continuous columns
    all_columns_names= train.columns.drop(['eid','age_imaging_derived','Center']).tolist() #Names of all columns used for training

    # Set Optuna logging level to INFO to show progress of trials
    optuna.logging.set_verbosity(optuna.logging.INFO)
    
    # Define the sampler for hyperparameter selection
    # The seed ensures reproducibility of the sampling process
    sampler = optuna.samplers.TPESampler(seed=97)
    
    # Define the pruner to stop unpromising trials early based on intermediate results
    # MedianPruner stops trials if their performance is below the median of completed trials
    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=5,   # allow the first 5 trials to complete before pruning begins
        n_warmup_steps=5,     # let each trial run for 5 steps/folds before checking for pruning
        interval_steps=1      # check pruning condition after every fold after warmup
        )
    
    # Create an Optuna study to optimize the objective function
    # direction="maximize" ensures Optuna tries to maximize the average R² score of the valiodation folds
    study = optuna.create_study(
            direction="maximize",
            sampler=sampler,
            pruner=pruner)

    # Run the hyperparameter optimization
    study.optimize(objective, n_trials=200)

    # Print the best hyperparameters and best score
    print("Best hyperparameters:", study.best_params)
    print("Best R^2 score:", study.best_value)
    
    # Convert the trial results into a pandas DataFrame   
    df_results = study.trials_dataframe(attrs=("number",          #trial number
                                "value",          # objective value (mean R2)
                                "params"         # hyperparameters for each trial
                             )) 
    
    df_results.to_csv("/path_to_the_file/optuna_results_catboost.csv", index=False)

