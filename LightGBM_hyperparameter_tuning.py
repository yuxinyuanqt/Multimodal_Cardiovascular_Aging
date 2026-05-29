#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os 
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import lightgbm as lgb
from sklearn.preprocessing import StandardScaler
import optuna


# In[ ]:


def objective(trial):
    # Hyperparameter suggestions
    param = {
        'boosting_type': 'gbdt',
        'num_leaves': trial.suggest_int('num_leaves', 2, 256),
        'max_depth': trial.suggest_int('max_depth', 1, 4),
        'subsample': trial.suggest_float('subsample', 0.1, 1.0),
        'learning_rate': trial.suggest_float('learning_rate', 1e-5, 1, log=True),
        'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1, log=True),
        'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1, log=True),
        'colsample_bytree': trial.suggest_float('colsample_bytree', 0.1, 1.0),
        'min_child_samples': trial.suggest_int('min_child_samples', 5, 500),
        'min_child_weight': trial.suggest_float('min_child_weight', 1e-5, 100, log=True),
        'eval_metric' : 'mae',
        'device' : 'gpu',
        'gpu_platform_id' : 0,
        'gpu_device_id' : 0,
        'random_state' : 97, 
        'verbosity': -1
    }
  
    
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    val_r2_scores = []
    fold_idx=0
    for train_index, val_index in kf.split(X):
        
        fold_idx += 1
        X_train, X_val = X[train_index], X[val_index]
        y_train, y_val = y[train_index],y[val_index]

        scaler = StandardScaler()
        
        X_train_scaled = X_train.copy()
        X_val_scaled = X_val.copy()

        X_train_scaled[:, continuous_columns_indices] = scaler.fit_transform(X_train[:, continuous_columns_indices])
        X_val_scaled[:, continuous_columns_indices] = scaler.transform(X_val[:, continuous_columns_indices]) 
        
        model = lgb.LGBMRegressor(**param)
        callbacks = [lgb.early_stopping(100, verbose=False)]
        model.fit(
                X=X_train_scaled,
                y=y_train,
                eval_set=(X_val_scaled, y_val),
                callbacks=callbacks,
            )
        
        val_pred = model.predict(X_val_scaled)
 
        if np.isclose(np.std(val_pred), 0)or np.any(np.isnan(val_pred)):
                val_r2 = 0
        else:
                val_r2 = r2_score(y_val, val_pred)
              
        val_r2_scores.append(val_r2)
        
        trial.report(np.mean(val_r2_scores), fold_idx)
        if trial.should_prune():
            raise optuna.TrialPruned()
            
    
    return np.mean(val_r2_scores)


# In[ ]:


if __name__ == '__main__':
    
    df_lgbm = pd.read_csv("/path_to_the_file/df_nocvd_cohort_train_test.csv")
    
    train = df_lgbm[df_lgbm['Center'] == 0].copy() 

    X = train.drop(columns=['eid','age_imaging_derived','Center']).values
    y = train['age_imaging_derived'].values
    
    sex_column_index = 0  
    continuous_columns_indices = list(range(1, X.shape[1])) 
    all_columns_names= train.columns.drop(['eid','age_imaging_derived','Center']).tolist() 
   
    optuna.logging.set_verbosity(optuna.logging.INFO)
    sampler = optuna.samplers.TPESampler(seed=97)
    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=5,  
        n_warmup_steps=5,    
        interval_steps=1      
        )

    study = optuna.create_study(
            direction="maximize", 
            sampler=sampler,
            pruner=pruner)

    study.optimize(objective, n_trials=200)

    print("Best hyperparameters:", study.best_params)
    print("Best R^2 score:", study.best_value)
    