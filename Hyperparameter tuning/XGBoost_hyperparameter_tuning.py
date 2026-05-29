#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os 
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import optuna
from sklearn.preprocessing import StandardScaler
import xgboost as xgb


# In[ ]:


def objective(trial):
        
    # parameter ranges to search    
    param = {
            'subsample': trial.suggest_float('subsample', 0.1, 1.0),
            'colsample_bytree': trial.suggest_float('colsample_bytree', 0.1, 1.0),
            'learning_rate': trial.suggest_float('learning_rate', 1e-5, 1, log=True),
            'max_depth': trial.suggest_int('max_depth', 3, 12),
            'min_child_weight': trial.suggest_float('min_child_weight', 1, 100),
            'gamma': trial.suggest_float('gamma', 0, 1),
            'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1, log=True),
            'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1, log=True),
            'n_estimators': trial.suggest_int('n_estimators', 50, 1000),# default:100
            'eval_metric' : 'mae',
            'tree_method': 'hist',
            'device':'cuda',
             'seed' : 97,
             'verbosity':0
        }

                        
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    val_r2_scores = []
    fold_idx=0
    for train_index, val_index in kf.split(X):
        fold_idx += 1
        X_train, X_val = X[train_index], X[val_index]
        y_train, y_val = y[train_index], y[val_index]
        
        scaler = StandardScaler()
        continuous_columns_indices = list(range(1, X_train.shape[1])) 

        X_train_scaled = X_train.copy()
        X_val_scaled = X_val.copy()
    
        X_train_scaled[:, continuous_columns_indices] = scaler.fit_transform(X_train[:, continuous_columns_indices])
        X_val_scaled[:, continuous_columns_indices] = scaler.transform(X_val[:, continuous_columns_indices])
        
        es = xgb.callback.EarlyStopping(
        rounds=100,
        save_best=True,
        maximize=False,
        metric_name="mae")
        
        model=xgb.XGBRegressor(**param,callbacks=[es])
        
        model.fit(X_train_scaled,
                  y_train,
                  eval_set = [(X_val_scaled, y_val)],
                  verbose=False
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
    
    df_xgboost = pd.read_csv("/path_to_the_file/df_nocvd_cohort_train_test.csv")
    
    train = df_xgboost[df_xgboost['Center'] == 0].copy() 
  
    X = train.drop(columns=['eid','age_imaging_derived','Center']).values
    y = train['age_imaging_derived'].values

    optuna.logging.set_verbosity(optuna.logging.INFO)
    sampler = optuna.samplers.TPESampler(seed=97)
    pruner = optuna.pruners.MedianPruner(
        n_startup_trials=5,   
        n_warmup_steps=5,    
        interval_steps=1
        )


    study = optuna.create_study(
            direction="maximize", ## maximize R2 value
            sampler=sampler,
            pruner=pruner)

    study.optimize(objective, n_trials=200)

    print("Best hyperparameters:", study.best_params)
    print("Best R^2 score:", study.best_value)
