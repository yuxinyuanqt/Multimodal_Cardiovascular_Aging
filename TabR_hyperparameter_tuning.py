#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os 
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
import torch
import optuna
from pytorch_tabr import TabRRegressor
import warnings
import sys
import contextlib

# In[ ]:


@contextlib.contextmanager
def suppress_output():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
def restore_output():
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__


# In[ ]:


def objective(trial):
    
    #hyperparameter suggestions         
    lr=trial.suggest_float('lr', 5e-5, 1e-3,log=True)
    weight_decay=trial.suggest_float("weight_decay", 1e-6, 1e-3,log=True)
    encoder_n_blocks=trial.suggest_int("encoder_n_blocks", 0,1)
    predictor_n_blocks= trial.suggest_int("predictor_n_blocks", 1,2)
    d_main=trial.suggest_int("d_main", 96,400)
    d_multiplier=trial.suggest_float("d_multiplier", 0.01, 100.0)
    context_dropout=trial.suggest_float("context_dropout",0, 0.8)
    dropout0=trial.suggest_float("dropout0",0, 0.8)
    context_size=trial.suggest_int("context_size", 96,400)
    
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    fold_idx=0
    val_r2_scores = []
    
    for train_index, val_index in kf.split(X):

        X_train, X_val = X[train_index], X[val_index]
        y_train, y_val = y[train_index], y[val_index]
        
        
        scaler = StandardScaler()
        
        X_train_scaled = X_train.copy()
        X_val_scaled = X_val.copy()

        X_train_scaled[:, continuous_columns_indices] = scaler.fit_transform(X_train[:, continuous_columns_indices])
        X_val_scaled[:, continuous_columns_indices] = scaler.transform(X_val[:, continuous_columns_indices])
               
        model = TabRRegressor(
                bin_indices=[0],
                device_name="cuda",
                context_dropout=context_dropout,
                dropout0=dropout0,
                num_embeddings=None,
                optimizer_fn= torch.optim.AdamW,
                optimizer_params=dict(lr=lr,weight_decay=weight_decay),
                mixer_normalization='auto',
                d_main= d_main,
                encoder_n_blocks = encoder_n_blocks,
                predictor_n_blocks= predictor_n_blocks,
                d_multiplier=d_multiplier,
                normalization="LayerNorm",
                activation="ReLU",
                seed=97,
                verbose=False
                )
        
        with suppress_output():
            model.fit(
                    X_train_scaled, y_train.reshape(-1, 1),
                    eval_set=[(X_val_scaled, y_val.reshape(-1, 1))],
                    max_epochs=100,
                    patience=20,
                    batch_size=256
            )
        
        val_pred = model.predict(X_val_scaled)
        
        if np.isclose(np.std(val_pred), 0)or np.any(np.isnan(val_pred)):
                val_r2 = 0
        else:
                val_r2 = r2_score(y_val.flatten(), val_pred.flatten())
                         
        val_r2_scores.append(val_r2)
        trial.report(np.mean(val_r2_scores), fold_idx)
        if trial.should_prune():
            raise optuna.TrialPruned()
        

          
    
    torch.cuda.empty_cache()   
    
    return np.mean(val_r2_scores)
    
        
        


if __name__ == '__main__':

    torch.cuda.set_device(0)
    
    df_tabr = pd.read_csv("/path_to_the_file/df_nocvd_cohort_train_test.csv")
    
    train = df_tabr[df_tabr['Center'] == 0].copy() 
  
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
            direction="maximize", 
            sampler=sampler,
            pruner=pruner)

    study.optimize(objective, n_trials=100)

    print("Best hyperparameters:", study.best_params)
    print("Best R^2 score:", study.best_value)
