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
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


# In[ ]:


class EarlyStopping:
   
    def __init__(self, patience, verbose=True, save_model=False, filepath=None):
        """
        Args:
            patience (int, optional): Number of epochs with no improvement to wait before stopping.
            verbose (bool, optional): Print information about stopping.
            save_model (bool, optional): Whether to save the model when early stopping occurs. 
            filepath (str, optional): Path to save the model file.
        """
        self.patience = patience
        self.verbose = verbose
        self.best_val_loss = None
        self.best_val_r2 = None
        self.best_epoch = None  
        self.early_stop = False
        self.wait_count = 0  # keeps track of the number of epochs since the last validation loss improvement
        self.save_model = save_model
        self.filepath = filepath  

    def __call__(self, val_loss,val_r2, model, epoch): 
        """
        Checks if validation loss has improved and updates stopping criteria. Optionally saves the model.

        Args:
            val_loss: Current validation loss.
            val_r2 : Current validation R2 metric.
            model (nn.Module): The PyTorch model to save (if `save_model` is True).
            epoch (int): Current epoch number.
        """ 
        if self.best_val_loss is None or val_loss < self.best_val_loss:
        
                self.best_val_loss = val_loss
                self.best_val_r2 = val_r2
                self.best_epoch = epoch  
                self.wait_count = 0
                if self.save_model:
                    torch.save({'model_state_dict': model.state_dict()}, self.filepath)  
        else:
            self.wait_count += 1
            if self.wait_count >= self.patience:
                self.early_stop = True
                if self.verbose:
                    print("Early stopping triggered!")


# In[ ]:

# Define a fully connected deep neural network (DNN) 
class Dnn(nn.Module):
    """
    A customizable feedforward deep neural network with multiple hidden layers,
    batch normalization, dropout, and configurable activation functions.

    Parameters:
    -----------
    hid_dims : list of int
        Number of neurons in each hidden layer, e.g., [128, 64, 32]
    p_drop : float
        Dropout probability applied after each hidden layer
    activation : str
        Activation function to use: "ReLU", "ReLU6", or "leaky_ReLU"
    input_size : int, default=57
        Number of input features
    out_dims : int, default=1
        Number of output units (e.g., 1 for age prediction)
    """
    def __init__(self,hid_dims,p_drop,activation, input_size=57, out_dims=1):
        super(Dnn, self).__init__()
        
        n_layers = len(hid_dims)
        self.feature_extractor = nn.Sequential()
        
        if activation=="ReLU":
            act_fun = nn.ReLU()
        elif activation == "ReLU6":
            act_fun = nn.ReLU6()
        elif activation == "leaky_ReLU":
            act_fun = nn.LeakyReLU()
            
        for i in range(n_layers):
            in_channel = input_size if i == 0 else hid_dims[i - 1]
            out_channel = hid_dims[i]
            self.feature_extractor.add_module('lin_%d' % i, self.fully_connected_layers(in_channel, out_channel,act_fun,p_drop))
            
        self.pred_head = nn.Linear(hid_dims[-1], out_dims)

    @staticmethod
    def fully_connected_layers(in_channel, out_channel,act_fun,p_drop):
        return nn.Sequential(
            nn.Linear(in_channel, out_channel),
            nn.BatchNorm1d(out_channel, affine=False),
            act_fun,
            nn.Dropout(p_drop)
        )

    def forward(self, input):
        x_f = self.feature_extractor(input)
        out = self.pred_head(x_f)
        return out


# In[ ]:


# Training process
def Train_Dnn(X, y,kf,lr,weight_decay,hid_dims,activation,dropout,device):
    """
    Train a deep neural network using k-fold cross-validation and return the average R² score.

    Parameters
    ----------
    X : np.ndarray
        Input features, shape (n_samples, n_features)
    y : np.ndarray
        Target variable (chronological age), shape (n_samples,)
    kf : sklearn.model_selection._BaseKFold
        K-fold cross-validation splitter
    lr : float
        Learning rate for Adam optimizer
    weight_decay : float
        L2 regularization parameter
    hid_dims : list of int
        Hidden layer sizes, e.g., [128, 64, 32]
    activation : str
        Activation function for hidden layers ("ReLU", "ReLU6", "leaky_ReLU")
    dropout : float
        Dropout probability
    device : str or torch.device
        Device to run the model on ("cpu" or "cuda")
    
    Returns
    -------
    avg_val_r2 : float
        Average R² score across all cross-validation folds
    """

    val_r2_scores = []
    
    for train_index, val_index in kf.split(X):
        
        X_train, X_val = X[train_index], X[val_index]
        y_train, y_val = y[train_index], y[val_index]
       
        scaler = StandardScaler()
        
        X_train_scaled = X_train.copy()
        X_val_scaled = X_val.copy()
       

        X_train_scaled[:, continuous_columns_indices] = scaler.fit_transform(X_train[:, continuous_columns_indices])
        X_val_scaled[:, continuous_columns_indices] = scaler.transform(X_val[:, continuous_columns_indices])
        
     
        train_data = TensorDataset(torch.FloatTensor(X_train_scaled), torch.FloatTensor(y_train))
        train_loader = DataLoader(train_data, batch_size=256, shuffle=True)
        
        val_data = TensorDataset(torch.FloatTensor(X_val_scaled), torch.FloatTensor(y_val))
        val_loader = DataLoader(val_data,batch_size=256,shuffle=True)
        
        
        # Initialize the model
        model = Dnn(hid_dims=hid_dims,p_drop=dropout,activation=activation).to(device)
        criterion = nn.MSELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=lr,weight_decay=weight_decay)
        early_stopping = EarlyStopping(patience=20, verbose=True, save_model=False, filepath=None)

        # Train the model for up to 100 epochs
        for epoch in range(100):
            # Training loop
            model.train()
            train_loss = 0.0
            for inputs, labels in train_loader:
                inputs,labels = inputs.to(device),labels.to(device)
                optimizer.zero_grad()
                inputs=inputs.to(device)
                labels=labels.to(device)
                outputs = model(inputs).squeeze()
                loss = criterion(outputs, labels)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item()
                

            #Validation loop
            model.eval()
            val_loss = 0.0
            all_val_preds = []
            all_val_labels = []
            with torch.no_grad():
                for inputs, labels in val_loader:
                    inputs,labels = inputs.to(device),labels.to(device)
                    outputs = model(inputs).squeeze()
                    loss = criterion(outputs, labels)
                    val_loss += loss.item()

                    all_val_preds.append(outputs.cpu())
                    all_val_labels.append(labels.cpu())
                
               
            all_val_preds = torch.cat(all_val_preds).numpy()
            all_val_labels = torch.cat(all_val_labels).numpy()
            
            val_r2 = r2_score(all_val_labels, all_val_preds)
                        
            #print(f"Epoch {epoch}/{100 - 1}, Train Loss: {train_loss:.4f},  Val Loss: {val_loss:.4f}, Val R2 :{val_r2:.4f}")
            
            early_stopping(val_loss,val_r2,model,epoch)
            if early_stopping.early_stop:
                break
        
        best_epoch = early_stopping.best_epoch
        best_val_r2 = early_stopping.best_val_r2
        best_val_loss = early_stopping.best_val_loss
        
        val_r2_scores.append(best_val_r2)
        
        
            
    
    avg_val_r2 = np.mean(val_r2_scores)
    
    return avg_val_r2


# In[ ]:


def objective(trial):
    
    device = 'cuda'
    
    #hyperparameter suggestions
    num_layers = trial.suggest_int('num_layers', 3, 7)
    
    hidden_dims = {}
    choices = [32, 64, 128, 256, 512, 1024, 2048]
    hidden_dims = [trial.suggest_categorical(f'hidim_{i}',choices) for i in range(num_layers)] 
    
    dropout = trial.suggest_float('p_drop', 0.0, 0.9)
    activation=trial.suggest_categorical('act_fun',['ReLU','ReLU6','leaky_ReLU'])
    lr = trial.suggest_float('lr', 1e-5, 1e-2,log=True)
    weight_decay = trial.suggest_float('weight_decay', 1e-7, 1e-2,log=True)
    
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    
    avg_val_r2=Train_Dnn(X,y,kf,lr,weight_decay,hidden_dims,activation,dropout,device)
    
    return avg_val_r2


# In[ ]:


if __name__ == '__main__':
    
    torch.cuda.set_device(0)
    
    df_dnn = pd.read_csv("/path_to_the_file/df_nocvd_cohort_train_test.csv")
    
    train = df_dnn[df_dnn['Center'] == 0].copy() 
  
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
    