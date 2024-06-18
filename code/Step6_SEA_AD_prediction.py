#! /usr/bin/env /ua/choi267/anaconda3/bin/python3

import torch
import torch.nn as nn
import torch.nn.functional as F

import pandas as pd

from torch.utils.data import Dataset
from torch.utils.data.sampler import Sampler
from torch.utils.data import Dataset, DataLoader, random_split, SubsetRandomSampler
from torch.utils.data.dataset import ConcatDataset

import matplotlib.pyplot as plt
from matplotlib import pyplot as plt
import numpy as np

from sklearn.model_selection import KFold

import torch.optim as optim
import torch.optim.lr_scheduler as lr_scheduler
import argparse
import matplotlib
import torch.nn as nn

from tqdm import tqdm

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

from sklearn.metrics import mean_squared_error, r2_score

matplotlib.style.use('ggplot')


import copy

import sys
import pandas as pd

# Load your data
data1 = pd.read_csv('/ua/choi267/Oligo project/Step2/Interaction_median_1/sub_median_1_DEG_TFs_df.csv', header=0)
data2 = pd.read_csv('/ua/choi267/Oligo project/Step2/Interaction_median_1/sub_median_1_DEG_TGs_df.csv', header=0)

# Remove columns corresponding to "MEF2B" and "ZNF8" from data1
data1_filtered = data1.drop(columns=["MEF2B", "ZNF8"])
df3 = pd.concat([data1_filtered, data2.iloc[:, 403]], axis=1)

# Extract features and labels
data = df3.loc[:, df3.columns != df3.columns[df3.shape[1] - 1]].values
labels = df3.loc[:, df3.columns == df3.columns[df3.shape[1] - 1]].values

# Normalize the features using z-score normalization
scaler = StandardScaler()
data_normalized = scaler.fit_transform(data)

data = data_normalized

gene = "MBP"

path = "/ua/Dissertation_JC/Oligo/Step6/Validation/"
model_path = "/ua/Dissertation_JC/Oligo/Step6/Validation/" + gene + ".trained.pt"
history_path = "/ua/Dissertation_JC/Oligo/Step6/Validation" + gene +  ".history.csv"
train_path = "/ua/Dissertation_JC/Oligo/Step6/Validation/" + gene + "_train.png"
val_path = "/ua/Dissertation_JC/Oligo/Step6/Validation/" + gene + "_val.png"
MSE_path = "/ua/Dissertation_JC/Oligo/Step6/Validation/" + gene + "_MSE.csv"
prediction_path = "/ua/Dissertation_JC/Oligo/Step6/Validation/"+ gene + "_prediction.png"

# Separate data and labels
features = data.shape[1]

# Define your network
class Net(nn.Module):
    def __init__(self, features):
        super(Net, self).__init__()
        self.layers = nn.Sequential(
            nn.Linear(features, features * 2),
            nn.BatchNorm1d(features * 2),
            nn.LeakyReLU(),

            nn.Linear(features * 2, features * 3),
            nn.BatchNorm1d(features * 3),
            nn.LeakyReLU(),

            nn.Linear(features * 3, features * 4),
            nn.BatchNorm1d(features * 4),
            nn.LeakyReLU(),

            nn.Linear(features * 4, features * 3),
            nn.BatchNorm1d(features * 3),
            nn.LeakyReLU(),

            nn.Linear(features * 3, features * 2),
            nn.BatchNorm1d(features * 2),
            nn.LeakyReLU(),

            nn.Linear(features * 2, features),
            nn.BatchNorm1d(features),
            nn.LeakyReLU(),

            nn.Linear(features, features//2),
            nn.BatchNorm1d(features//2),
            nn.LeakyReLU(),

            nn.Linear(features//2, 1),
            nn.BatchNorm1d(1),
            nn.LeakyReLU()
        )

    def forward(self, x):
        x = x.float()
        x = self.layers(x)
        x = x.view(x.size(0), -1)
        return x

# Define your custom dataset
class CustomDataset(Dataset):
    def __init__(self, data, labels):
        self.data = torch.tensor(data, dtype=torch.float32)
        self.labels = torch.tensor(labels, dtype=torch.float32)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx], self.labels[idx]

# Define your fit and validate functions
def fit(model, dataloader, optimizer, criterion, device):
    model.train()
    running_loss = 0.0

    for batch_idx, (data, labels) in enumerate(tqdm(dataloader, total=len(dataloader))):
        optimizer.zero_grad()
        data, labels = data.to(device), labels.to(device)
        output = model(data)
        loss = criterion(output, labels)
        loss.backward()
        optimizer.step()
        running_loss += loss.item()

    train_loss = running_loss / len(dataloader)
    return train_loss

def validate(model, dataloader, criterion, device):
    model.eval()
    running_loss = 0.0

    with torch.no_grad():
        for data, labels in tqdm(dataloader, total=len(dataloader)):
            data = data.to(device)
            labels = labels.to(device)
            output = model(data)
            loss = criterion(output, labels)
            running_loss += loss.item()

    val_loss = running_loss / len(dataloader)
    return val_loss

# Set your device
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

# Define your parameters
batch_size = 32
epochs = 100
lr = 0.0001
k = 5  # Number of folds for cross-validation

# Initialize K-Fold cross-validation
kf = KFold(n_splits=k, shuffle=True, random_state=42)

# Define your loss function
criterion = nn.MSELoss()

# Split the data into train and validation (hold-out) sets
train_data, val_data, train_labels, val_labels = train_test_split(data, labels, test_size=0.25, random_state=42, shuffle=True)

# Custom datasets for hold-out validation data
val_dataset = CustomDataset(val_data, val_labels)
val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False)

# Initialize best validation loss
best_val_loss = float('inf')
best_model_weights = None

# Initialize history dictionary
history = {
    'fold_num': [],
    'epoch_num': [],
    'train_loss': [],
    'val_train_loss': [],
    'val_loss': []
}

# Initialize lists to store actual and predicted values
best_actual_values = []
best_predicted_values = []

# K-Fold Cross-Validation
for fold, (train_idx, val_idx) in enumerate(kf.split(train_data)):
    print('Fold {}'.format(fold + 1))
    
    # Separate train and validation data for this fold
    train_data_fold, val_data_fold = train_data[train_idx], train_data[val_idx]
    train_labels_fold, val_labels_fold = train_labels[train_idx], train_labels[val_idx]

    # DataLoader for training data
    train_dataset_fold = CustomDataset(train_data_fold, train_labels_fold)
    train_loader_fold = DataLoader(train_dataset_fold, batch_size=batch_size, shuffle=True)

    # DataLoader for validation data
    val_dataset_fold = CustomDataset(val_data_fold, val_labels_fold)
    val_loader_fold = DataLoader(val_dataset_fold, batch_size=batch_size, shuffle=False)

    # Initialize the model, optimizer, and scheduler for each fold
    model = Net(features).to(device)
    optimizer = optim.Adam(model.parameters(), lr=0.01)
    lambda1 = lambda epoch: epoch / 1000
    # scheduler = lr_scheduler.LambdaLR(optimizer, lambda1)
    scheduler = lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', patience=5, factor=0.5, verbose=True)
    
    # Train the model
    for epoch in range(epochs):
        print(f"Epoch {epoch+1} of {epochs}")

        # Train the model and get the training loss
        train_epoch_loss = fit(model, train_loader_fold, optimizer, criterion, device)/ len(train_loader_fold)
        print(f"Train Loss: {train_epoch_loss:.4f}")
        
        # Validate the model on training data (within fold) and get the loss
        val_train_epoch_loss = validate(model, train_loader_fold, criterion, device)/ len(train_loader_fold)
        print(f"Val Train Loss: {val_train_epoch_loss:.4f}")

        # Validate the model on validation data (within fold) and get the loss
        val_epoch_loss = validate(model, val_loader_fold, criterion, device)/ len(val_loader_fold)
        print(f"Val Loss: {val_epoch_loss:.4f}")

        # Check if current validation loss is better than the best validation loss
        if val_epoch_loss < best_val_loss:
            best_val_loss = val_epoch_loss
            best_model_weights = copy.deepcopy(model.state_dict())
            # Evaluate the best model on the validation set of this fold
            model.eval()
            with torch.no_grad():
                best_actual_values = []
                best_predicted_values = []
                for data, labels in tqdm(val_loader_fold, total=len(val_loader_fold)):
                    data = data.to(device)
                    labels = labels.to(device)
                    output = model(data)
                    best_actual_values.extend(labels.cpu().numpy())
                    best_predicted_values.extend(output.cpu().numpy())
        
        # Append fold and epoch numbers to history
        history['fold_num'].append(fold + 1)
        history['epoch_num'].append(epoch + 1)

        # Append training and validation losses to history
        history['train_loss'].append(train_epoch_loss)
        history['val_train_loss'].append(val_train_epoch_loss)
        history['val_loss'].append(val_epoch_loss)

        # Step the scheduler
        scheduler.step(val_epoch_loss)

# Save the actual and predicted values for the best model
np.savetxt("best_actual_values.csv", np.array(best_actual_values), delimiter=",")
np.savetxt("best_predicted_values.csv", np.array(best_predicted_values), delimiter=",")

# Calculate Mean Squared Error (MSE)
mse = mean_squared_error(best_actual_values, best_predicted_values)

# Calculate R-squared (R2)
r2 = r2_score(best_actual_values, best_predicted_values)

# Save MSE and R2
np.savetxt(path + "mse.csv", np.array([mse]), delimiter=",", header="MSE")
np.savetxt(path + "r2.csv", np.array([r2]), delimiter=",", header="R2")

# Plot actual vs predicted values
plt.figure(figsize=(8, 6))
plt.scatter(best_actual_values, best_predicted_values, color='blue', alpha=0.5)
plt.xlabel('Actual Values')
plt.ylabel('Predicted Values')
plt.title('Actual vs Predicted Values')
plt.grid(True)
plt.savefig(prediction_path)

# Load the best model weights
model.load_state_dict(best_model_weights)

# Save the model
torch.save(model.state_dict(), model_path)

# Save the training history
history_df = pd.DataFrame(history)
history_df.to_csv(history_path, index=False)
# Load the best model weights
model.load_state_dict(best_model_weights)

# Save the model
torch.save(model.state_dict(), model_path)

# Save the training history
history_df = pd.DataFrame(history)
history_df.to_csv(history_path, index=False)


# Report train losses
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams["figure.autolayout"] = True

fig, ((ax1, ax2)) = plt.subplots(2,1)

# history_df = pd.DataFrame(history)
history_df.set_index('epoch_num', inplace=True)
history_df.groupby('fold_num')['train_loss'].plot(ax=ax1, title ='MSE (train)')
history_df.groupby('fold_num')['val_train_loss'].plot(ax=ax2, title ='MSE (val-train)')

plt.legend()
plt.savefig(train_path)


# Report evaluation losses
plt.rcParams["figure.figsize"] = [10, 8]
plt.rcParams["figure.autolayout"] = True

fig, ax = plt.subplots(figsize=(8, 6))

history_df = pd.DataFrame(history)
history_df.set_index('epoch_num', inplace=True)
history_df.groupby('fold_num')['val_loss'].plot(ax=ax, title ='MSE (validation)')

plt.legend()
plt.savefig(val_path)


# Report evaluation metric
data_1 = history_df[history_df.fold_num == 1]
data_1_min = np.min(data_1['val_loss'])
data_2 = history_df[history_df.fold_num == 2]
data_2_min = np.min(data_2['val_loss'])
data_3 = history_df[history_df.fold_num == 3]
data_3_min = np.min(data_3['val_loss'])
data_4 = history_df[history_df.fold_num == 4]
data_4_min = np.min(data_4['val_loss'])
data_5 = history_df[history_df.fold_num == 5]
data_5_min = np.min(data_5['val_loss'])

MSE_mean = np.mean((data_1_min,data_2_min,data_3_min,data_4_min,data_5_min))
MSE_sd = np.std((data_1_min,data_2_min,data_3_min,data_4_min,data_5_min)) * 1.96

print(MSE_mean)
print(MSE_sd)

MSE_df = pd.DataFrame([MSE_mean,MSE_sd])
MSE_df.to_csv(MSE_path + '_df',
                  sep = ',',header=True,index=True)
np.savetxt(MSE_path, MSE_df, delimiter=',')

MSE_history = pd.DataFrame([(data_1_min,data_2_min,data_3_min,data_4_min,data_5_min)])
MSE_history.to_csv(MSE_path + '_history_df',
                  sep = ',',header=True,index=True)

np.savetxt(MSE_path, MSE_history, delimiter=',')

