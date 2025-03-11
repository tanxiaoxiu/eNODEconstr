import torch
import numpy as np
import os
import h5py
from tqdm import tqdm
import torch.nn as nn
import torch.optim as optim
from torchdiffeq import odeint
from torch.utils.data import Dataset, DataLoader
import time
import random
import pandas as pd
from copy import deepcopy


N = 10
M = 10
learning_rate = 0.001
num_epochs = 1000
device = torch.device("cpu")

t5 = torch.tensor([0, 2, 6, 14, 24], dtype=torch.float32)
t10 = torch.tensor([0, 1, 2, 6, 9, 12, 15, 18, 21, 24], dtype=torch.float32)
t15 =  torch.tensor([0, 1, 2, 3, 4, 5, 7, 9, 11, 13, 15, 17, 19, 21, 24], dtype=torch.float32)
t20 = torch.tensor([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24], dtype=torch.float32)
t25 = torch.tensor([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], dtype=torch.float32)

class CustomDataset(Dataset):
    def __init__(self, tensor):
        self.tensor = tensor
    def __len__(self):
        return len(self.tensor)
    def __getitem__(self, idx):
        input_data = self.tensor[idx, 0, :]
        output_data = self.tensor[idx, :, :]
        return input_data, output_data

def loss_fn(pred_y, y):
    return torch.mean(torch.sum(torch.square(y - pred_y)))

class ODEFunc(torch.nn.Module):
    def __init__(self, N, M, mu_base):
        super(ODEFunc, self).__init__()

        self.N = N 
        self.M = M
        self.mu = torch.nn.Parameter(torch.rand(N, M, device=device))
        self.mu_base = mu_base
        self.mu_mask = (mu_base != 0).float()
        # Other parameters
        self.lambda_ = torch.nn.Parameter(torch.rand(1, device=device))
        self.m = torch.nn.Parameter(torch.rand(N, device=device))
        self.rho = torch.nn.Parameter(torch.rand(M, device=device))
        self.omega = torch.nn.Parameter(torch.rand(M, device=device))
        # Byproduct matrix
        nonzero_byproduct=0.5
        self.l = torch.nn.Parameter(self.generate_byproduct_matrix(M, nonzero_byproduct, device))
    
    @staticmethod
    def generate_byproduct_matrix(M, nonzero_byproduct, device):
        num_elements = M * M
        num_nonzeros = int(num_elements * nonzero_byproduct)
        values = np.concatenate([np.random.uniform(0, 1, num_nonzeros), np.zeros(num_elements - num_nonzeros)])
        np.random.shuffle(values)
        D = torch.tensor(values.reshape(M, M), dtype=torch.float32, device=device)
        col_sums = torch.sum(D, axis=0)
        col_sums[col_sums == 0] = 1
        D_norm = D / col_sums
        return D_norm

    def forward(self, t, ys):
        dydts_list = []
        for y in ys: 
            N, M = self.N, self.M
            C, R = y[:N], y[N:]
            C = C.unsqueeze(0)
            R = R.unsqueeze(0)

            dCdt = C * (R @ self.mu.t() * (1 - self.lambda_)) - C * self.m
            dRdt = self.rho - R * self.omega - (C @ self.mu) * R + self.lambda_ * ((C @ self.mu) * R) @ self.l.t()
            dydt = torch.cat([dCdt.squeeze(0), dRdt.squeeze(0)])
            dydt = dydt * (y > 0).float()
            dydts_list.append(dydt)
        dydts = torch.stack(dydts_list)
        return dydts

    def constrained(self):
        with torch.no_grad():
            self.lambda_.data = torch.relu(self.lambda_.data)
            self.m.data = torch.relu(self.m.data)
            self.rho.data = torch.relu(self.rho.data)
            self.omega.data = torch.relu(self.omega.data)
            self.l.data = torch.relu(self.l.data)

def set_seed(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True  
    torch.backends.cudnn.benchmark = False 

max_seed_attempts = 10 
factor = 1

for s in ["s10", "s15", "s20"]:
    for t in ["t5", "t10", "t15", "t20", "t25"]:
        for iter_num in range(1, 11):  
            current_seed = 1  
            successful_train = False  
            while not successful_train and current_seed < max_seed_attempts:
                try:
                    set_seed(current_seed)  
                    if t == "t5":
                        batch_t = t5.to(device)
                    elif t == "t10":
                        batch_t = t10.to(device)
                    elif t == "t15":
                        batch_t = t15.to(device)
                    elif t == "t20":
                        batch_t = t20.to(device)
                    elif t == "t25":
                        batch_t = t25.to(device)
            
                    input_file_path = f'~/eNODEconstr/simulation/generation_data/n{N}m{M}/{s}/{t}/iter{iter_num}'
                    output_dir_path = f'~/eNODEconstr/simulation/eNODEconstr/n{N}m{M}/training/{s}/{t}/iter{iter_num}'

                    mu_file = os.path.join(input_file_path, 'mu_binary_matrix.csv')
                    mu_base = pd.read_csv(mu_file, header=None, index_col=False)
                    mu_base = torch.tensor(mu_base.values, dtype=torch.float32, device=device)
                    os.makedirs(output_dir_path, exist_ok=True)
                    data_file = os.path.join(input_file_path, 'true_initial_state.h5')
                    with h5py.File(data_file, 'r') as hdf:
                        data = hdf['true_initial_state'][:]
                        tensor_data = torch.tensor(data, dtype=torch.float32).permute(0, 2, 1)
                    dataset = CustomDataset(tensor_data)
                    batch_s = int(int(s[1:])* 0.2)
                    dataloader = DataLoader(dataset, batch_size= batch_s, shuffle=False)
                    model = ODEFunc(N, M, mu_base).to(device)
                    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
                    best_loss = float('inf')
                    best_model = None
                    loss_history = []
                    avg_loss_file_path = os.path.join(output_dir_path, 'average_loss.txt')
                    with open(avg_loss_file_path, 'w') as file:
                        file.write("")

                    start_time = time.time()
                    for epoch in range(num_epochs):
                        for batch_idx, (input_data, output_data) in enumerate(dataloader):
                            input_data, output_data = input_data.to(device), output_data.to(device)
                            optimizer.zero_grad()
                            input_data = input_data.float()
                            pred_y = odeint(model, input_data, batch_t)
                            pred_y = pred_y.transpose(1, 0)
                            loss = loss_fn(pred_y, output_data)
                            l1_loss = torch.sum(torch.abs(torch.mul(model.mu, (1-model.mu_mask))))   
                            try:
                                with torch.no_grad():
                                    l1_penalty = factor * loss.item() / l1_loss.item()    
                            except:
                                l1_penalty = 0
                            loss += l1_penalty * l1_loss
                            loss.backward()
                            optimizer.step()
                            model.constrained()                          
                            if loss.item() < best_loss:
                                best_loss = deepcopy(loss.item())
                                best_model = deepcopy(model.state_dict())
                            loss_history.append(loss.item())

                        avg_loss = np.mean(loss_history[-len(dataloader):])  
                        print(f"Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}")
                        with open(avg_loss_file_path, 'a') as file:
                            file.write(f"Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}\n")
                   
                    end_time = time.time()
                    total_time = end_time - start_time
                    torch.save(best_model, os.path.join(output_dir_path, 'best_model.pth'))
                    np.savetxt(os.path.join(output_dir_path, 'loss_history.txt'), loss_history)
                    runtime_file_path = os.path.join(output_dir_path, 'runtime.txt')
                    with open(runtime_file_path, 'w') as file:
                        file.write(f"Model: {s}/{t}/iter{iter_num}. Total runtime: {total_time:.2f} seconds\n")                  
                    print(f"Model trained for {s}/{t}/iter{iter_num}. Total runtime: {total_time:.2f} seconds")
                    successful_train = True
                    print(f"Iter {iter_num} trained successfully with seed {current_seed}.")                 
                except AssertionError as e:
                    if "underflow in dt" in str(e):
                        print(f"Underflow error on Iter {iter_num} with seed {current_seed}. Trying a different seed...")
                        current_seed += 1  
                    else:
                        raise             
            if not successful_train:
                print(f"Failed to train Iter {iter_num} after trying {max_seed_attempts} seeds.")
                break 

