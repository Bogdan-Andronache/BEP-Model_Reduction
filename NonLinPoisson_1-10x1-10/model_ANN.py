import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn

from ezyrb import POD, RBF, GPR, Database, ANN
from ezyrb import ReducedOrderModel as ROM

from tqdm import tqdm
import time

snapshots = np.load("snapshots.npy")
param = np.load("parameters.npy")

print(param.shape, snapshots.shape)

indices = np.random.permutation(param.shape[0])
snapshots = snapshots[indices]
param = param[indices]

testing_ratio = 0.20
split_idx = int((1 - testing_ratio) * param.shape[0])
train_snap, testing_snap = snapshots[:split_idx], snapshots[split_idx:]
train_param, testing_param = param[:split_idx], param[split_idx:]

db = Database(train_param, train_snap)

average_errors = []
start_time = time.time()

for rank in tqdm(range(1, 21), desc="Processing files", unit="file"):
    pod = POD('svd', rank=rank)

    rom = ROM(db, pod, ANN([30, 20, 30], nn.Tanh(), 1e-5))
    rom.fit()

    relative_errors = []

    for i, (test_param, true_snap) in enumerate(zip(testing_param, testing_snap)):
        # Get the predicted solution (predicted nodes)
        pred_sol = rom.predict(test_param).snapshots_matrix
        
        pred_sol_flat = pred_sol.flatten()
        true_snap_flat = true_snap.flatten()

        diff = pred_sol_flat - true_snap_flat
        rel_error = np.linalg.norm(diff) / np.linalg.norm(true_snap_flat)

        relative_errors.append(rel_error)

    avg_error = np.mean(relative_errors)
    average_errors.append(avg_error)
    
total_time = time.time() - start_time
print(f"\nAll files processed in {total_time:.2f} seconds.")

plt.figure(figsize=(8, 6))
plt.plot(average_errors, marker='o', linestyle='-', color='b')
plt.title("Average Testing Error vs Rank - Param Space: 1-10x1-10")
plt.xlabel("Rank")
plt.ylabel("Average Relative Error")
plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.show()
