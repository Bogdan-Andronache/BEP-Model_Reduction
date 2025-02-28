import numpy as np
import matplotlib.pyplot as plt
from ezyrb import POD, RBF, Database
from ezyrb import ReducedOrderModel as ROM

# Assuming snapshots.npy and parameters.npy are already loaded
snapshots = np.load("snapshots.npy")
param = np.load("parameters.npy")

# Shuffle data (same as your original code)
indices = np.random.permutation(param.shape[0])
snapshots = snapshots[indices]
param = param[indices]

testing_ratio = 0.30
split_idx = int((1 - testing_ratio) * param.shape[0])
train_snap, testing_snap = snapshots[:split_idx], snapshots[split_idx:]
train_param, testing_param = param[:split_idx], param[split_idx:]

db = Database(train_param, train_snap)

pod = POD('svd', rank=10)
rom = ROM(db, pod, RBF())
rom.fit()

relative_errors_grid = np.zeros((testing_param.shape[0],))  # Store relative errors for each test

for i, (test_param, true_snap) in enumerate(zip(testing_param, testing_snap)):
    # Get the predicted solution (predicted nodes)
    pred_sol = rom.predict(test_param).snapshots_matrix

    # Reshape the prediction and true solution to 1D
    pred_sol_flat = pred_sol.flatten()
    true_snap_flat = true_snap.flatten()

    # Compute the difference (prediction - true solution)
    diff = pred_sol_flat - true_snap_flat

    # Compute the relative error using the L2 norm
    rel_error = np.linalg.norm(diff) / np.linalg.norm(true_snap_flat)

    # Store the relative error in the grid
    relative_errors_grid[i] = rel_error

# Assuming your testing_param has two parameters, reshape it into a grid
x_param = testing_param[:, 0]
y_param = testing_param[:, 1]

# Create a scatter plot with color scaling
plt.figure(figsize=(8, 6))
plt.scatter(x_param, y_param, c=relative_errors_grid, cmap='coolwarm', s=50)
plt.colorbar(label="Relative Error")  # Color bar for the error scale
plt.title("Relative Error Based on Testing Parameters")
plt.xlabel("Delta")
plt.ylabel("Gamma")
plt.grid(True)
plt.show()
