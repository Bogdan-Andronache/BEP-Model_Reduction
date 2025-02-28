import numpy as np
import matplotlib.pyplot as plt
from ezyrb import POD, RBF, GPR, Database
from ezyrb import ReducedOrderModel as ROM
from matplotlib.tri import Triangulation

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
pod = POD('svd', rank=50)

rom = ROM(db, pod, RBF())
rom.fit()

modes = rom.reduction.modes
print(modes.shape)

plt.figure(figsize=(8, 6))
plt.plot(rom.reduction.singular_values, marker='o', linestyle='-', color='b')
plt.title("Singular values")
plt.xlabel("Index")
plt.ylabel("Singular values")
plt.yscale('log')
plt.grid(True, which="both", ls="--")
plt.show()

grid_size = 41
x = np.linspace(0, 1, grid_size)
y = np.linspace(0, 1, grid_size)
X, Y = np.meshgrid(x, y)
points = np.vstack([X.flatten(), Y.flatten()]).T

triang = Triangulation(points[:, 0], points[:, 1])

fig, axes = plt.subplots(5, 10, figsize=(20, 24))

axes = axes.flatten()

for i in range(50):
    mode = modes[:, i].reshape(grid_size, grid_size)
    axes[i].tripcolor(triang, mode.flatten(), cmap='coolwarm')
    axes[i].set_title(f'Mode {i+1}')
    axes[i].set_xlabel('X')
    axes[i].set_ylabel('Y')
    axes[i].set_aspect('equal')

plt.tight_layout()
plt.show()
