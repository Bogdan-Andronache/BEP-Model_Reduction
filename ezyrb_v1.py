import numpy as np
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

from ezyrb import POD, RBF, GPR, Database
from ezyrb import ReducedOrderModel as ROM

snapshots = np.load("snapshots.npy")
values = np.loadtxt("parameters.txt")[:snapshots.shape[0]]

indices = np.random.permutation(snapshots.shape[0])
snapshots = snapshots[indices]
values = values[indices]

param = np.array([[v] for v in values])

testing_ratio = 0.2
split_idx = int((1 - testing_ratio) * 120)
train_snap, test_snap = snapshots[:split_idx], snapshots[split_idx:]
train_param, test_param = param[:split_idx], param[split_idx:]

db = Database(train_param, train_snap)
pod = POD('svd')

rom = ROM(db, pod, GPR())
rom.fit()

for new_alpha in (14.0402, 16.4948, 17.2618, 18.8727, 19.3329):
    pred_sol = rom.predict([new_alpha]).snapshots_matrix
    np.save(f"pred{new_alpha}.npy", pred_sol)

for new_alpha in test_param:
    pred_sol = rom.predict(new_alpha).snapshots_matrix
    np.save(f"pred{new_alpha[0]}.npy", pred_sol)