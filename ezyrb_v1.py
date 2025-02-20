import numpy as np
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

from ezyrb import POD, RBF, Database
from ezyrb import ReducedOrderModel as ROM

snapshots = np.load("snapshots.npy")
values = np.linspace(0.1, 10, snapshots.shape[0])
param = np.array([[v] for v in values])
print('Snapshots max: ', np.max(snapshots), 'min: ', np.min(snapshots))

db = Database(param, snapshots)
pod = POD('svd', rank=3)
rbf = RBF()

rom = ROM(db, pod, rbf)
rom.fit()

new_alpha = [5.0]
pred_sol = rom.predict(new_alpha).snapshots_matrix
print('Prediction max: ', np.max(pred_sol), 'min: ', np.min(pred_sol))
np.save("new_sol.npy", pred_sol)