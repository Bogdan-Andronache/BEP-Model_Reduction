import numpy as np
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

from ezyrb import POD, RBF, GPR, Database
from ezyrb import ReducedOrderModel as ROM

snapshots = np.load("poisson/non-lin_poisson/snapshots.npy")

delta_values = np.linspace(1, 10, 100)
gamma_values = np.linspace(0.1, 1, 100)

# Stack them together into a (100, 2) array
result = np.column_stack((delta_values, gamma_values))
print(snapshots.shape)

db = Database(result, snapshots)
pod = POD('svd', rank=20)

rom = ROM(db, pod, RBF())
rom.fit()

pred_sol = rom.predict([2.3, 0.45]).snapshots_matrix
np.save(f"pred.npy", pred_sol)