import numpy as np
import matplotlib.tri as mtri
import matplotlib.pyplot as plt

from ezyrb import POD, RBF, GPR, Database
from ezyrb import ReducedOrderModel as ROM

snapshots = np.load("snapshots.npy")

delta_values = np.linspace(1, 10, 100)
gamma_values = np.linspace(0.1, 1, 100)

# Stack them together into a (100, 2) array
result = np.column_stack((delta_values, gamma_values))
print(result.shape)

# db = Database(param, snapshots)
# pod = POD('svd')

# rom = ROM(db, pod, GPR())
# rom.fit()

# pred_sol = rom.predict([2.3, 0.45]).snapshots_matrix
# np.save(f"pred{new_alpha}.npy", pred_sol)