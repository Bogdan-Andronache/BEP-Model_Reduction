import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def plot_live_data():
    output_folder = "Data/Txt"
    
    values = [0, 14.0402, 16.4948, 17.2618, 18.8727, 19.3329]
    nr = values[1]
    #nr = 4.0683
    # Load the 'new_sol.npy' file and reshape it to (41, 41)
    pred_path = f"pred.npy"
    new_sol = np.load(pred_path).reshape((41, 41))
    
    sol_path = f"sol{nr}.txt"
    # Load the 'sol.00.txt' file for comparison
    sol_file = os.path.join(output_folder, sol_path)
    if not os.path.exists(sol_file):
        print(f"File {sol_file} does not exist.")
        return
    sol_data = np.loadtxt(sol_file)
    size = int(np.sqrt(sol_data.shape[0]))
    sol_data = sol_data.reshape((size, size))
    
    # Flip and rotate new_sol to match the orientation of sol_data
    new_sol = np.flipud(new_sol)
    new_sol = np.rot90(new_sol, k=1)

    plt.ion()  # Turn on interactive mode
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))  # 3 subplots: new_sol, sol, and difference
    
    cax1 = ax1.imshow(new_sol, cmap="coolwarm", origin="lower", aspect="equal")
    cbar1 = fig.colorbar(cax1, ax=ax1, label="Value")
    ax1.set_xticks([]); ax1.set_yticks([])  # Hide ticks on both axes
    ax1.set_title(f"Heatmap: {pred_path}")
    
    # Middle heatmap for 'sol.00.txt'
    cax2 = ax2.imshow(sol_data, cmap="coolwarm", origin="lower", aspect="equal")
    cbar2 = fig.colorbar(cax2, ax=ax2, label="Value")
    ax2.set_xticks([]); ax2.set_yticks([])  # Hide ticks on both axes
    ax2.set_title(f"Heatmap: {sol_path}")
    
    # Right heatmap for the difference between new_sol and sol.00.txt
    diff = np.abs(new_sol - sol_data)
    cax3 = ax3.imshow(diff, cmap="gray", origin="lower", aspect="equal")
    cbar3 = fig.colorbar(cax3, ax=ax3, label="Difference")
    ax3.set_xticks([]); ax3.set_yticks([])  # Hide ticks on both axes
    ax3.set_title(f"Difference from {sol_path}")
    
    # Dynamically adjust color scale for the difference heatmap
    diff_vmin, diff_vmax = np.min(diff), np.max(diff)
    cax3.set_norm(Normalize(vmin=diff_vmin, vmax=diff_vmax))  # Normalize based on diff min/max
    
    plt.ioff()  # Turn off interactive mode
    plt.show()

def main():
    plot_live_data()

if __name__ == "__main__":
    main()
