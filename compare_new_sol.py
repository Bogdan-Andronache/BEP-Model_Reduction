import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def plot_live_data():
    output_folder = "Data/Txt"
    
    # Load the 'new_sol.npy' file and reshape it to (41, 41)
    new_sol = np.load('new_sol.npy').reshape((41, 41))
    
    # Load the 'sol15.00.txt' file for comparison
    sol15_file = os.path.join(output_folder, "sol5.00.txt")
    if not os.path.exists(sol15_file):
        print(f"File {sol15_file} does not exist.")
        return
    sol15_data = np.loadtxt(sol15_file)
    size = int(np.sqrt(sol15_data.shape[0]))
    sol15_data = sol15_data.reshape((size, size))
    sol15_data = np.flipud(sol15_data)
    sol15_data = np.rot90(sol15_data, k=1)
    
    # Flip and rotate new_sol to match the orientation of sol15_data
    new_sol = np.flipud(new_sol)
    new_sol = np.rot90(new_sol, k=1)

    # Ensure the shapes match after reshaping
    if new_sol.shape != sol15_data.shape:
        print(f"Shape mismatch: new_sol shape {new_sol.shape} != sol15_data shape {sol15_data.shape}")
        return

    plt.ion()  # Turn on interactive mode
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(24, 8))  # 3 subplots: new_sol, sol15, and difference
    
    # Left heatmap for 'new_sol.npy'
    cax1 = ax1.imshow(new_sol, cmap="coolwarm", origin="lower", aspect="equal")
    cbar1 = fig.colorbar(cax1, ax=ax1, label="Value")
    ax1.set_xticks([]); ax1.set_yticks([])  # Hide ticks on both axes
    ax1.set_title("Heatmap: new_sol.npy")
    
    # Middle heatmap for 'sol15.00.txt'
    cax2 = ax2.imshow(sol15_data, cmap="coolwarm", origin="lower", aspect="equal")
    cbar2 = fig.colorbar(cax2, ax=ax2, label="Value")
    ax2.set_xticks([]); ax2.set_yticks([])  # Hide ticks on both axes
    ax2.set_title("Heatmap: sol5.00.txt")
    
    # Right heatmap for the difference between new_sol and sol15.00.txt
    diff = np.abs(new_sol - sol15_data)
    cax3 = ax3.imshow(diff, cmap="gray", origin="lower", aspect="equal")
    cbar3 = fig.colorbar(cax3, ax=ax3, label="Difference")
    ax3.set_xticks([]); ax3.set_yticks([])  # Hide ticks on both axes
    ax3.set_title("Difference from sol5.00.txt")
    
    # Dynamically adjust color scale for the difference heatmap
    diff_vmin, diff_vmax = np.min(diff), np.max(diff)
    cax3.set_norm(Normalize(vmin=diff_vmin, vmax=diff_vmax))  # Normalize based on diff min/max
    
    plt.ioff()  # Turn off interactive mode
    plt.show()

def main():
    plot_live_data()

if __name__ == "__main__":
    main()
