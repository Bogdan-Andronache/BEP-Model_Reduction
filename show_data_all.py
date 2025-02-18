import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import time

def plot_live_data():
    output_folder = "Data/Txt"
    
    sol_numbers = [f"{i:.2f}" for i in np.arange(0.1, 100.01, 0.1)]  # Generate file names from 0.10 to 100.00
    txt_files = [os.path.join(output_folder, f"sol{num}.txt") for num in sol_numbers]
    txt_files = [f for f in txt_files if os.path.exists(f)]  # Filter existing files
    
    if not txt_files:
        print("No valid files to display.")
        return
    
    plt.ion()  # Turn on interactive mode
    fig, ax = plt.subplots(figsize=(8, 8))
    cax = ax.imshow(np.zeros((10, 10)), cmap="coolwarm", origin="lower", aspect="equal")
    cbar = fig.colorbar(cax, label="Value")
    ax.set_xticks([])
    ax.set_yticks([])
    
    for txt_file in txt_files:
        data = np.loadtxt(txt_file)
        size = int(np.sqrt(data.shape[0]))
        data = data.reshape((size, size))
        data = np.flipud(data)
        data = np.rot90(data, k=-1)
        
        vmin, vmax = np.min(data), np.max(data)  # Adjust color scale dynamically
        cax.set_norm(Normalize(vmin=vmin, vmax=vmax))
        cax.set_data(data)
        ax.set_title(f"Heatmap: {os.path.basename(txt_file)}")
        plt.pause(0.1)  # Pause to create a live animation
    
    plt.ioff()  # Turn off interactive mode
    plt.show()

def main():
    plot_live_data()

if __name__ == "__main__":
    main()