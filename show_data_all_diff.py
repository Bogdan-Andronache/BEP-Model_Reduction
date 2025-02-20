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
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # Initialize images for both heatmaps
    cax1 = ax1.imshow(np.zeros((10, 10)), cmap="coolwarm", origin="lower", aspect="equal")
    cbar1 = fig.colorbar(cax1, ax=ax1, label="Value")
    ax1.set_xticks([]); ax1.set_yticks([])  # Hide ticks on both axes
    
    cax2 = ax2.imshow(np.zeros((10, 10)), cmap="gray", origin="lower", aspect="equal")  # Use grayscale for difference
    cbar2 = fig.colorbar(cax2, ax=ax2, label="Difference")
    ax2.set_xticks([]); ax2.set_yticks([])

    prev_data = None  # To store the previous data for calculating the difference
    
    for txt_file in txt_files:
        data = np.loadtxt(txt_file)
        size = int(np.sqrt(data.shape[0]))
        data = data.reshape((size, size))
        data = np.flipud(data)
        data = np.rot90(data, k=1)
        
        # Color scale for the data heatmap
        vmin, vmax = np.min(data), np.max(data)  # Adjust color scale dynamically
        cax1.set_norm(Normalize(vmin=vmin, vmax=vmax))
        cax1.set_data(data)
        ax1.set_title(f"Heatmap: {os.path.basename(txt_file)}")
        
        if prev_data is not None:
            # Calculate the difference from the previous frame
            diff = np.abs(data - prev_data)
            
            # Color scale for the difference heatmap (black-and-white)
            diff_vmin, diff_vmax = np.min(diff), np.max(diff)
            cax2.set_norm(Normalize(vmin=diff_vmin, vmax=diff_vmax))  # Normalize based on diff min/max
            cax2.set_data(diff)
        
        prev_data = data  # Update the previous data
        
        # Extract the solution number from the filename, e.g., 'sol1.00.txt' → 1.00
        solution_num = float(os.path.basename(txt_file)[3:-4])  # Remove 'sol' and '.txt'
        
        # Adjust the pause time to make the transition slower from 0.1 to 15
        if solution_num <= 15:
            plt.pause(0.1)  # Longer pause for solution numbers between 0.1 and 15
        else:
            plt.pause(0.01)  # Standard pause for solution numbers greater than 15
    
    plt.ioff()  # Turn off interactive mode
    plt.show()

def main():
    plot_live_data()

if __name__ == "__main__":
    main()
