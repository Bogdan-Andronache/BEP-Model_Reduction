import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import sys

def plot_data(txt_files):
    fig, axes = plt.subplots(1, len(txt_files), figsize=(8 * len(txt_files), 8))
    if len(txt_files) == 1:
        axes = [axes]  # Ensure axes is iterable
    
    for ax, txt_file in zip(axes, txt_files):
        data = np.loadtxt(txt_file)
        size = int(np.sqrt(data.shape[0]))
        data = data.reshape((size, size))
        data = np.flipud(data)
        data = np.rot90(data, k=1)

        # Set vmin and vmax based on the data's min and max values
        vmin, vmax = data.min(), data.max()
        norm = Normalize(vmin=vmin, vmax=vmax)

        cax = ax.imshow(data, cmap="coolwarm", origin="lower", aspect="equal", norm=norm)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"Heatmap: {os.path.basename(txt_file)}")
    
    fig.colorbar(cax, ax=axes, fraction=0.046, pad=0.04, label="Value")
    plt.show()

def main():
    output_folder = "Data/Txt"
    
    if len(sys.argv) > 1:
        txt_files = []
        for arg in sys.argv[1:]:
            sol_number = float(arg)
            formatted_sol_number = f"{sol_number:.2f}"
            txt_file = os.path.join(output_folder, f"sol{formatted_sol_number}.txt")
            if os.path.exists(txt_file):
                txt_files.append(txt_file)
            else:
                print(f"File not found: {txt_file}")
        
        if txt_files:
            plot_data(txt_files)
        else:
            print("No valid files to display.")
    else:
        print("Please provide at least one solution number as an argument.")

if __name__ == "__main__":
    main()
