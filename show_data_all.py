import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import re

def extract_number(filename):
    match = re.search(r"sol\s*(\d+\.\d{4})\.txt", filename)  # Match "solX.XXXX.txt"
    return float(match.group(1)) if match else float('inf')

def plot_and_save_gif():
    output_folder = "Data/Txt"
    
    txt_files = [f for f in os.listdir(output_folder) if os.path.isfile(os.path.join(output_folder, f))]
    txt_files = sorted(txt_files, key=extract_number)
    
    if not txt_files:
        print("No valid files to display.")
        return
    
    fig, ax = plt.subplots(figsize=(8, 8))
    cax = ax.imshow(np.zeros((10, 10)), cmap="coolwarm", origin="lower", aspect="equal")
    cbar = fig.colorbar(cax, label="Value")
    ax.set_xticks([])
    ax.set_yticks([])

    def update(frame):
        txt_file = os.path.join(output_folder, txt_files[frame])
        data = np.loadtxt(txt_file)
        size = int(np.sqrt(data.shape[0]))
        data = data.reshape((size, size))
        data = np.flipud(data)
        data = np.rot90(data, k=1)

        vmin, vmax = np.min(data), np.max(data)  # Adjust color scale dynamically
        cax.set_norm(Normalize(vmin=vmin, vmax=vmax))
        cax.set_data(data)
        ax.set_title(f"Heatmap: {os.path.basename(txt_file)}")
        return cax,

    # Create the animation
    ani = animation.FuncAnimation(fig, update, frames=len(txt_files), interval=100, blit=False)

    # Show animation
    plt.show()

    # Save animation as GIF
    ani.save("heatmap_animation.gif", writer='pillow', fps=5)

    print("GIF saved as heatmap_animation.gif")

    plt.close(fig)  # Close the figure after saving

def main():
    plot_and_save_gif()

if __name__ == "__main__":
    main()
