import numpy as np
import os
import re

def extract_number(filename):
    match = re.search(r"sol_(\d+\.\d+)_([\d\.]+)\.txt", filename)  # Match "sol_X.Y_Z.Z.txt"
    if match:
        # Extract both numbers as floats
        num1 = float(match.group(1))
        num2 = float(match.group(2))
        return (num1, num2)
    return (float('inf'), float('inf'))  # Return a large value if match fails

def make_snapshots():
    output_folder = "Data/Txt"
    
    txt_files = [f for f in os.listdir(output_folder) if os.path.isfile(os.path.join(output_folder, f))]
    txt_files = sorted(txt_files, key=extract_number)  # Sort based on both numbers
    
    # Display the sorted filenames
    print("Sorted filenames:")
    for filename in txt_files:
        print(filename)
    
    data_list = []
    for txt_file in txt_files:
        txt_file = output_folder + "/" + txt_file
        data = np.loadtxt(txt_file)
        data_list.append(data)
    
    snapshots = np.array(data_list)
    print(f"Combined data shape: {snapshots.shape}")
    return snapshots

def save_snapshots(snapshots):
    np.save("snapshots.npy", snapshots)
    print("Data saved to snapshots.npy")
    
    # Save to snapshots.txt
    np.savetxt("snapshots.txt", snapshots.reshape(-1, snapshots.shape[-1]), fmt='%.6f')  # Flatten if needed
    print("Data saved to snapshots.txt")

def main():
    snapshots = make_snapshots()
    save_snapshots(snapshots)

if __name__ == "__main__":
    main()
