import numpy as np
import os
import re

def extract_number(filename):
    match = re.search(r"sol\s*(\d+\.\d{4})\.txt", filename)  # Match "solX.XXXX.txt"
    return float(match.group(1)) if match else float('inf')

def make_snapshots():
    output_folder = "Data/Txt"
    
    txt_files = [f for f in os.listdir(output_folder) if os.path.isfile(os.path.join(output_folder, f))]
    txt_files = sorted(txt_files, key=extract_number)
    
    data_list = []
    for txt_file in txt_files:
        txt_file = output_folder + "/" + txt_file
        data = np.loadtxt(txt_file)
        data_list.append(data)
    
    snapshots = np.array(data_list)
    print(f"Combined data shape: {snapshots.shape}")
    return snapshots

def main():
    snapshots = make_snapshots()
    np.save("snapshots.npy", snapshots)
    print("Data saved to snapshots.npy")

if __name__ == "__main__":
    main()