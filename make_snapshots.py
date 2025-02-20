import numpy as np
import os

def make_snapshots(min, max):
    output_folder = "Data/Txt"
    
    sol_numbers = []
    for i in np.arange(min, max + 0.1, 0.1):
        sol_numbers.append(f"sol{i:.2f}.txt")
    
    txt_files = [os.path.join(output_folder, num) for num in sol_numbers]
    txt_files = [f for f in txt_files if os.path.exists(f)]
    
    data_list = []
    for txt_file in txt_files:
        data = np.loadtxt(txt_file)
        data_list.append(data)
    
    snapshots = np.array(data_list)
    print(f"Combined data shape: {snapshots.shape}")
    return snapshots

def main():
    snapshots = make_snapshots(0.1, 10.0)
    np.save("snapshots.npy", snapshots)
    print("Data saved to snapshots.npy")

if __name__ == "__main__":
    main()