import subprocess
import numpy as np
import pygetwindow as gw
import pyautogui
import time

def create_input_file(filename, alpha):
    # Create an input file with only the alpha parameter
    with open(filename, 'w') as f:
        f.write("&comppar\n")
        f.write(f"alpha = {alpha:.4f}\n")
        f.write("/\n")

def find_bash_window():
    # Find the bash window with HWND = 657978
    for window in gw.getAllWindows():
        if "657978" in window.title:
            return window
    return None

def run_fortran_simulation(input_file):
    # Send `poisson1 < input_file` command to the bash window
    bash_window = find_bash_window()
    if bash_window:
        bash_window.activate()
        time.sleep(1)  # Give time to switch windows
        pyautogui.write(f'poisson1 < {input_file}\n')
        pyautogui.press('enter')
    else:
        print("Bash window HWND = 657978 not found.")

def append_output_files(input_files, output_combined):
    # Combine output files into a single file
    with open(output_combined, 'w') as f_combined:
        for input_file in input_files:
            output_file = input_file.replace("input_", "output_")
            try:
                with open(output_file, 'r') as f_output:
                    f_combined.write(f_output.read())
            except FileNotFoundError:
                print(f"Warning: {output_file} not found. Skipping...")

def write_alpha_values(alpha_values, filename):
    # Save alpha values to a file
    with open(filename, 'w') as f:
        for alpha in alpha_values:
            f.write(f"{alpha:.4f}\n")

# Generate alpha values from 0.1 to 10.0 with step 0.1
alpha_values = np.arange(0.1, 10.1, 0.1)
write_alpha_values(alpha_values, "parameters.txt")
print("Alpha values have been saved to parameters.txt.")

if __name__ == "__main__":
    input_files = []
    for alpha in alpha_values:
        filename = f"input_{alpha:.4f}.txt"
        create_input_file(filename, alpha)
        input_files.append(filename)
    
    # Run simulations for each input file in HWND = 657978 bash
    for input_file in input_files:
        run_fortran_simulation(input_file)
    
    # Append output from all simulations into a single file
    output_combined = "snapshots.txt"
    append_output_files(input_files, output_combined)
    print("Snapshots created.")