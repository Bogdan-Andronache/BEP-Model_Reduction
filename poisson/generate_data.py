import ctypes
import time
import win32gui
import win32con
import win32api
import re
import numpy as np
from tqdm import tqdm

# Set the window handle (HWND) of the terminal
HWND = 657978  # Change this to the actual HWND of your terminal

# Bring the window to the foreground
def focus_terminal():
    win32gui.ShowWindow(HWND, win32con.SW_RESTORE)  # Restore if minimized
    win32gui.SetForegroundWindow(HWND)

# Function to send keystrokes
def send_keystrokes(text):
    for char in text:
        if char == '<':  
            vk_code = 0xBC  # Virtual key code for '<' (same key as ',')
            shift_pressed = True
        else:
            vk_code = win32api.VkKeyScan(char) & 0xFF
            shift_pressed = (win32api.VkKeyScan(char) >> 8) & 1  

        if shift_pressed:
            ctypes.windll.user32.keybd_event(win32con.VK_SHIFT, 0, 0, 0)

        ctypes.windll.user32.keybd_event(vk_code, 0, 0, 0)
        time.sleep(0.02)
        ctypes.windll.user32.keybd_event(vk_code, 0, win32con.KEYEVENTF_KEYUP, 0)

        if shift_pressed:
            ctypes.windll.user32.keybd_event(win32con.VK_SHIFT, 0, win32con.KEYEVENTF_KEYUP, 0)

        time.sleep(0.001)

# Function to automate the process
def run_fortran_simulation(input_file):
        focus_terminal()  # Bring the terminal to the foreground
        send_keystrokes(f"poisson1 <{input_file}")  # Run the command
        send_keystrokes("\n")

        time.sleep(5)

def create_input_file(filename, alpha_value):
   # Create an input file with specified alpha value in namelist format
   with open(filename, 'w') as f:
       f.write("&comppar\n")
       f.write(f"alpha = {alpha_value:.4f}\n")
       f.write("/\n")

def append_output_files(input_files, output_combined):
   # Open the combined output file for writing
   with open(output_combined, 'w') as f_combined:
       # Iterate over each input file
       for input_file in input_files:
           # Read the contents of the current output file and append to combined file
           with open(input_file.replace("input_", "output_"), 'r') as f_output:
               f_combined.write(f_output.read())

def write_alpha_values(alpha_values, filename):
   # Write the alpha values to a specified file
   with open(filename, 'w') as f:
       for alpha in alpha_values:
           f.write(f"{alpha:.4f}\n")

# Run the automation

alpha_values = np.linspace(1.0, 20.1, 250)

# Write the alpha values to a new file
write_alpha_values(alpha_values, "parameters.txt")
print("Alpha values have been saved to parameters.txt.")

if __name__ == "__main__":
   # Create input files with different alpha values
   input_files = []
   for alpha in alpha_values:
       filename = f"input_{alpha:.4f}.txt"
       create_input_file(filename, alpha)
       input_files.append(filename)

   start_time = time.time()

   # Run simulations for each input file
   for input_file in tqdm(input_files, desc="Processing files", unit="file"):
       run_fortran_simulation(input_file)

   total_time = time.time() - start_time
   print(f"\nAll files processed in {total_time:.2f} seconds.")
