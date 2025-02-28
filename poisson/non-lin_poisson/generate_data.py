# Python script to generate text files with varying delta and gamma
import ctypes
import time
import win32gui
import win32con
import win32api
import re
import numpy as np
from tqdm import tqdm
import os

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
        send_keystrokes(f"poisson_nonlin1 <{input_file}")  # Run the command
        send_keystrokes("\n")

        time.sleep(2)

if __name__ == "__main__":
    delta_values = [i for i in range(1, 11)]  # Delta values from 1 to 10
    gamma_values = [round(0.1 * i, 1) for i in range(1, 11)]  # Gamma values from 0.1 to 1.0

    input_files = []
    # Create 100 txt files with the given pattern
    file_counter = 0
    for delta in delta_values:
        for gamma in gamma_values:
            file_counter += 1
            # Prepare file name (you can modify this as needed)
            file_name = f"input_{delta}_{gamma}.txt"
            input_files.append(file_name)
            
            # Template content with delta and gamma substituted
            content = f"""
&comppar
outputfile = 'sol_{delta}_{gamma}.vtk'
delta = {delta:.1f}
gamma = {gamma:.1f}
nx = 20
ny = 20
eps = 1e-9
nPicard = 10,
itermax = 100
/
"""
            # Write to the file
            with open(file_name, 'w') as file:
                file.write(content)


    start_time = time.time()
    print(input_files)

   # Run simulations for each input file
    for input_file in tqdm(input_files, desc="Processing files", unit="file"):
       run_fortran_simulation(input_file)

    total_time = time.time() - start_time
    print(f"\nAll files processed in {total_time:.2f} seconds.")
