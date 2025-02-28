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
HWND = 396972  # Change this to the actual HWND of your terminal

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
        send_keystrokes(f"poisson_nonlin1 <inputtxt/{input_file}")  # Run the command
        send_keystrokes("\n")

        time.sleep(3)

if __name__ == "__main__":
   start_time = time.time()

   input_files = [f for f in os.listdir('non-lin_poisson/non-lin_poisson/poisson/inputtxt') if os.path.isfile(os.path.join('non-lin_poisson/non-lin_poisson/poisson/inputtxt', f))]
   
   for input_file in tqdm(input_files, desc="Processing files", unit="file"):
       run_fortran_simulation(input_file)

   total_time = time.time() - start_time
   print(f"\nAll files processed in {total_time:.2f} seconds.")
