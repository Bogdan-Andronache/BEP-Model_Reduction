import ctypes
import time
import win32gui
import win32con
import win32api
import re

# Set the window handle (HWND) of the terminal
HWND = 657978  # Change this to the actual HWND of your terminal

# Bring the window to the foreground
def focus_terminal():
    win32gui.ShowWindow(HWND, win32con.SW_RESTORE)  # Restore if minimized
    win32gui.SetForegroundWindow(HWND)
    time.sleep(0.1)  # Allow time for the window to be ready

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

        time.sleep(0.02)

# Function to modify input.txt and update alpha
def update_input_file(alpha_value):
    with open("poisson/input.txt", "r") as file:
        lines = file.readlines()

    # Modify the line containing alpha
    new_lines = [re.sub(r'alpha\s*=\s*\d+\.\d+', f'alpha = {alpha_value:.1f}', line) for line in lines]

    # Write the new content back
    with open("poisson/input.txt", "w") as file:
        file.writelines(new_lines)

# Function to automate the process
def automate_poisson1():
    alpha = 0.1  # Start value

    while alpha <= 10.0:
        update_input_file(alpha)  # Update the input.txt file

        focus_terminal()  # Bring the terminal to the foreground
        send_keystrokes("poisson1 <input.txt")  # Run the command
        send_keystrokes("\n")  # Press Enter

        time.sleep(.2)  # Allow time for execution

        alpha += 0.1  # Increment alpha

# Run the automation
automate_poisson1()
