import win32gui

def list_all_windows():
    def callback(hwnd, windows):
        if win32gui.IsWindowVisible(hwnd):
            windows.append((hwnd, win32gui.GetWindowText(hwnd)))
        return True

    windows = []
    win32gui.EnumWindows(callback, windows)
    return windows

for hwnd, title in list_all_windows():
    print(f"HWND: {hwnd}, Title: {title}")