import platform
import os
import sys
import struct
import subprocess

def print_arch():
    print("[DEV] platform.machine():", platform.machine())
    print("[DEV] platform.processor():", platform.processor())
    print("[DEV] platform.platform():", platform.platform())
    print("[DEV] platform.architecture():", platform.architecture())
    print("[DEV] Python:", sys.version.split()[0])
    print("[DEV] Pointer bits:", struct.calcsize("P") * 8)

    try:
        print("os.uname():", os.uname())
    except AttributeError:
        pass

    for cmd in (["uname", "-m"], ["uname", "-a"]):
        try:
            out = subprocess.check_output(cmd, text=True).strip()
            print(" ".join(cmd) + ":", out)
        except Exception:
            pass

