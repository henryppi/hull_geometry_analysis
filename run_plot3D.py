import sys
if sys.platform == "linux" or sys.platform == "linux2":
    print("Running on Linux")
    run_os = 'linux'
elif sys.platform == "darwin":
    print("Running on macOS")
    run_os = 'macos'

