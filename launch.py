#!/usr/bin/env python3
"""
ChemicalDataCli Launcher
========================

This script automatically handles virtual environment activation and 
launches the Chemical Data CLI application.

Usage:
    python launch.py
    
Or simply double-click this file if Python is properly configured.
"""

import os
import sys
import subprocess
import platform

def main():
    """Main launcher function."""
    print("🧪 Chemical Data CLI Launcher")
    print("=" * 40)
    
    # Get the directory where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    
    print(f"📁 Working directory: {script_dir}")
    
    # Check for virtual environment
    venv_path = os.path.join(script_dir, "chemdata_env")
    if not os.path.exists(venv_path):
        print("❌ Virtual environment 'chemdata_env' not found!")
        print("\nTo create it, run:")
        print("python -m venv chemdata_env")
        input("\nPress Enter to exit...")
        return
    
    # Check for main.py
    main_py = os.path.join(script_dir, "main.py")
    if not os.path.exists(main_py):
        print("❌ main.py not found!")
        input("\nPress Enter to exit...")
        return
    
    # Determine the correct Python executable in the virtual environment
    system = platform.system()
    if system == "Windows":
        python_exe = os.path.join(venv_path, "Scripts", "python.exe")
        pip_exe = os.path.join(venv_path, "Scripts", "pip.exe")
    else:  # Unix/Linux/Mac
        python_exe = os.path.join(venv_path, "bin", "python")
        pip_exe = os.path.join(venv_path, "bin", "pip")
    
    if not os.path.exists(python_exe):
        print(f"❌ Python executable not found at: {python_exe}")
        input("\nPress Enter to exit...")
        return
    
    print("✅ Virtual environment found")
    
    # Check if required packages are installed
    print("🔍 Checking required packages...")
    try:
        result = subprocess.run([python_exe, "-c", "import chemicals"], 
                              capture_output=True, text=True)
        if result.returncode != 0:
            print("⚠️  Required packages not found. Installing...")
            
            # Install requirements
            requirements_file = os.path.join(script_dir, "requirements.txt")
            if os.path.exists(requirements_file):
                subprocess.run([pip_exe, "install", "-r", requirements_file])
            else:
                print("Installing basic requirements...")
                subprocess.run([pip_exe, "install", "chemicals", "thermo", "fluids", "ht"])
            
            print("✅ Requirements installed")
        else:
            print("✅ Required packages found")
    
    except Exception as e:
        print(f"⚠️  Error checking packages: {e}")
        print("Proceeding anyway...")
    
    # Launch the application
    print("\n🚀 Launching Chemical Data CLI...")
    print("=" * 40)
    
    try:
        # Use the virtual environment's Python to run main.py
        subprocess.run([python_exe, main_py])
        
    except KeyboardInterrupt:
        print("\n\n👋 Application interrupted by user")
    except Exception as e:
        print(f"\n❌ Error running application: {e}")
    
    print("\n=" * 40)
    print("✅ Chemical Data CLI session ended")
    
    # Keep window open on Windows
    if system == "Windows":
        input("\nPress Enter to exit...")

if __name__ == "__main__":
    main()