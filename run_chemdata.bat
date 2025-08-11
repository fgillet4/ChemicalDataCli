@echo off
REM ChemicalDataCli Auto-Startup Script for Windows
REM This script automatically activates the virtual environment and runs the application

echo 🧪 Starting Chemical Data CLI...
echo 📁 Current directory: %CD%

REM Check if virtual environment exists
if not exist "chemdata_env" (
    echo ❌ Virtual environment 'chemdata_env' not found!
    echo Please create it first with: python -m venv chemdata_env
    pause
    exit /b 1
)

REM Check if main.py exists
if not exist "main.py" (
    echo ❌ main.py not found!
    echo Make sure you're in the correct directory
    pause
    exit /b 1
)

REM Activate virtual environment
echo 🔧 Activating virtual environment...
call chemdata_env\Scripts\activate.bat

REM Check if required packages are installed
echo 🔍 Checking required packages...
python -c "import chemicals" 2>nul
if errorlevel 1 (
    echo ⚠️  'chemicals' package not found. Installing requirements...
    if exist "requirements.txt" (
        pip install -r requirements.txt
    ) else (
        echo Installing basic requirements...
        pip install chemicals thermo fluids ht
    )
)

REM Run the application
echo 🚀 Launching Chemical Data CLI...
echo ==========================================
python main.py

REM Deactivate virtual environment when done
echo ==========================================
echo 👋 Deactivating virtual environment...
call chemdata_env\Scripts\deactivate.bat
echo ✅ Chemical Data CLI session ended.
pause