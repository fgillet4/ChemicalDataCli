#!/bin/bash

# ChemicalDataCli Auto-Startup Script
# This script automatically activates the virtual environment and runs the application

echo "🧪 Starting Chemical Data CLI..."
echo "📁 Current directory: $(pwd)"

# Check if virtual environment exists
if [ ! -d "chemdata_env" ]; then
    echo "❌ Virtual environment 'chemdata_env' not found!"
    echo "Please create it first with: python3 -m venv chemdata_env"
    exit 1
fi

# Check if main.py exists
if [ ! -f "main.py" ]; then
    echo "❌ main.py not found!"
    echo "Make sure you're in the correct directory"
    exit 1
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source chemdata_env/bin/activate

# Check if activation was successful
if [ -z "$VIRTUAL_ENV" ]; then
    echo "❌ Failed to activate virtual environment!"
    exit 1
fi

echo "✅ Virtual environment activated: $VIRTUAL_ENV"

# Check if required packages are installed
echo "🔍 Checking required packages..."
python -c "import chemicals" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "⚠️  'chemicals' package not found. Installing requirements..."
    if [ -f "requirements.txt" ]; then
        pip install -r requirements.txt
    else
        echo "Installing basic requirements..."
        pip install chemicals thermo fluids ht
    fi
fi

# Run the application
echo "🚀 Launching Chemical Data CLI..."
echo "=========================================="
python main.py

# Deactivate virtual environment when done
echo "=========================================="
echo "👋 Deactivating virtual environment..."
deactivate
echo "✅ Chemical Data CLI session ended."