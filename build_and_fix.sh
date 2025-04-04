#!/bin/bash
# Script to build MMonitor and apply fixes

set -e  # Exit on error

echo "Starting build process for MMonitor..."

# Make the fix script executable
chmod +x fix_mac_app.sh

# Run the Python build script
python3 build_mac.py

# Run the fix script
./fix_mac_app.sh

echo "Build process completed successfully!"
echo "The app is located at: dist/MMonitor.app" 