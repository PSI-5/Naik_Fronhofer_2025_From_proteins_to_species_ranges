#!/bin/bash

# install_dependencies.sh
# Install all Python dependencies for the thermal adaptation models project

echo "Installing dependencies for thermal adaptation models project..."

# Update pip
python3 -m pip install --upgrade pip

# Install core scientific computing packages
pip install numpy pandas scipy matplotlib seaborn

# Install specialized packages
pip install alive-progress  # For progress bars
pip install regex           # For advanced regular expressions

# Install optional packages (uncomment if needed)
# pip install jupyter        # For interactive notebooks

echo "Dependencies installed successfully!"
