#!/bin/bash

# Exit immediately if a command fails
set -e

echo "=== S_jq Matrix Batch Plotter ==="
python3 plot_csv.py
echo "=== Done! All plots saved to ./plots ==="
