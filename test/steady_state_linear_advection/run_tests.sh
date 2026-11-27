#!/bin/bash

# 1. Immediately move to the directory where this script lives
cd "$(dirname "$0")"

# Define the root of your project and the build folder path
PROJECT_ROOT="../.."
BUILD_DIR="$PROJECT_ROOT/build"

# 2. Clean up any previous partial builds
#rm -rf "$BUILD_DIR"

# 3. Configure and build CMake by explicitly targeting the project root
cmake -S "$PROJECT_ROOT" -B "$BUILD_DIR"
cmake --build "$BUILD_DIR"

# 4. Reset the CSV result file in the local script directory
rm -f l2_norm_vs_h.csv
echo "K,N,l2_norm" > l2_norm_vs_h.csv

# 5. Point to the newly generated executable binary inside the build directory
BINARY="$BUILD_DIR/linear_advection_steady_state"

# 6. Run the sequential testing loops
for p in {0..7}; do
  K=$((2**p))
  
  for N in {0..11}; do
    echo "Processing K=$K (2^$p), N=$N"
    # This blocks and waits until the binary finishes running
    echo "i, x, uh, u0" > "K${K}_N${N}.csv"
    "$BINARY" "$K" "$N"
  done
done

# 7. Generate output plots sequentially without overwriting logs
python3 plot_csv.py > run_output.txt
python3 plot_convergence.py >> run_output.txt
mv l2_norm_vs_deltax.svg ./svg
echo "All runs completed."
