#!/bin/bash

cd "$(dirname "$0")"

BINARY="../../../build/module_test_L2_projection"

# Loop p from 0 to 10 (K = 2^0 to 2^10)
for p in {0..7}; do
  K=$((2**p))
  
  for N in {0..11}; do
    echo "Processing K=$K (2^$p), N=$N"
    "$BINARY" "$K" "$N" > "K${K}_N${N}.csv"
  done
done

echo "All runs completed."
