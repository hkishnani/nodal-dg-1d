import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# 1. Read data from the CSV file
# Expected CSV format: columns named K, N, and l2_norm
try:
    df = pd.read_csv('l2_norm_vs_h.csv')
except FileNotFoundError:
    print("Error: 'l2_norm_vs_h.csv' not found. Please ensure it is in the same directory.")
    exit()

# Domain details: Uniform grid over L = 2*PI
domain_length = 2 * np.pi

# 2. Setup Figure Layout
plt.figure(figsize=(9, 6), dpi=100)
# Clean, distinct presentation palette (Modern Corporate style)
colors = ['#1d3557', '#457b9d', '#2a9d8f', '#76c893', '#e9c46a', '#e76f51', '#d62828']

# 3. Process and Plot per N (0 to 6)
for idx, N in enumerate(range(7)):
    # Filter data for the current N and sort by K ascending (which means dx descending)
    n_data = df[df['N'] == N].sort_values(by='K')
    
    if n_data.empty:
        continue
    
    # Calculate delta_x = 2*pi / K
    K_vals = n_data['K'].to_numpy()
    dx_vals = domain_length / K_vals
    l2_norms = n_data['l2_norm'].to_numpy()
    
    # Filter out numerical zeros to keep log-log scaling stable
    valid_idx = l2_norms > 1e-16
    dx_valid = dx_vals[valid_idx]
    l2_valid = l2_norms[valid_idx]
    
    if len(dx_valid) == 0:
        continue
        
    color = colors[idx]
    
    # Plot Empirical Data (Solid line with markers across ALL valid points)
    plt.loglog(dx_valid, l2_valid, marker='o', linestyle='-', 
               color=color, linewidth=2.0, markersize=6, label=f'N = {N}')
    
    # --- ASYMPTOTIC REGION MANIPULATION ---
    # Isolate exactly the last 3 valid data points (smallest delta_x values)
    num_asymptotic_pts = min(3, len(dx_valid))
    dx_asymptotic = dx_valid[-num_asymptotic_pts:]
    l2_asymptotic = l2_valid[-num_asymptotic_pts:]
    
    # Generate the Theoretical Rate: O(delta_x^(N+1))
    # Anchor the rate line at the first point of this 3-point asymptotic window
    expected_slope = N + 1
    anchor_dx = dx_asymptotic[0]
    anchor_l2 = l2_asymptotic[0]
    
    theoretical_baseline = anchor_l2 * (dx_asymptotic / anchor_dx)**expected_slope
    
    # Plot dashed theoretical rate ONLY over these last 3 points
    plt.loglog(dx_asymptotic, theoretical_baseline, linestyle='--', 
               color=color, linewidth=2.5, alpha=0.6)

# 4. Axes Formatting & Grid Setup
plt.xscale('log', base=2)
plt.yscale('log', base=10)

plt.xlabel(r'Element Size ($\Delta x = 2\pi / K$)', fontsize=11)
plt.ylabel(r'$L_2$ Norm Error', fontsize=11)
plt.title(r'$L_2$ Norm Convergence vs. Element Size ($\Delta x$)', fontsize=12, pad=15)

# Enforce clean grid lines for major/minor levels
plt.grid(True, which="both", linestyle=":", alpha=0.5, color='gray')

# Construct custom clean labels for the dual-state legend entries
handles, labels = plt.gca().get_legend_handles_labels()
data_handle = Line2D([], [], color='black', linestyle='-', marker='o', label='Empirical Data')
rate_handle = Line2D([], [], color='black', linestyle='--', linewidth=2.5, alpha=0.6, label='Asymptotic Rate $O(\Delta x^{N+1})$')

# Final plot compilation
plt.legend(handles=handles + [data_handle, rate_handle], loc='lower right', frameon=True, fontsize=9)
plt.tight_layout()

# Save options (useful for headless servers/SSH connections)
plt.savefig('l2_norm_vs_deltax.svg', format='svg')
plt.show()

