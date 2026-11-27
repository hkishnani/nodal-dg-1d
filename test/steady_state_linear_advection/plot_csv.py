#!/usr/bin/env python3

import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt

from concurrent.futures import ProcessPoolExecutor


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

DATA_DIR = "."

OUTPUT_DIR = "svg"

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ------------------------------------------------------------
# Extract K and N from filename
# ------------------------------------------------------------

def get_KN(filename):

    basename = os.path.basename(filename)

    match = re.search(r"K(\d+)_N(\d+)", basename)

    if match is None:
        raise ValueError(
            f"Could not determine K and N from {filename}"
        )

    K = int(match.group(1))
    N = int(match.group(2))

    return K, N


# ------------------------------------------------------------
# Find all CSV files
# ------------------------------------------------------------

CSV_FILES = glob.glob(
    os.path.join(DATA_DIR, "K*_N*.csv")
)

CSV_FILES.sort(
    key=lambda f: get_KN(f)
)


# ------------------------------------------------------------
# Find highest-fidelity reference solution
# ------------------------------------------------------------

REFERENCE_FILE = max(
    CSV_FILES,
    key=lambda f: get_KN(f)
)

K_REF, N_REF = get_KN(REFERENCE_FILE)


# ------------------------------------------------------------
# Read reference solution
#
# This is done once in the main process and passed to workers.
# ------------------------------------------------------------

REFERENCE_DF = pd.read_csv(REFERENCE_FILE)

REFERENCE_DF.columns = REFERENCE_DF.columns.str.strip()

X_REF = REFERENCE_DF["x"].to_numpy()
U0_REF = REFERENCE_DF["u0"].to_numpy()


# ------------------------------------------------------------
# Plot one CSV file
# ------------------------------------------------------------

def plot_csv(filename):

    df = pd.read_csv(filename)

    df.columns = df.columns.str.strip()

    K, N = get_KN(filename)

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    fig, ax = plt.subplots(figsize=(10, 6))

    # --------------------------------------------------------
    # Common high-fidelity reference solution
    # --------------------------------------------------------

    ax.plot(
        X_REF,
        U0_REF,
        linewidth=2.0,
        label=(
            rf"Reference $u_0$ "
            rf"($K={K_REF},\,N={N_REF}$)"
        )
    )

    # --------------------------------------------------------
    # DG solution
    # --------------------------------------------------------

    first_element = True

    for i, element in df.groupby("i"):

        ax.plot(
            element["x"],
            element["uh"],
            "o-",
            markersize=4,
            linewidth=1.2,
            label=r"DG $u_h$" if first_element else None
        )

        first_element = False

    # --------------------------------------------------------
    # Formatting
    # --------------------------------------------------------

    ax.set_xlabel("$x$")
    ax.set_ylabel("$u$")

    ax.set_title(
        rf"$L^2$ Projection: $K={K}$, $N={N}$"
    )

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    # --------------------------------------------------------
    # Save SVG
    # --------------------------------------------------------

    output_file = os.path.join(
        OUTPUT_DIR,
        f"K{K}_N{N}.svg"
    )

    fig.savefig(
        output_file,
        format="svg",
        bbox_inches="tight"
    )

    plt.close(fig)

    return f"K{K}_N{N}.svg"


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

if __name__ == "__main__":

    print()
    print(f"Found {len(CSV_FILES)} CSV files.")

    print(
        f"Reference: K{K_REF}_N{N_REF}.csv"
    )

    # --------------------------------------------------------
    # Number of CPU cores
    # --------------------------------------------------------

    num_workers = os.cpu_count()

    print(
        f"Using {num_workers} parallel workers."
    )

    print()

    # --------------------------------------------------------
    # Parallel plotting
    # --------------------------------------------------------

    with ProcessPoolExecutor(
        max_workers=num_workers
    ) as executor:

        results = executor.map(
            plot_csv,
            CSV_FILES
        )

        for result in results:
            print(f"Saved: {result}")

    print()
    print(
        f"All plots saved in: {OUTPUT_DIR}/"
    )