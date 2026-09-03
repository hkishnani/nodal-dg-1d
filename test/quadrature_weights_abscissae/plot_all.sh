#!/usr/bin/env bash
set -e

# ==============================================================================
# USER INPUT CHECK
# ==============================================================================
if [ -z "$1" ]; then
    echo "Usage: $0 <target_directory> [backend: matplotlib|plotly] [extension: png|svg|html]"
    echo "Example: $0 ../quadrature_weights_abscissae/ matplotlib svg GL"
    exit 1
fi

TARGET_DIR="${1:-./}"              # Default to current Directory
BACKEND="${2:-matplotlib}"         # Default to 'matplotlib'
EXT="${3:-svg}"                    # Default to 'png'

if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' does not exist."
    exit 1
fi

# ==============================================================================
# CONFIGURATION
# ==============================================================================
SCRIPT="plot_csv.py"               # Path to Python script
KIND="both"                        # 'line', 'scatter', or 'both'

# Math / Presentation Labels
X_LABEL='$\zeta$'
Y_LABEL='$w(\zeta)$'

# ==============================================================================
# EXECUTION
# ==============================================================================
echo "Processing CSV files in: $TARGET_DIR"
echo "Backend: $BACKEND | Format: $EXT"
echo "--------------------------------------------------------"

for csv_file in "$TARGET_DIR"/*.csv; do
    [ -e "$csv_file" ] || { echo "No CSV files found in $TARGET_DIR"; exit 1; }

    filename=$(basename -- "$csv_file")
    basename_no_ext="${filename%.*}"
    
    # Dynamically select the Y column based on file prefix
    if [[ "$QUAD_TYPE" == GL ]]; then
        X_COL="zeta_GL"
        Y_COL="w_GL"
        TITLE_PREFIX="Gauss-Legendre quadrature weights vs abscissae"
    elif [[ "$QUAD_TYPE" == GLL ]]; then
        X_COL="zeta_GLL"
        Y_COL="w_GLL"
        TITLE_PREFIX="Gauss-Lobatto-Legendre weights and abscissae"
    else
        # Fallback to column index 1 (2nd column) if standard names don't match
        X_COL="0"
        Y_COL="1"
        TITLE_PREFIX="Lagrange Basis Polynomials"
    fi

    out_file="${TARGET_DIR}/${basename_no_ext}.${EXT}"
    title="${TITLE_PREFIX} (${basename_no_ext})"

    echo "Plotting: $filename ($Y_COL) -> ${basename_no_ext}.${EXT}"

    ARGS=(
        --file "$csv_file"
        --x "$X_COL"
        --y "$Y_COL"
        --kind "$KIND"
        --backend "$BACKEND"
        --output "$out_file"
        --xlabel "$X_LABEL"
        --ylabel "$Y_LABEL"
        --title "$title"
    )

    python3 "$SCRIPT" "${ARGS[@]}"
done

echo "--------------------------------------------------------"
echo "Done! All plots generated inside '$TARGET_DIR/'."
