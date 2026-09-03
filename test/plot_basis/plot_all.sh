#!/usr/bin/env bash
set -e

# ==============================================================================
# USER INPUT CHECK
# ==============================================================================
if [ -z "$1" ]; then
    echo "Usage: $0 <target_directory> [backend: matplotlib|plotly] [extension: png|svg|html]"
    echo "Example: $0 plot_basis matplotlib png"
    exit 1
fi

TARGET_DIR="$1"
BACKEND="${2:-matplotlib}"         # Default to 'matplotlib'
EXT="${3:-png}"                    # Default to 'png'

if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' does not exist."
    exit 1
fi

# ==============================================================================
# CONFIGURATION
# ==============================================================================
SCRIPT="plot_csv.py"               # Path to Python script
X_COL="eta"                        # X-axis column name
GROUP_BY="j"                       # Grouping column
KIND="both"                        # 'line', 'scatter', or 'both'

# Math / Presentation Labels
X_LABEL='$\eta$'
Y_LABEL='$l_j(\eta)$'

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
    if [[ "$filename" == GL_basis_* ]]; then
        Y_COL="l_j_GL"
        TITLE_PREFIX="Gauss-Legendre Basis Polynomials"
    elif [[ "$filename" == GLL_basis_* ]]; then
        Y_COL="l_j_GLL"
        TITLE_PREFIX="Gauss-Lobatto-Legendre Basis Polynomials"
    else
        # Fallback to column index 3 (4th column) if standard names don't match
        Y_COL="3"
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

    if [ -n "$GROUP_BY" ]; then
        ARGS+=(--group-by "$GROUP_BY")
    fi

    python3 "$SCRIPT" "${ARGS[@]}"
done

echo "--------------------------------------------------------"
echo "Done! All plots generated inside '$TARGET_DIR/'."