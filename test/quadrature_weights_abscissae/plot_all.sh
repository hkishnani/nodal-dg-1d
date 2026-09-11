#!/usr/bin/env bash
set -e

MAX_JOBS=16

# ==============================================================================
# USER INPUT CHECK
# ==============================================================================
if [ -z "$1" ]; then
    echo "Usage: $0 <target_directory> [backend: matplotlib|plotly] [extension: png|svg|html]"
    echo "Example: $0 ../quadrature_weights_abscissae/ matplotlib svg"
    exit 1
fi

TARGET_DIR="${1:-./}"              # Default to current Directory
BACKEND="${2:-matplotlib}"         # Default to 'matplotlib'
EXT="${3:-svg}"                    # Default to 'svg'

if [ ! -d "$TARGET_DIR" ]; then
    echo "Error: Directory '$TARGET_DIR' does not exist."
    exit 1
fi

# ==============================================================================
# CONFIGURATION
# ==============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/plot_csv.py" # Path to Python script

# Math / Presentation Labels
X_LABEL='$\zeta$'
Y_LABEL='$w(\zeta)$'

# ==============================================================================
# EXECUTION
# ==============================================================================
echo "Processing CSV files in: $TARGET_DIR"
echo "Backend: $BACKEND | Format: $EXT"
echo "--------------------------------------------------------"

# Count CSV files
csv_files=("$TARGET_DIR"/*.csv)
if [ ${#csv_files[@]} -eq 0 ] || [ ! -e "${csv_files[0]}" ]; then
    echo "No CSV files found in $TARGET_DIR"
    exit 1
fi

for csv_file in "${csv_files[@]}"; do
    [ -e "$csv_file" ] || continue

    filename=$(basename -- "$csv_file")
    basename_no_ext="${filename%.*}"
    
    # Dynamically select the Y column based on file prefix
    if [[ "$filename" == *"GL"* && "$filename" != *"GLL"* ]]; then
        X_COL="zeta_GL"
        Y_COL="w_GL"
        TITLE_PREFIX="Gauss-Legendre Quadrature Weights vs Abscissae"
    elif [[ "$filename" == *"GLL"* ]]; then
        X_COL="zeta_GLL"
        Y_COL="w_GLL"
        TITLE_PREFIX="Gauss-Lobatto-Legendre Weights vs Abscissae"
    else
        # Fallback to column index 0 and 1
        X_COL="0"
        Y_COL="1"
        TITLE_PREFIX="Weights vs Abscissae"
    fi

    out_file="${TARGET_DIR}/${basename_no_ext}.${EXT}"
    title="${TITLE_PREFIX} (${basename_no_ext})"

    ARGS=(
        --file "$csv_file"
        --x "$X_COL"
        --y "$Y_COL"
        --backend "$BACKEND"
        --output "$out_file"
        --xlabel "$X_LABEL"
        --ylabel "$Y_LABEL"
        --title "$title"
    )

    # Execute Python script in the background for parallel processing
    python3 "$SCRIPT" "${ARGS[@]}" &

    # Throttle job concurrency
    while [ $(jobs -p | wc -l) -ge $MAX_JOBS ]; do
        sleep 0.1
    done
done

# Wait for remaining background processes
wait

echo "--------------------------------------------------------"
echo "Done! All plots generated inside '$TARGET_DIR/'."