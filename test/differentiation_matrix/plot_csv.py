from concurrent.futures import ProcessPoolExecutor
import os
import re
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

def process_single_csv(args):
    filepath, basis_dir, quadrature_dir, out_path, cmap, save_format = args
    try:
        matrix = np.loadtxt(filepath, delimiter=",", ndmin=2)
        rows, cols = matrix.shape

        match = re.search(
            r"Dqj_matrix_N_(\d+)_Nq_(\d+)_basis_(GLL|GL)_quadrature_(GLL|GL)",
            filepath.name,
        )
        if match is None:
            raise ValueError(
                "Expected filename pattern Dqj_matrix_N_<N>_Nq_<Nq>_basis_<GL|GLL>_quadrature_<GL|GLL>.csv"
            )

        polynomial_order, quadrature_order, basis_family, quadrature_family = match.groups()
        basis_order = int(polynomial_order) + 1
        basis_path = basis_dir / f"{basis_family}_basis_Nq_{basis_order}.csv"
        if not basis_path.exists():
            raise FileNotFoundError(f"Matching basis CSV not found: {basis_path}")

        quadrature_path = (
            quadrature_dir
            / f"{quadrature_family}_abscissae_weights{quadrature_order}.csv"
        )
        if not quadrature_path.exists():
            raise FileNotFoundError(
                f"Matching quadrature CSV not found: {quadrature_path}"
            )

        basis_df = pd.read_csv(basis_path)
        basis_value_column = f"l_j_{basis_family}"
        required_columns = {"eta", "j", basis_value_column}
        missing_columns = required_columns.difference(basis_df.columns)
        if missing_columns:
            raise ValueError(
                f"Basis CSV {basis_path.name} is missing columns: {sorted(missing_columns)}"
            )

        quadrature_df = pd.read_csv(quadrature_path)
        abscissae_column = f"zeta_{quadrature_family}"
        weights_column = f"w_{quadrature_family}"
        required_columns = {abscissae_column, weights_column}
        missing_columns = required_columns.difference(quadrature_df.columns)
        if missing_columns:
            raise ValueError(
                f"Quadrature CSV {quadrature_path.name} is missing columns: "
                f"{sorted(missing_columns)}"
            )

        fig_width = max(12, cols * 0.7 + 5)
        fig_height = max(7, rows * 0.7)
        fig, (matrix_ax, basis_ax) = plt.subplots(
            1,
            2,
            figsize=(fig_width, fig_height),
            gridspec_kw={"width_ratios": [1.2, 1.35]},
        )

        # 1. RENDER COLORFUL BACKGROUND (USING MATRIX VALUES)
        norm = Normalize(vmin=matrix.min(), vmax=matrix.max())
        matrix_ax.imshow(matrix, cmap=cmap, norm=norm, aspect="equal")

        # 2. ADD TEXT ANNOTATIONS (With dynamic text color for contrast)
        def get_text_color(val, cmap_instance, norm_instance):
            rgb = cmap_instance(norm_instance(val))[:3]
            luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
            return "white" if luminance < 0.5 else "black"

        font_size = max(5, min(10, int(80 / max(rows, cols))))

        for i in range(rows):
            for j in range(cols):
                val = matrix[i, j]
                
                if abs(val) > 0 and abs(val) < 1e-3:
                    text_str = f"{val:.2e}"
                else:
                    text_str = f"{val:.4f}".rstrip("0").rstrip(".")

                txt_color = get_text_color(val, plt.get_cmap(cmap), norm)
                
                matrix_ax.text(
                    j,
                    i,
                    text_str,
                    ha="center",
                    va="center",
                    fontsize=font_size,
                    color=txt_color,
                )

        # 3. Plot Formatting (grid, titles, ticks)
        matrix_ax.set_xticks(np.arange(cols) - 0.5, minor=True)
        matrix_ax.set_yticks(np.arange(rows) - 0.5, minor=True)
        matrix_ax.grid(which="minor", color="gray", linestyle="-", linewidth=0.7)
        matrix_ax.tick_params(which="minor", size=0)
        matrix_ax.set_xticks(np.arange(cols))
        matrix_ax.set_yticks(np.arange(rows))
        
        title = filepath.stem.replace("Dqj_matrix_", "").replace("_", " ")
        matrix_ax.set_title(f"Dqj Matrix: {title}", fontsize=10, pad=10)
        matrix_ax.set_xlabel("Basis Index (j)")
        matrix_ax.set_ylabel("Basis Index (i)")

        # Plot all Lagrange basis polynomials associated with this Dqj matrix.
        for basis_index, group in basis_df.groupby("j", sort=True):
            group = group.sort_values("eta")
            basis_ax.plot(
                group["eta"],
                group[basis_value_column],
                linewidth=1.8,
                marker="o",
                markersize=2.5,
                label=fr"$l_{{{int(basis_index)}}}(\eta)$",
            )

        basis_y_max = basis_ax.get_ylim()[1]
        max_weight = quadrature_df[weights_column].max()
        weight_scale = 0.85 * basis_y_max / max_weight
        bar_width = 0.08 if len(quadrature_df) > 1 else 0.12
        quadrature_bars = basis_ax.bar(
            quadrature_df[abscissae_column],
            quadrature_df[weights_column] * weight_scale,
            width=bar_width,
            alpha=0.28,
            color="tab:orange",
            edgecolor="tab:orange",
            linewidth=0.8,
            label="_nolegend_",
            zorder=0,
        )

        basis_ax.set_title(
            f"{basis_family} Basis + {quadrature_family} Weights "
            f"(N = {polynomial_order}, Nq = {quadrature_order})",
            fontsize=10,
            pad=10,
        )
        basis_ax.set_xlabel(r"$\eta$")
        basis_ax.set_ylabel(r"$l_j(\eta)$")
        quadrature_ax = basis_ax.secondary_yaxis(
            "right",
            functions=(
                lambda basis_value: basis_value / weight_scale,
                lambda weight: weight * weight_scale,
            ),
        )
        quadrature_ax.set_ylabel(r"Quadrature weight $w$")
        unique_weights = np.sort(quadrature_df[weights_column].unique())
        quadrature_ax.set_yticks(unique_weights)
        quadrature_xaxis = basis_ax.secondary_xaxis(
            "top",
            functions=(lambda eta: eta, lambda zeta: zeta),
        )
        quadrature_xaxis.set_xticks(quadrature_df[abscissae_column])
        quadrature_xaxis.set_xlabel(r"Quadrature abscissae $\zeta$")
        quadrature_xaxis.tick_params(axis="x", labelsize=7, pad=1)
        for tick_label in quadrature_xaxis.get_xticklabels():
            tick_label.set_rotation(45)
            tick_label.set_horizontalalignment("left")
        basis_ax.grid(True, linestyle="--", alpha=0.6, linewidth=0.8)
        basis_ax.axhline(0, color="black", linewidth=0.7)
        basis_handles, basis_labels = basis_ax.get_legend_handles_labels()
        quadrature_handles, quadrature_labels = [quadrature_bars], [
            f"{quadrature_family} quadrature weights ($Nq={quadrature_order}$)"
        ]
        fig.legend(
            basis_handles + quadrature_handles,
            basis_labels + quadrature_labels,
            loc="upper center",
            bbox_to_anchor=(0.28, 0.86),
            fontsize=7,
            ncol=4,
            framealpha=0.8,
        )

        plot_path = out_path / f"{filepath.stem}.{save_format}"
        fig.tight_layout(rect=[0, 0, 1, 0.81])
        fig.savefig(plot_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        print(f"Error processing {filepath.name}: {e}")

def plot_all_Dqj_matrices(
    source_dir=".",
    basis_dir="../plot_basis",
    quadrature_dir="../quadrature_weights_abscissae",
    output_dir="plots",
    cmap="magma",
    save_format="svg",
    max_workers=32,
):
    src_path = Path(source_dir)
    basis_path = Path(basis_dir)
    quadrature_path = Path(quadrature_dir)
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    csv_files = list(src_path.glob("Dqj_matrix_*.csv"))

    if not csv_files:
        print(f"No Dqj matrix CSV files found in '{src_path.resolve()}'.")
        return

    available_cores = os.cpu_count() or max_workers
    workers = min(max_workers, available_cores)

    print(
        f"Processing {len(csv_files)} files in parallel using {workers} CPU workers..."
    )

    tasks = [
        (f, basis_path, quadrature_path, out_path, cmap, save_format)
        for f in csv_files
    ]

    with ProcessPoolExecutor(max_workers=workers) as executor:
        list(executor.map(process_single_csv, tasks))

    print(f"Successfully generated all plots in '{out_path.name}/'.")

if __name__ == "__main__":
    plot_all_Dqj_matrices(max_workers=32, cmap="magma")
