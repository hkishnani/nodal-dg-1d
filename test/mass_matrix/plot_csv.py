#!/usr/bin/env python3

import sys
import argparse
from pathlib import Path
import pandas as pd

def resolve_column(df: pd.DataFrame, col_arg: str) -> str:
    """Resolve a column argument (either name or 0-based index) to a valid DataFrame column name."""
    if col_arg in df.columns:
        return col_arg
    try:
        idx = int(col_arg)
        if 0 <= idx < len(df.columns):
            return df.columns[idx]
    except ValueError:
        pass
    raise ValueError(f"Column '{col_arg}' not found in CSV. Available columns: {list(df.columns)}")

def plot_matplotlib(df, x_col, y_col, group_col, kind, out_path, title, xlabel, ylabel):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Presentation/PPT Styling Defaults
    plt.rcParams.update({
        "mathtext.fontset": "cm",  # Computer Modern font for LaTeX math
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "figure.titlesize": 18
    })

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    show_lines = kind in ["line", "both"]
    marker = "o" if kind in ["scatter", "both"] else None
    linestyle = "-" if show_lines else "None"

    if group_col:
        groups = df.groupby(group_col, sort=False)
        for val, group in groups:
            group_sorted = group.sort_values(by=x_col)
            label = f"{group_col} = {val}"
            ax.plot(
                group_sorted[x_col],
                group_sorted[y_col],
                label=label,
                linestyle=linestyle,
                marker=marker,
                linewidth=2.5 if show_lines else 0,
                markersize=6.0,
                zorder=2
            )
    else:
        df_sorted = df.sort_values(by=x_col)
        ax.plot(
            df_sorted[x_col],
            df_sorted[y_col],
            label=ylabel or y_col,
            linestyle=linestyle,
            marker=marker,
            linewidth=2.5 if show_lines else 0,
            markersize=6.0,
            zorder=2
        )

    clean_x = xlabel if xlabel else x_col
    clean_y = ylabel if ylabel else y_col
    default_title = f"{clean_y} vs {clean_x}" + (f" (grouped by {group_col})" if group_col else "")
    clean_title = title if title else default_title

    ax.set_xlabel(clean_x, fontweight="bold")
    ax.set_ylabel(clean_y, fontweight="bold")
    ax.set_title(clean_title, pad=12, fontweight="bold")

    ax.grid(True, linestyle="--", alpha=0.6, linewidth=0.8)
    # Fixed facealpha -> framealpha
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_plotly(df, x_col, y_col, group_col, kind, out_path, title, xlabel, ylabel):
    import plotly.express as px

    mode_map = {
        "line": "lines",
        "scatter": "markers",
        "both": "lines+markers"
    }

    if group_col:
        df_sorted = df.sort_values(by=[group_col, x_col])
    else:
        df_sorted = df.sort_values(by=x_col)

    clean_x = xlabel if xlabel else x_col
    clean_y = ylabel if ylabel else y_col
    final_title = title if title else (f"{clean_y} vs {clean_x}" + (f" (grouped by {group_col})" if group_col else ""))

    if kind == "scatter" and not group_col:
        fig = px.scatter(df_sorted, x=x_col, y=y_col, title=final_title)
    elif kind == "scatter" and group_col:
        fig = px.scatter(df_sorted, x=x_col, y=y_col, color=group_col, title=final_title)
    elif group_col:
        fig = px.line(df_sorted, x=x_col, y=y_col, color=group_col, title=final_title)
        fig.update_traces(mode=mode_map[kind])
    else:
        fig = px.line(df_sorted, x=x_col, y=y_col, title=final_title)
        fig.update_traces(mode=mode_map[kind])

    fig.update_layout(
        template="plotly_white",
        font=dict(size=14, family="Arial, sans-serif"),
        title=dict(font=dict(size=20, color="black")),
        xaxis=dict(
            title=dict(text=clean_x, font=dict(size=16, color="black")),
            showgrid=True, gridwidth=1, gridcolor="LightGray"
        ),
        yaxis=dict(
            title=dict(text=clean_y, font=dict(size=16, color="black")),
            showgrid=True, gridwidth=1, gridcolor="LightGray"
        )
    )
    
    fig.update_traces(line=dict(width=3), marker=dict(size=8))

    html_out = out_path.with_suffix(".html")
    fig.write_html(str(html_out), include_mathjax="cdn")
    return html_out

def main():
    parser = argparse.ArgumentParser(
        description="Plot X vs Y from any CSV file with presentation styling and TeX math rendering."
    )
    parser.add_argument("--file", "-f", type=str, required=True, help="Path to input CSV file.")
    parser.add_argument("--x", type=str, required=True, help="Column name or index for X-axis.")
    parser.add_argument("--y", type=str, required=True, help="Column name or index for Y-axis.")
    parser.add_argument("--group-by", "-g", type=str, default=None, help="Column name or index for group legends.")
    parser.add_argument("--title", type=str, default=None, help="Custom figure title.")
    parser.add_argument("--xlabel", type=str, default=None, help="Custom X-axis label.")
    parser.add_argument("--ylabel", type=str, default=None, help="Custom Y-axis label.")
    parser.add_argument(
        "--kind", "-k",
        choices=["line", "scatter", "both"],
        default="line",
        help="Plot style: 'line', 'scatter', or 'both' (default: line)."
    )
    parser.add_argument(
        "--backend", "-b",
        choices=["matplotlib", "plotly"],
        default="matplotlib",
        help="Backend: 'matplotlib' (PNG/SVG) or 'plotly' (HTML) (default: matplotlib)."
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default=None,
        help="Output filename."
    )
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    
    csv_path = Path(args.file).resolve()
    if not csv_path.exists():
        raise FileNotFoundError(f"File not found: {csv_path}")

    default_filename = "plot.html" if args.backend == "plotly" else "plot.png"
    out_arg = Path(args.output) if args.output else Path(default_filename)
    
    if out_arg.is_absolute():
        out_path = out_arg
    else:
        out_path = csv_path.parent / out_arg.name

    df = pd.read_csv(csv_path)
    
    x_col = resolve_column(df, args.x)
    y_col = resolve_column(df, args.y)
    group_col = resolve_column(df, args.group_by) if args.group_by is not None else None

    if args.backend == "matplotlib":
        plot_matplotlib(df, x_col, y_col, group_col, args.kind, out_path, args.title, args.xlabel, args.ylabel)
        print(f"Successfully saved Matplotlib presentation plot -> {out_path}")
    else:
        final_out = plot_plotly(df, x_col, y_col, group_col, args.kind, out_path, args.title, args.xlabel, args.ylabel)
        print(f"Successfully saved Plotly interactive HTML -> {final_out}")

if __name__ == "__main__":
    main()