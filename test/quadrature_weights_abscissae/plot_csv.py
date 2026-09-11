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

def format_decimal(val) -> str:
    """Format numeric values as normal decimal notation with up to 5 decimal places."""
    try:
        return f"{float(val):.5f}"
    except (ValueError, TypeError):
        return str(val)

def plot_matplotlib(df, x_col, y_col, group_col, out_path, title, xlabel, ylabel):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Presentation/PPT Styling Defaults
    plt.rcParams.update({
        "mathtext.fontset": "cm",
        "font.size": 12,
        "axes.labelsize": 14,
        "axes.titlesize": 16,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "legend.fontsize": 11,
        "figure.titlesize": 18
    })

    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    def annotate_points(sub_df):
        """Helper to annotate decimal coordinates near each bar top."""
        for x, y in zip(sub_df[x_col], sub_df[y_col]):
            if pd.notnull(x) and pd.notnull(y):
                x_str = format_decimal(x)
                y_str = format_decimal(y)
                
                ax.annotate(
                    f"({x_str}, {y_str})",
                    (x, y),
                    xytext=(5, 5),
                    textcoords="offset points",
                    ha="left",
                    va="bottom",
                    fontsize=8,
                    alpha=0.8,
                    zorder=3
                )

    if group_col:
        groups = df.groupby(group_col, sort=False)
        for val, group in groups:
            group_sorted = group.sort_values(by=x_col)
            label = f"{group_col} = {val}"
            # Render as bar graph
            ax.bar(
                group_sorted[x_col],
                group_sorted[y_col],
                width=0.015,
                label=label,
                alpha=0.7,
                edgecolor="black",
                zorder=2
            )
            annotate_points(group_sorted)
    else:
        df_sorted = df.sort_values(by=x_col)
        # Render as bar graph
        ax.bar(
            df_sorted[x_col],
            df_sorted[y_col],
            width=0.015,
            label=ylabel or y_col,
            alpha=0.8,
            color="royalblue",
            edgecolor="black",
            zorder=2
        )
        annotate_points(df_sorted)

    clean_x = xlabel if xlabel else x_col
    clean_y = ylabel if ylabel else y_col
    default_title = f"{clean_y} vs {clean_x}" + (f" (grouped by {group_col})" if group_col else "")
    clean_title = title if title else default_title

    ax.set_xlabel(clean_x, fontweight="bold")
    ax.set_ylabel(clean_y, fontweight="bold")
    ax.set_title(clean_title, pad=12, fontweight="bold")
    
    # Fix X-axis limits to [-1.0, 1.0] (padded slightly so bars don't clip at the edges)
    ax.set_xlim(-1.05, 1.05)

    ax.grid(True, linestyle="--", alpha=0.6, linewidth=0.8, zorder=0)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", frameon=True, framealpha=0.95)

    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

def plot_plotly(df, x_col, y_col, group_col, out_path, title, xlabel, ylabel):
    import plotly.express as px

    if group_col:
        df_sorted = df.sort_values(by=[group_col, x_col])
    else:
        df_sorted = df.sort_values(by=x_col)

    clean_x = xlabel if xlabel else x_col
    clean_y = ylabel if ylabel else y_col
    final_title = title if title else (f"{clean_y} vs {clean_x}" + (f" (grouped by {group_col})" if group_col else ""))

    # Render as bar graph using px.bar
    if group_col:
        fig = px.bar(df_sorted, x=x_col, y=y_col, color=group_col, barmode="group", title=final_title)
    else:
        fig = px.bar(df_sorted, x=x_col, y=y_col, title=final_title)

    # Add annotations for coordinates
    for _, row in df_sorted.iterrows():
        x, y = row[x_col], row[y_col]
        if pd.notnull(x) and pd.notnull(y):
            x_str = format_decimal(x)
            y_str = format_decimal(y)
            
            fig.add_annotation(
                x=x,
                y=y,
                text=f"({x_str}, {y_str})",
                showarrow=False,
                xanchor="left",
                yanchor="bottom",
                xshift=5,
                yshift=5,
                font=dict(size=10, color="gray")
            )

    fig.update_layout(
        template="plotly_white",
        font=dict(size=14, family="Arial, sans-serif"),
        title=dict(font=dict(size=20, color="black")),
        xaxis=dict(
            title=dict(text=clean_x, font=dict(size=16, color="black")),
            showgrid=True, gridwidth=1, gridcolor="LightGray",
            type="linear", # Enforces proper spacing for precise abscissae coordinates
            range=[-1.05, 1.05] # Lock X-axis bounds
        ),
        yaxis=dict(
            title=dict(text=clean_y, font=dict(size=16, color="black")),
            showgrid=True, gridwidth=1, gridcolor="LightGray"
        ),
        bargap=0.85 # Makes bars thin to reflect specific coordinate locations clearly
    )

    html_out = out_path.with_suffix(".html")
    fig.write_html(str(html_out), include_mathjax="cdn")
    return html_out

def main():
    parser = argparse.ArgumentParser(
        description="Plot X vs Y bar graphs from a CSV file with presentation styling."
    )
    parser.add_argument("--file", "-f", type=str, required=True, help="Path to input CSV file.")
    parser.add_argument("--x", type=str, required=True, help="Column name or index for X-axis (Abscissae).")
    parser.add_argument("--y", type=str, required=True, help="Column name or index for Y-axis (Weights).")
    parser.add_argument("--group-by", "-g", type=str, default=None, help="Column name or index for group legends.")
    parser.add_argument("--title", type=str, default=None, help="Custom figure title.")
    parser.add_argument("--xlabel", type=str, default=None, help="Custom X-axis label.")
    parser.add_argument("--ylabel", type=str, default=None, help="Custom Y-axis label.")
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
        plot_matplotlib(df, x_col, y_col, group_col, out_path, args.title, args.xlabel, args.ylabel)
        print(f"Successfully saved Matplotlib presentation plot -> {out_path}", flush=True)
    else:
        final_out = plot_plotly(df, x_col, y_col, group_col, out_path, args.title, args.xlabel, args.ylabel)
        print(f"Successfully saved Plotly interactive HTML -> {final_out}", flush=True)

if __name__ == "__main__":
    main()