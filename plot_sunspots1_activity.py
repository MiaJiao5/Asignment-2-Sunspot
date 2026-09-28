# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib"]
# ///

r"""
Summarize daily SILSO sunspot counts by year and save a relationship matrix.

    uv run .\plot_sunspots.py
"""

import csv
import math
import random
import statistics
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

FILE = "SN_d_tot_V2.0.csv"
PICTURE = "sunspot_relationship_activity.png"
MIN_VALID_DAYS = 300

HERE = Path(__file__).parent
DATA = HERE / "data" / FILE
OUT = HERE / "out"

FEATURES = [
    ("Annual mean", "Mean daily sunspot number"),
    ("Annual median", "Median daily sunspot number"),
    ("Annual max", "Maximum daily sunspot number"),
    ("Daily SD", "Standard deviation of daily counts"),
    ("Daily IQR", "Interquartile range of daily counts"),
    ("90th pct.", "90th percentile of daily counts"),
    ("Spotless %", "Percentage of valid days with zero sunspots"),
    ("High activity %", "Percentage of valid days with count >= 100"),
]

PURPLE = "#8f8298"
LIGHT_PURPLE = "#c9bdcf"
GREEN = "#78998b"
LIGHT_GREEN = "#b9cabe"
NEUTRAL = "#f3f0eb"
INK = "#37343b"
MUTED = "#77727b"
CORRELATION_CMAP = LinearSegmentedColormap.from_list(
    "sage_plum", [PURPLE, NEUTRAL, GREEN]
)


def annual_metrics(path):
    daily_counts = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle, delimiter=";"):
            if len(row) < 5:
                continue
            try:
                year = int(row[0].strip())
                count = float(row[4].strip())
            except ValueError:
                continue
            if count >= 0:
                daily_counts.setdefault(year, []).append(count)

    years = []
    values = {label: [] for label, _ in FEATURES}
    for year, counts in sorted(daily_counts.items()):
        if len(counts) < MIN_VALID_DAYS:
            continue
        quartiles = statistics.quantiles(counts, n=4, method="inclusive")
        deciles = statistics.quantiles(counts, n=10, method="inclusive")
        years.append(year)
        values["Annual mean"].append(statistics.fmean(counts))
        values["Annual median"].append(statistics.median(counts))
        values["Annual max"].append(max(counts))
        values["Daily SD"].append(statistics.pstdev(counts))
        values["Daily IQR"].append(quartiles[2] - quartiles[0])
        values["90th pct."].append(deciles[8])
        values["Spotless %"].append(100 * sum(count == 0 for count in counts) / len(counts))
        values["High activity %"].append(100 * sum(count >= 100 for count in counts) / len(counts))
    return years, values, len(daily_counts)


def average_ranks(values):
    ordered = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(ordered):
        end = start + 1
        while end < len(ordered) and ordered[end][1] == ordered[start][1]:
            end += 1
        rank = ((start + 1) + end) / 2
        for position in range(start, end):
            ranks[ordered[position][0]] = rank
        start = end
    return ranks


def spearman_correlation(left, right):
    left_ranks = average_ranks(left)
    right_ranks = average_ranks(right)
    left_mean = statistics.fmean(left_ranks)
    right_mean = statistics.fmean(right_ranks)
    covariance = sum(
        (left_rank - left_mean) * (right_rank - right_mean)
        for left_rank, right_rank in zip(left_ranks, right_ranks)
    )
    left_sum = sum((rank - left_mean) ** 2 for rank in left_ranks)
    right_sum = sum((rank - right_mean) ** 2 for rank in right_ranks)
    denominator = math.sqrt(left_sum * right_sum)
    return covariance / denominator if denominator else 0.0


def kernel_density(values, points=100):
    spread = statistics.pstdev(values)
    if spread == 0:
        return [], []
    bandwidth = 1.06 * spread * len(values) ** (-0.2)
    low = min(values) - 3 * bandwidth
    high = max(values) + 3 * bandwidth
    x_values = [low + (high - low) * index / (points - 1) for index in range(points)]
    scale = len(values) * bandwidth * math.sqrt(2 * math.pi)
    y_values = [
        sum(math.exp(-0.5 * ((x_value - value) / bandwidth) ** 2) for value in values) / scale
        for x_value in x_values
    ]
    return x_values, y_values


def padded_limits(values):
    low = min(values)
    high = max(values)
    padding = (high - low) * 0.06 or 1
    return low - padding, high + padding


def main():
    years, metrics, total_years = annual_metrics(DATA)
    if not years:
        raise ValueError(f"No years have at least {MIN_VALID_DAYS} valid daily observations")

    labels = [label for label, _ in FEATURES]
    descriptions = [description for _, description in FEATURES]
    print(f"{DATA.name}: {len(years)} usable years ({years[0]}-{years[-1]})")
    print(f"Excluded {total_years - len(years)} years with fewer than {MIN_VALID_DAYS} valid days")

    standardized = {}
    for label in labels:
        values = metrics[label]
        center = statistics.fmean(values)
        spread = statistics.pstdev(values)
        standardized[label] = [(value - center) / spread if spread else 0 for value in values]

    feature_count = len(labels)
    fig = plt.figure(figsize=(18, 17), facecolor="#fbfaf8")
    grid = fig.add_gridspec(
        feature_count, feature_count + 2,
        width_ratios=[1.5, *([1] * feature_count), 0.11],
        left=0.055, right=0.95, top=0.86, bottom=0.24,
        wspace=0.12, hspace=0.12,
    )

    violin_ax = fig.add_subplot(grid[:, 0])
    violin_values = [standardized[label] for label in labels]
    violins = violin_ax.violinplot(
        violin_values,
        positions=range(len(labels)),
        orientation="horizontal",
        showmedians=True,
        showextrema=False,
        widths=0.78,
    )
    violin_colors = [
        LIGHT_GREEN if index % 2 == 0 else LIGHT_PURPLE
        for index in range(feature_count)
    ]
    for index, body in enumerate(violins["bodies"]):
        body.set_facecolor(violin_colors[index])
        body.set_edgecolor(GREEN if index % 2 == 0 else PURPLE)
        body.set_linewidth(1.1)
        body.set_alpha(0.9)
    violins["cmedians"].set_color(INK)
    randomizer = random.Random(12)
    for index, values in enumerate(violin_values):
        jitter = [index + randomizer.uniform(-0.13, 0.13) for _ in values]
        violin_ax.scatter(values, jitter, s=7, color=PURPLE, alpha=0.23, linewidths=0)
    violin_ax.set_yticks(range(len(labels)), labels=labels)
    violin_ax.invert_yaxis()
    violin_ax.set_xlim(-3.5, 3.5)
    violin_ax.set_xlabel("Standardized annual value", color=INK)
    violin_ax.set_title("Annual distributions", loc="left", color=INK, pad=12, fontsize=11)
    violin_ax.grid(axis="x", color="#e6e1e5", linewidth=0.7)
    violin_ax.set_axisbelow(True)
    violin_ax.tick_params(colors=MUTED, labelsize=8)
    violin_ax.spines[["top", "right", "left"]].set_visible(False)
    violin_ax.spines["bottom"].set_color("#d8d2d7")

    axes = []
    limits = {label: padded_limits(metrics[label]) for label in labels}
    for row_index, row_label in enumerate(labels):
        axis_row = []
        for column_index, column_label in enumerate(labels):
            ax = fig.add_subplot(grid[row_index, column_index + 1])
            ax.set_facecolor("#fffefd")
            ax.tick_params(colors=MUTED, labelsize=7, length=2)
            for spine in ax.spines.values():
                spine.set_color("#e4dfe2")
                spine.set_linewidth(0.8)

            if row_index == column_index:
                values = metrics[row_label]
                ax.hist(values, bins=18, density=True, color=LIGHT_GREEN, alpha=0.7, edgecolor="none")
                density_x, density_y = kernel_density(values)
                if density_x:
                    ax.plot(density_x, density_y, color=PURPLE, linewidth=1.5)
                ax.set_xlim(limits[column_label])
                ax.set_yticks([])
            elif row_index < column_index:
                coefficient = spearman_correlation(metrics[column_label], metrics[row_label])
                ax.set_facecolor(CORRELATION_CMAP((coefficient + 1) / 2))
                ax.text(
                    0.5, 0.5, f"{coefficient:+.2f}",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=12, color=INK, fontweight="normal",
                )
                ax.set_xticks([])
                ax.set_yticks([])
            else:
                ax.scatter(
                    metrics[column_label], metrics[row_label],
                    s=15, color=GREEN, alpha=0.48, edgecolors="none",
                )
                ax.set_xlim(limits[column_label])
                ax.set_ylim(limits[row_label])
                ax.grid(color="#eee9eb", linewidth=0.55)
                ax.set_axisbelow(True)

            if row_index == 0:
                ax.set_title(column_label, fontsize=8, color=INK, pad=7)
            if column_index == 0:
                ax.set_ylabel(row_label, fontsize=8, color=INK)
            else:
                ax.set_yticklabels([])
            if row_index == len(labels) - 1:
                ax.set_xlabel(column_label, fontsize=7, color=INK, labelpad=4)
            else:
                ax.set_xticklabels([])
            axis_row.append(ax)
        axes.append(axis_row)

    colorbar_ax = fig.add_subplot(grid[:, -1])
    colorbar = fig.colorbar(
        ScalarMappable(norm=Normalize(-1, 1), cmap=CORRELATION_CMAP),
        cax=colorbar_ax,
    )
    colorbar.set_label("Spearman correlation", color=INK, fontsize=9)
    colorbar.ax.tick_params(colors=MUTED, labelsize=8)

    fig.suptitle("SUNSPOT ACTIVITY", x=0.07, y=0.96, ha="left", fontsize=19, color=INK, fontweight="bold")
    fig.text(
        0.07, 0.915,
        f"Annual metrics from {DATA.name} · {years[0]}–{years[-1]} · {len(years)} usable years",
        ha="left", fontsize=10, color=MUTED,
    )
    fig.text(
        0.055, 0.190,
        "Diagonal: distributions · upper: Spearman correlation · lower: one dot per year · left: standardized distributions",
        fontsize=8, color=INK,
    )
    fig.text(
        0.055, 0.155,
        "Metrics: mean, median, max, daily SD, IQR, 90th percentile, spotless-day %, and high-activity-day % (count >= 100).",
        fontsize=8, color=MUTED,
    )
    fig.text(
        0.055, 0.120,
        "Purple = negative correlation; green = positive. Daily -1 values and years with fewer than 300 valid days are excluded.",
        fontsize=8, color=MUTED,
    )

    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / PICTURE, dpi=180, facecolor=fig.get_facecolor())
    print(f"saved out/{PICTURE}")
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
