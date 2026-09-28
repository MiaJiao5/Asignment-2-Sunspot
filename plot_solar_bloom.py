# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib"]
# ///

r"""Turn the SILSO daily sunspot record into a chronological solar garden.

    uv run .\plot_solar_bloom.py
"""

import csv
import math
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Circle

FILE = "SN_d_tot_V2.0.csv"
PICTURE = "sunspot_solar_bloom.png"
MIN_YEAR_DAYS = 300
MIN_MONTH_DAYS = 15
HERE = Path(__file__).parent
DATA = HERE / "data" / FILE
OUT = HERE / "out"

PAPER = "#08090d"
INK = "#fff0c2"
MUTED = "#b99c6b"
STEM = "#78894c"
PETAL_EDGE = "#ffd978"
CORE = "#f3a51b"
SUNSET = LinearSegmentedColormap.from_list(
    "sunspot_sunset", ["#ffe27a", "#ffc13b", "#ff8b24", "#ed4d1c", "#9d2532"]
)


def read_data(path):
    by_month = defaultdict(list)
    by_year = defaultdict(list)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.reader(handle, delimiter=";"):
            if len(row) < 5:
                continue
            try:
                year, month = int(row[0].strip()), int(row[1].strip())
                count = float(row[4].strip())
            except ValueError:
                continue
            if count >= 0:
                by_month[(year, month)].append(count)
                by_year[year].append(count)

    years = sorted(year for year, counts in by_year.items() if len(counts) >= MIN_YEAR_DAYS)
    monthly = {}
    quiet_share = {}
    annual_mean = {}
    for year in years:
        year_counts = by_year[year]
        quiet_share[year] = sum(value == 0 for value in year_counts) / len(year_counts)
        annual_mean[year] = statistics.fmean(year_counts)
        for month in range(1, 13):
            counts = by_month.get((year, month), [])
            if len(counts) >= MIN_MONTH_DAYS:
                monthly[(year, month)] = (statistics.fmean(counts), statistics.pstdev(counts))
    return years, monthly, quiet_share, annual_mean, len(by_year)


def percentile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    weight = position - low
    return ordered[low] * (1 - weight) + ordered[high] * weight


def activity_peaks(years, monthly, minimum_gap=96):
    first_year, last_year = years[0], years[-1]
    series = [
        monthly.get((year, month), (None, None))[0]
        for year in range(first_year, last_year + 1)
        for month in range(1, 13)
    ]
    smoothed = [None] * len(series)
    weights = [0.5] + [1.0] * 11 + [0.5]
    for index in range(6, len(series) - 6):
        window = series[index - 6:index + 7]
        valid = [(value, weight) for value, weight in zip(window, weights) if value is not None]
        if len(valid) >= 11:
            smoothed[index] = sum(value * weight for value, weight in valid) / sum(
                weight for _, weight in valid
            )

    candidates = [
        index for index in range(1, len(smoothed) - 1)
        if smoothed[index] is not None
        and smoothed[index - 1] is not None
        and smoothed[index + 1] is not None
        and smoothed[index] >= smoothed[index - 1]
        and smoothed[index] > smoothed[index + 1]
    ]
    selected = []
    for index in sorted(candidates, key=lambda item: smoothed[item], reverse=True):
        if all(abs(index - other) >= minimum_gap for other in selected):
            selected.append(index)
    selected.sort()
    return [
        (first_year + index // 12, index % 12 + 1, smoothed[index], index)
        for index in selected
    ]


def petal(center_x, center_y, angle, length, half_width):
    dx, dy = math.cos(angle), math.sin(angle)
    nx, ny = -dy, dx
    upper, lower = [], []
    for step in range(13):
        t = step / 12
        fullness = math.sin(math.pi * t) ** 0.8
        x, y = center_x + length * t * dx, center_y + length * t * dy
        upper.append((x + half_width * fullness * nx, y + half_width * fullness * ny))
        lower.append((x - half_width * fullness * nx, y - half_width * fullness * ny))
    return upper + lower[::-1]


def main():
    years, monthly, quiet_share, annual_mean, all_years = read_data(DATA)
    if not years or not monthly:
        raise ValueError("No sufficiently complete sunspot observations were found")

    means = [value[0] for value in monthly.values()]
    spreads = [value[1] for value in monthly.values()]
    cycle_peaks = activity_peaks(years, monthly)
    if not cycle_peaks:
        raise ValueError("No smoothed solar activity peaks were detected")
    cycle_values = [peak[2] for peak in cycle_peaks]
    mean_cap = percentile(means + cycle_values, 0.95) or 1
    spread_cap = percentile(spreads, 0.95) or 1
    figure = plt.figure(figsize=(9, 16), facecolor=PAPER)
    garden = figure.add_axes([0.035, 0.20, 0.93, 0.66], facecolor=PAPER)
    garden.set(xlim=(-5.8, 5.8), ylim=(-7.55, 7.0), aspect="equal")
    garden.axis("off")

    center_y = 1.05
    garden.plot([0, 0.10, -0.04, 0.04, 0], [-1.3, -2.9, -4.4, -5.8, -7.0],
                color=STEM, linewidth=13, solid_capstyle="round", zorder=0)
    leaves = [
        petal(0.03, -3.35, 0.32, 1.95, 0.34),
        petal(-0.02, -4.65, 2.83, 1.82, 0.31),
    ]
    garden.add_collection(PolyCollection(
        leaves, facecolors=["#607c3d", "#789447"], edgecolors="#c7a94e",
        linewidths=1.0, closed=True, zorder=1,
    ))
    garden.add_collection(LineCollection(
        [((0.03, -3.35), (1.80, -2.74)), ((-0.02, -4.65), (-1.62, -3.75))],
        colors="#d5bd70", linewidths=1.2, zorder=2,
    ))

    ray_polygons, ray_values, ray_veins = [], [], []
    for index, (peak_year, peak_month, value, peak_slot) in enumerate(cycle_peaks):
        angle = math.pi / 2 - index * 2 * math.pi / len(cycle_peaks)
        dx, dy = math.cos(angle), math.sin(angle)
        root_x, root_y = 2.48 * dx, center_y + 2.48 * dy
        if len(cycle_peaks) > 1:
            next_peak = cycle_peaks[min(index + 1, len(cycle_peaks) - 1)][3]
            previous_peak = cycle_peaks[max(index - 1, 0)][3]
            period = (next_peak - previous_peak) / (2 if index not in (0, len(cycle_peaks) - 1) else 1)
        else:
            period = 132
        cycle_length = min(period / 132, 1.35)
        intensity = math.sqrt(min(value / mean_cap, 1))
        length = 0.72 + 1.10 * intensity
        width = 0.12 + 0.10 * math.sqrt(cycle_length)
        ray_polygons.append(petal(root_x, root_y, angle, length, width))
        ray_values.append(value)
        ray_veins.append(((root_x + length * 0.14 * dx, root_y + length * 0.14 * dy),
                          (root_x + length * 0.91 * dx, root_y + length * 0.91 * dy)))
        if index % 2 == 0:
            label_radius = 2.48 + length + 0.19
            garden.text(label_radius * dx, center_y + label_radius * dy, str(peak_year),
                        ha="center", va="center", color="#e8d6a8", fontsize=6, zorder=7)

    normalization = Normalize(0, mean_cap, clip=True)
    rays = PolyCollection(
        ray_polygons, cmap=SUNSET, norm=normalization,
        edgecolors=PETAL_EDGE, linewidths=0.9, closed=True, zorder=2,
    )
    rays.set_array(ray_values)
    garden.add_collection(rays)
    garden.add_collection(LineCollection(
        ray_veins, colors="#ffe6a1", linewidths=0.8, alpha=0.9, zorder=3
    ))

    garden.add_patch(Circle((0, center_y), 2.52, facecolor="#21120d",
                            edgecolor="#e77a25", linewidth=1.2, zorder=3.2))
    first_year = years[0]
    total_slots = (years[-1] - first_year + 1) * 12
    seed_x, seed_y, seed_means, seed_sizes = [], [], [], []
    golden_angle = math.pi * (3 - math.sqrt(5))
    seed_positions = {}
    for (year, month), (mean, spread) in sorted(monthly.items()):
        slot = (year - first_year) * 12 + month - 1
        radius = 2.42 * math.sqrt((slot + 0.5) / total_slots)
        angle = slot * golden_angle
        position = (radius * math.cos(angle), center_y + radius * math.sin(angle))
        seed_positions[slot] = position
        seed_x.append(position[0])
        seed_y.append(position[1])
        seed_means.append(mean)
        seed_sizes.append(2.0 + 8.0 * math.sqrt(min(spread / spread_cap, 1)))

    for stride, color, alpha in ((34, "#ffcf5c", 0.18), (55, "#ff773c", 0.15)):
        spiral_segments = [
            (position, seed_positions[slot + stride])
            for slot, position in seed_positions.items()
            if slot + stride in seed_positions
        ]
        garden.add_collection(LineCollection(
            spiral_segments, colors=color, linewidths=0.38, alpha=alpha, zorder=3.5
        ))

    seeds = garden.scatter(
        seed_x, seed_y, c=seed_means, s=seed_sizes,
        cmap=SUNSET, norm=normalization, edgecolors="none", alpha=0.96, zorder=4,
    )
    garden.scatter([0], [center_y], s=31, color="#ffe18a",
                   edgecolors="#fff3c8", linewidths=0.8, zorder=5)
    garden.plot([0, 0.10, -0.04, 0.04], [-1.42, -2.9, -4.4, -5.8],
                color="#bacb77", linewidth=1.0, alpha=0.75, zorder=5)

    figure.text(0.075, 0.955, "SOLAR BLOOM", color="#ffcf5c", fontsize=28,
                fontfamily="serif", fontweight="bold")
    figure.text(0.08, 0.929,
                f"A data-grown sunflower · {years[0]}—{years[-1]} · {len(monthly):,} observed months",
                color="#e5d2a7", fontsize=10)
    figure.text(0.08, 0.155, "THE SEED DISK", color="#ffcf5c", fontsize=8, fontweight="bold")
    figure.text(0.08, 0.132,
                "Each point is one month, placed chronologically on a golden-angle spiral; oldest near the center.",
                color="#e5d2a7", fontsize=8)
    figure.text(0.08, 0.103,
                "COLOR = monthly mean  ·  POINT SIZE = daily variability  ·  OUTER PETAL = 13-month-smoothed peak",
                color="#e5d2a7", fontsize=7.5)
    figure.text(0.08, 0.074,
                "Fine gold/coral curves trace Fibonacci seed spirals. Stem and leaves are ornamental; gaps are not zero.",
                color="#b99c6b", fontsize=7.5)
    figure.text(0.08, 0.045,
                f"SILSO · -1 excluded · sparse years omitted: {all_years - len(years)} · large early data gap: 1829–1848",
                color="#b99c6b", fontsize=7)

    figure.text(0.69, 0.155, "MEAN COUNT / DAY", color="#ffcf5c", fontsize=7, fontweight="bold")
    colorbar_ax = figure.add_axes([0.69, 0.132, 0.24, 0.012])
    colorbar = figure.colorbar(seeds, cax=colorbar_ax, orientation="horizontal", extend="max")
    colorbar.set_ticks([0, mean_cap / 2, mean_cap])
    colorbar.set_ticklabels(["0", f"{mean_cap / 2:.0f}", f"{mean_cap:.0f}+"])
    colorbar.ax.tick_params(labelsize=6, colors="#e5d2a7", length=2, pad=2)
    colorbar.outline.set_edgecolor("#79603d")

    OUT.mkdir(exist_ok=True)
    output = OUT / PICTURE
    figure.savefig(output, dpi=160, facecolor=PAPER)
    print(f"{FILE}: {len(years)} annual flowers ({years[0]}–{years[-1]})")
    print(f"{len(seed_means)} monthly seeds; {len(cycle_peaks)} smoothed solar peak petals")
    print(f"Excluded {all_years - len(years)} years with fewer than {MIN_YEAR_DAYS} valid days")
    print(f"saved {output.relative_to(HERE)}")
    plt.show()
    plt.close(figure)


if __name__ == "__main__":
    main()