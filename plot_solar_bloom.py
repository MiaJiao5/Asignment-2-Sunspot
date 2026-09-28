# /// script
# requires-python = ">=3.10"
# dependencies = ["matplotlib", "plotly"]
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
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_hex
from matplotlib.patches import Circle
import plotly.graph_objects as go

FILE = "SN_d_tot_V2.0.csv"
PICTURE = "sunspot_solar_bloom.png"
MIN_YEAR_DAYS = 300
MIN_MONTH_DAYS = 15
HERE = Path(__file__).parent
DATA = HERE / "data" / FILE
OUT = HERE / "out"
SITE = HERE / "site"

PAPER = "#f3ead8"
INK = "#342820"
MUTED = "#78685a"
PETAL_EDGE = "#b28b59"
CORE = "#bd742f"
SUNSET = LinearSegmentedColormap.from_list(
    "sunspot_sunset", ["#927454", "#c69b4a", "#e3a13b", "#cf5b32", "#812b43"]
)
QUIET_CMAP = LinearSegmentedColormap.from_list(
    "quiet_days", ["#c28a4e", "#83905e", "#45634e"]
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


def plotly_colorscale(colormap):
    return [[step / 8, to_hex(colormap(step / 8))] for step in range(9)]


def write_interactive_site(years, monthly, quiet_share, annual_mean, cycle_peaks,
                           mean_cap, spread_cap, annual_cap, quiet_cap):
    figure = go.Figure()
    normalization = Normalize(0, mean_cap, clip=True)
    first_year = years[0]
    total_slots = (years[-1] - first_year + 1) * 12
    seed_positions = {}
    seed_x, seed_y, seed_values, seed_sizes, seed_data = [], [], [], [], []

    for (year, month), (mean, spread) in sorted(monthly.items()):
        slot = (year - first_year) * 12 + month - 1
        radius = 2.70 * math.sqrt((slot + 0.5) / total_slots)
        angle = slot * math.pi * (3 - math.sqrt(5))
        position = (radius * math.cos(angle), radius * math.sin(angle))
        seed_positions[slot] = position
        seed_x.append(position[0])
        seed_y.append(position[1])
        seed_values.append(mean)
        seed_sizes.append(4 + 7 * math.sqrt(min(spread / spread_cap, 1)))
        seed_data.append([year, month, mean, spread])

    for radius, dash, alpha in ((2.88, "solid", 0.75), (3.66, "dot", 0.46), (5.06, "dot", 0.32)):
        figure.add_shape(
            type="circle", x0=-radius, y0=-radius, x1=radius, y1=radius,
            line={"color": "#a7865f", "width": 1, "dash": dash}, opacity=alpha,
        )

    figure.add_shape(
        type="circle", x0=-2.82, y0=-2.82, x1=2.82, y1=2.82,
        fillcolor="#e8ddc8", line={"color": "#9f7655", "width": 1}, layer="below",
    )
    for ring_radius in (0.82, 1.56, 2.28):
        figure.add_shape(
            type="circle", x0=-ring_radius, y0=-ring_radius,
            x1=ring_radius, y1=ring_radius,
            line={"color": "#8f6c48", "width": 0.6, "dash": "dot"}, opacity=0.32,
        )

    for stride, color, alpha in ((34, "#d79f3f", 0.25), (55, "#c65d37", 0.21)):
        spiral_x, spiral_y = [], []
        for slot, position in seed_positions.items():
            target = seed_positions.get(slot + stride)
            if target is not None:
                spiral_x.extend([position[0], target[0], None])
                spiral_y.extend([position[1], target[1], None])
        figure.add_trace(go.Scattergl(
            x=spiral_x, y=spiral_y, mode="lines", line={"color": color, "width": 0.55},
            opacity=alpha, hoverinfo="skip", showlegend=False,
        ))

    quiet_scale = plotly_colorscale(QUIET_CMAP)
    for index, year in enumerate(years):
        angle = math.pi / 2 - index * 2 * math.pi / len(years)
        start = 2.96
        end = start + 0.18 + 0.48 * math.sqrt(min(annual_mean[year] / annual_cap, 1))
        x0, y0 = start * math.cos(angle), start * math.sin(angle)
        x1, y1 = end * math.cos(angle), end * math.sin(angle)
        color = to_hex(QUIET_CMAP(quiet_share[year] / quiet_cap))
        width = 0.55 + 1.4 * quiet_share[year] / quiet_cap
        figure.add_shape(
            type="line", x0=x0, y0=y0, x1=x1, y1=y1,
            line={"color": color, "width": width},
        )

    years_x, years_y, years_data, years_sizes = [], [], [], []
    for index, year in enumerate(years):
        angle = math.pi / 2 - index * 2 * math.pi / len(years)
        radius = 2.96 + 0.18 + 0.48 * math.sqrt(min(annual_mean[year] / annual_cap, 1))
        years_x.append(radius * math.cos(angle))
        years_y.append(radius * math.sin(angle))
        years_data.append([year, annual_mean[year], quiet_share[year]])
        years_sizes.append(4 + 4 * math.sqrt(min(annual_mean[year] / annual_cap, 1)))

    figure.add_trace(go.Scattergl(
        x=years_x, y=years_y, mode="markers", name="Annual summaries",
        marker={"size": years_sizes, "color": [quiet_share[year] for year in years],
                "colorscale": quiet_scale, "cmin": 0, "cmax": quiet_cap,
                "showscale": False, "line": {"color": "#f3ead8", "width": 0.6}},
        customdata=years_data,
        hovertemplate="<b>%{customdata[0]}</b><br>Annual mean: %{customdata[1]:.1f}/day"
                      "<br>Spotless-day share: %{customdata[2]:.1%}<extra></extra>",
    ))

    peak_x, peak_y, peak_data, peak_values = [], [], [], []
    for index, (peak_year, peak_month, value, peak_slot) in enumerate(cycle_peaks):
        angle = math.pi / 2 - index * 2 * math.pi / len(cycle_peaks)
        if len(cycle_peaks) > 1:
            next_peak = cycle_peaks[min(index + 1, len(cycle_peaks) - 1)][3]
            previous_peak = cycle_peaks[max(index - 1, 0)][3]
            period = (next_peak - previous_peak) / (2 if index not in (0, len(cycle_peaks) - 1) else 1)
        else:
            period = 132
        cycle_length = min(period / 132, 1.35)
        intensity = math.sqrt(min(value / mean_cap, 1))
        length = 0.48 + 0.68 * intensity
        curve_x, curve_y = [], []
        for step in range(25):
            t = step / 24
            radius = 3.72 + length * t
            bend = (0.035 + 0.025 * cycle_length) * math.sin(math.pi * t) * (1 if index % 2 else -1)
            curve_x.append(radius * math.cos(angle + bend))
            curve_y.append(radius * math.sin(angle + bend))
        figure.add_trace(go.Scattergl(
            x=curve_x, y=curve_y, mode="lines",
            line={"color": to_hex(SUNSET(normalization(value)), keep_alpha=False),
                  "width": 0.6 + 1.2 * cycle_length / 1.35},
            hoverinfo="skip", showlegend=False,
        ))
        peak_x.append(curve_x[-1])
        peak_y.append(curve_y[-1])
        peak_data.append([peak_year, peak_month, value, period / 12])
        peak_values.append(value)

    figure.add_trace(go.Scattergl(
        x=seed_x, y=seed_y, mode="markers", name="Monthly observations",
        marker={"size": seed_sizes, "color": seed_values,
                "colorscale": plotly_colorscale(SUNSET), "cmin": 0, "cmax": mean_cap,
                "showscale": False, "opacity": 0.95,
                "line": {"color": "#f7e9cf", "width": 0.2}},
        customdata=seed_data,
        hovertemplate="<b>%{customdata[0]}-%{customdata[1]:02d}</b>"
                      "<br>Monthly mean: %{customdata[2]:.1f}/day"
                      "<br>Daily SD: %{customdata[3]:.1f}<extra></extra>",
    ))
    figure.add_trace(go.Scattergl(
        x=peak_x, y=peak_y, mode="markers", name="Smoothed activity peaks",
        marker={"size": 10, "color": peak_values,
                "colorscale": plotly_colorscale(SUNSET), "cmin": 0, "cmax": mean_cap,
                "showscale": False, "symbol": "circle-open",
                "line": {"color": "#5a382b", "width": 1.2}},
        customdata=peak_data,
        hovertemplate="<b>Peak %{customdata[0]}-%{customdata[1]:02d}</b>"
                      "<br>13-month smoothed mean: %{customdata[2]:.1f}/day"
                      "<br>Spacing estimate: %{customdata[3]:.1f} years<extra></extra>",
    ))

    figure.add_trace(go.Scattergl(
        x=[0], y=[0], mode="markers+text", text=[str(years[0])],
        textposition="middle center", textfont={"size": 10, "color": "#342820"},
        marker={"size": 18, "color": "#f3ead8", "line": {"color": "#bd742f", "width": 1.5}},
        name="Earliest year", hovertemplate=f"Earliest usable year: {years[0]}<extra></extra>",
    ))

    figure.update_layout(
        height=1000, autosize=True,
        paper_bgcolor="#f3ead8", plot_bgcolor="#f3ead8",
        margin={"l": 15, "r": 15, "t": 20, "b": 45},
        font={"family": "Georgia, serif", "color": "#342820"},
        hoverlabel={"bgcolor": "#342820", "font": {"color": "#fff0c2"}},
        hovermode="closest", dragmode="pan", uirevision="sunspot-mandala",
        legend={"orientation": "h", "y": -0.04, "x": 0.5, "xanchor": "center"},
        xaxis={"range": [-5.5, 5.5], "visible": False, "constrain": "domain"},
        yaxis={"range": [-5.5, 5.5], "visible": False, "scaleanchor": "x", "scaleratio": 1},
    )

    chart = figure.to_html(
        full_html=False, include_plotlyjs=True,
        config={"responsive": True, "scrollZoom": True, "displaylogo": False,
                "modeBarButtonsToRemove": ["lasso2d", "select2d"]},
    )
    peak_years = ", ".join(str(peak[0]) for peak in cycle_peaks)
    html = f'''<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <meta name="description" content="An interactive, data-driven sunspot phyllotaxis mandala built from SILSO daily observations.">
  <title>Sunspot Phyllotaxis | SILSO</title>
  <style>
    :root {{ color-scheme: light; --paper: #f3ead8; --ink: #342820; --muted: #78685a; --accent: #bd742f; }}
    * {{ box-sizing: border-box; }}
    body {{ margin: 0; background: var(--paper); color: var(--ink); font-family: Georgia, "Palatino Linotype", serif; }}
    main {{ width: min(100% - 36px, 1120px); margin: 0 auto; padding: 32px 0 36px; }}
    header {{ border-bottom: 1px solid #cbb99d; padding: 0 4px 18px; }}
    .eyebrow {{ margin: 0 0 8px; color: var(--accent); font: 700 11px/1.4 system-ui, sans-serif; }}
    h1 {{ margin: 0; font-size: clamp(30px, 5vw, 52px); font-weight: 500; line-height: 1.05; }}
    .subtitle {{ margin: 10px 0 0; color: var(--muted); font: 14px/1.5 system-ui, sans-serif; }}
    .stats {{ display: flex; flex-wrap: wrap; gap: 8px 24px; margin-top: 15px; color: var(--ink); font: 12px/1.5 system-ui, sans-serif; }}
    .stats strong {{ color: var(--accent); }}
    .motion-controls {{ display: flex; align-items: center; flex-wrap: wrap; gap: 10px 14px; margin: 14px 4px 0; color: var(--muted); font: 12px/1.4 system-ui, sans-serif; }}
    .motion-controls button {{ border: 1px solid #a7865f; border-radius: 3px; padding: 7px 11px; background: transparent; color: var(--ink); font: inherit; cursor: pointer; }}
    .motion-controls button:hover {{ background: #e8ddc8; }}
    .motion-controls input {{ width: 130px; accent-color: var(--accent); vertical-align: middle; }}
    .motion-controls output {{ min-width: 3.5em; color: var(--ink); }}
    .chart {{ margin: 8px auto 0; width: 100%; }}
    .guide {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px; border-top: 1px solid #cbb99d; padding: 18px 4px 0; }}
    .guide h2 {{ margin: 0 0 7px; color: var(--accent); font: 700 11px/1.4 system-ui, sans-serif; }}
    .guide p {{ margin: 0; color: var(--muted); font: 12px/1.6 system-ui, sans-serif; }}
    .swatch {{ height: 7px; margin: 9px 0 5px; border-radius: 10px; }}
    .sun-scale {{ background: linear-gradient(90deg, #927454, #c69b4a, #e3a13b, #cf5b32, #812b43); }}
    .quiet-scale {{ background: linear-gradient(90deg, #c28a4e, #83905e, #45634e); }}
    footer {{ margin: 20px 4px 0; color: var(--muted); font: 11px/1.6 system-ui, sans-serif; }}
    @media (max-width: 680px) {{ main {{ width: min(100% - 20px, 560px); padding-top: 20px; }} .guide {{ grid-template-columns: 1fr; gap: 16px; }} }}
  </style>
</head>
<body>
  <main>
    <header>
      <p class="eyebrow">SILSO · SOLAR OBSERVATORY</p>
      <h1>Sunspot Phyllotaxis</h1>
      <p class="subtitle">A living geometry drawn from daily sunspot observations, {years[0]}–{years[-1]}.</p>
      <div class="stats"><span><strong>{len(monthly):,}</strong> monthly records</span><span><strong>{len(years)}</strong> usable years</span><span><strong>{len(cycle_peaks)}</strong> smoothed peaks</span></div>
    </header>
        <div class="motion-controls" aria-label="Visualization rotation controls">
            <button id="rotation-toggle" type="button" aria-pressed="true">Pause rotation</button>
            <label for="rotation-speed">Rotation speed</label>
            <input id="rotation-speed" type="range" min="0.5" max="10" step="0.5" value="6">
            <output id="rotation-speed-value" for="rotation-speed">6.0°/s</output>
        </div>
    <section class="chart" aria-label="Interactive sunspot data visualization">{chart}</section>
        <script>
            (() => {{
                const graph = document.querySelector('.js-plotly-plot');
                const toggle = document.getElementById('rotation-toggle');
                const speed = document.getElementById('rotation-speed');
                const speedValue = document.getElementById('rotation-speed-value');
                const reduceMotion = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
                const originalTraces = graph.data.map(trace => ({{
                    x: Array.from(trace.x || []),
                    y: Array.from(trace.y || [])
                }}));
                const originalShapes = (graph.layout.shapes || []).map(shape => ({{ ...shape }}));
                const traceIndices = graph.data.map((_, index) => index);
                let rotating = !reduceMotion;
                let angle = 0;
                let lastUpdate = 0;
                let updating = false;
                let animationTimer;
                const frameDelay = 120;

                function rotatePoint(x, y, cosine, sine) {{
                    return [x * cosine - y * sine, x * sine + y * cosine];
                }}

                function renderRotation() {{
                    const cosine = Math.cos(angle);
                    const sine = Math.sin(angle);
                    const traces = originalTraces.map(trace => {{
                        const x = new Array(trace.x.length);
                        const y = new Array(trace.y.length);
                        for (let index = 0; index < x.length; index += 1) {{
                            const originalX = trace.x[index];
                            const originalY = trace.y[index];
                            if (Number.isFinite(originalX) && Number.isFinite(originalY)) {{
                                x[index] = originalX * cosine - originalY * sine;
                                y[index] = originalX * sine + originalY * cosine;
                            }} else {{
                                x[index] = originalX;
                                y[index] = originalY;
                            }}
                        }}
                        return {{ x, y }};
                    }});
                    const shapes = originalShapes.map(shape => {{
                        if (shape.type === 'circle') return shape;
                        const start = rotatePoint(shape.x0, shape.y0, cosine, sine);
                        const end = rotatePoint(shape.x1, shape.y1, cosine, sine);
                        return {{ ...shape, x0: start[0], y0: start[1], x1: end[0], y1: end[1] }};
                    }});
                    updating = true;
                    Plotly.update(
                        graph,
                        {{ x: traces.map(trace => trace.x), y: traces.map(trace => trace.y) }},
                        {{ shapes }},
                        traceIndices
                    ).finally(() => {{ updating = false; }});
                }}

                function animate() {{
                    if (!rotating) return;
                    const time = performance.now();
                    const elapsed = lastUpdate ? Math.min(time - lastUpdate, 250) : frameDelay;
                    lastUpdate = time;
                    angle += Number(speed.value) * Math.PI / 180 * elapsed / 1000;
                    if (!updating) {{
                        renderRotation();
                    }}
                    animationTimer = window.setTimeout(animate, frameDelay);
                }}

                toggle.setAttribute('aria-pressed', String(rotating));
                toggle.textContent = rotating ? 'Pause rotation' : 'Start rotation';
                toggle.addEventListener('click', () => {{
                    rotating = !rotating;
                    toggle.setAttribute('aria-pressed', String(rotating));
                    toggle.textContent = rotating ? 'Pause rotation' : 'Start rotation';
                    lastUpdate = 0;
                    if (rotating) {{
                        animate();
                    }} else {{
                        window.clearTimeout(animationTimer);
                        animationTimer = undefined;
                    }}
                }});
                speed.addEventListener('input', () => {{
                    speedValue.value = `${{Number(speed.value).toFixed(1)}}°/s`;
                }});
                if (rotating) animationTimer = window.setTimeout(animate, frameDelay);
            }})();
        </script>
    <section class="guide">
      <div><h2>MONTHLY OBSERVATIONS</h2><p>Each seed is one month on a golden-angle spiral. Hover for its mean and daily variation; use the legend to toggle data layers.</p><div class="swatch sun-scale"></div><p>Monthly mean sunspot count per day</p></div>
      <div><h2>ANNUAL RING</h2><p>Tick length and radial distance show annual mean; green hue shows the share of spotless days. Outer rays mark 13-month-smoothed activity peaks.</p><div class="swatch quiet-scale"></div><p>Annual spotless-day share</p></div>
    </section>
    <footer>Source: SILSO daily total sunspot numbers · Missing values (-1) are excluded. Months with fewer than {MIN_MONTH_DAYS} valid observations and years with fewer than {MIN_YEAR_DAYS} valid days are omitted. The record has a substantial observation gap from 1829 to 1848. Detected peak years: {peak_years}.</footer>
  </main>
</body>
</html>'''
    SITE.mkdir(exist_ok=True)
    page = SITE / "index.html"
    page.write_text(html, encoding="utf-8")
    return page


def main():
    years, monthly, quiet_share, annual_mean, all_years = read_data(DATA)
    if not years or not monthly:
        raise ValueError("No sufficiently complete sunspot observations were found")

    means = [value[0] for value in monthly.values()]
    cycle_peaks = activity_peaks(years, monthly)
    if not cycle_peaks:
        raise ValueError("No smoothed solar activity peaks were detected")
    cycle_values = [peak[2] for peak in cycle_peaks]
    mean_cap = percentile(means + cycle_values, 0.95) or 1
    spread_cap = percentile([value[1] for value in monthly.values()], 0.95) or 1
    annual_cap = percentile(list(annual_mean.values()), 0.95) or 1
    quiet_cap = max(quiet_share.values()) or 1
    figure = plt.figure(figsize=(9, 12), facecolor=PAPER)
    garden = figure.add_axes([0.10, 0.235, 0.80, 0.58], facecolor=PAPER)
    garden.set(xlim=(-5.5, 5.5), ylim=(-5.5, 5.5), aspect="equal")
    garden.axis("off")

    normalization = Normalize(0, mean_cap, clip=True)
    for radius, alpha, style in ((2.88, 0.75, "solid"), (3.66, 0.46, (0, (2, 4))), (5.06, 0.32, (0, (1, 5)))):
        garden.add_patch(Circle(
            (0, 0), radius, fill=False, edgecolor="#a7865f",
            linewidth=0.55, alpha=alpha, linestyle=style, zorder=1,
        ))

    ray_paths, ray_values, ray_widths, ray_tips = [], [], [], []
    for index, (peak_year, peak_month, value, peak_slot) in enumerate(cycle_peaks):
        angle = math.pi / 2 - index * 2 * math.pi / len(cycle_peaks)
        dx, dy = math.cos(angle), math.sin(angle)
        if len(cycle_peaks) > 1:
            next_peak = cycle_peaks[min(index + 1, len(cycle_peaks) - 1)][3]
            previous_peak = cycle_peaks[max(index - 1, 0)][3]
            period = (next_peak - previous_peak) / (2 if index not in (0, len(cycle_peaks) - 1) else 1)
        else:
            period = 132
        cycle_length = min(period / 132, 1.35)
        intensity = math.sqrt(min(value / mean_cap, 1))
        length = 0.48 + 0.68 * intensity
        curve = []
        for step in range(25):
            t = step / 24
            radius = 3.72 + length * t
            bend = (0.035 + 0.025 * cycle_length) * math.sin(math.pi * t) * (1 if index % 2 else -1)
            curve.append((radius * math.cos(angle + bend), radius * math.sin(angle + bend)))
        ray_paths.append(curve)
        ray_values.append(value)
        ray_widths.append(0.6 + 1.2 * cycle_length / 1.35)
        ray_tips.append(curve[-1])
        if index % 2 == 0:
            label_radius = 5.22
            garden.text(label_radius * dx, label_radius * dy, str(peak_year),
                        ha="center", va="center", color=MUTED, fontsize=6, zorder=7)

    rays = LineCollection(
        ray_paths, cmap=SUNSET, norm=normalization,
        linewidths=ray_widths, alpha=0.86, capstyle="round", zorder=2,
    )
    rays.set_array(ray_values)
    garden.add_collection(rays)
    garden.scatter([point[0] for point in ray_tips], [point[1] for point in ray_tips],
                   c=ray_values, s=18, cmap=SUNSET, norm=normalization,
                   edgecolors=INK, linewidths=0.45, zorder=3)

    garden.add_patch(Circle((0, 0), 2.82, facecolor="#e8ddc8",
                            edgecolor="#9f7655", linewidth=0.7, zorder=3.1))
    first_year = years[0]
    total_slots = (years[-1] - first_year + 1) * 12
    seed_x, seed_y, seed_means, seed_sizes = [], [], [], []
    golden_angle = math.pi * (3 - math.sqrt(5))
    seed_positions = {}
    for (year, month), (mean, spread) in sorted(monthly.items()):
        slot = (year - first_year) * 12 + month - 1
        radius = 2.70 * math.sqrt((slot + 0.5) / total_slots)
        angle = slot * golden_angle
        position = (radius * math.cos(angle), radius * math.sin(angle))
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
            spiral_segments, colors=color, linewidths=0.42, alpha=alpha, zorder=3.4
        ))

    seeds = garden.scatter(
        seed_x, seed_y, c=seed_means, s=seed_sizes,
        cmap=SUNSET, norm=normalization, edgecolors="#f7e9cf", linewidths=0.12,
        alpha=0.98, zorder=4,
    )

    annual_segments = []
    annual_quiet = []
    annual_widths = []
    for index, year in enumerate(years):
        angle = math.pi / 2 - index * 2 * math.pi / len(years)
        radius_start = 2.96
        radius_end = radius_start + 0.18 + 0.48 * math.sqrt(min(annual_mean[year] / annual_cap, 1))
        annual_segments.append((
            (radius_start * math.cos(angle), radius_start * math.sin(angle)),
            (radius_end * math.cos(angle), radius_end * math.sin(angle)),
        ))
        annual_quiet.append(quiet_share[year])
        annual_widths.append(0.55 + 1.4 * quiet_share[year] / quiet_cap)
    annual_ring = LineCollection(
        annual_segments, cmap=QUIET_CMAP, norm=Normalize(0, quiet_cap),
        linewidths=annual_widths, capstyle="round", zorder=5,
    )
    annual_ring.set_array(annual_quiet)
    garden.add_collection(annual_ring)

    for ring_radius in (0.82, 1.56, 2.28):
        garden.add_patch(Circle((0, 0), ring_radius, fill=False,
                                edgecolor="#8f6c48", linewidth=0.35,
                                alpha=0.32, linestyle=(0, (1, 4)), zorder=5))

    garden.add_patch(Circle((0, 0), 0.23, facecolor=PAPER,
                            edgecolor=CORE, linewidth=0.9, zorder=6))
    garden.add_patch(Circle((0, 0), 0.16, fill=False,
                            edgecolor=PETAL_EDGE, linewidth=0.55, zorder=7))
    garden.text(0, 0, str(years[0]), ha="center", va="center",
                color=INK, fontsize=5.5, fontweight="bold", zorder=8)

    for tick_index in range(72):
        angle = math.pi / 2 - tick_index * 2 * math.pi / 72
        inner = 5.00 if tick_index % 6 == 0 else 5.025
        outer = 5.12
        garden.plot([inner * math.cos(angle), outer * math.cos(angle)],
                    [inner * math.sin(angle), outer * math.sin(angle)],
                    color=MUTED, linewidth=0.38 if tick_index % 6 else 0.75,
                    alpha=0.58, zorder=6)

    figure.text(0.075, 0.955, "SUNSPOT / PHYLLOTAXIS", color=INK, fontsize=22,
                fontfamily="serif", fontweight="bold")
    figure.text(0.08, 0.929,
                f"A solar mandala built from {len(monthly):,} monthly observations · {years[0]}—{years[-1]}",
                color=MUTED, fontsize=9)
    figure.text(0.08, 0.205, "THREE DATA LAYERS", color=CORE, fontsize=8, fontweight="bold")
    figure.text(0.08, 0.178,
                "Each seed is a month; time spirals outward from the center. Fine arcs trace 34- and 55-step Fibonacci families.",
                color=INK, fontsize=7.1)
    figure.text(0.08, 0.15,
                "Seed color = monthly mean · size = daily variation · annual tick length = yearly mean · green hue = spotless-day share.",
                color=INK, fontsize=7.1)
    figure.text(0.08, 0.122,
                "Outer rays = 13-month-smoothed activity peaks; ray width reflects peak spacing.",
                color=INK, fontsize=7.1)
    figure.text(0.08, 0.094,
                f"SILSO daily totals · -1 excluded · {all_years - len(years)} sparse years omitted · record gap: 1829–1848",
                color=MUTED, fontsize=6.8)

    figure.text(0.08, 0.061, "MONTHLY MEAN", color=CORE, fontsize=6.5, fontweight="bold")
    monthly_bar_ax = figure.add_axes([0.08, 0.041, 0.34, 0.012])
    monthly_bar = figure.colorbar(seeds, cax=monthly_bar_ax, orientation="horizontal", extend="max")
    monthly_bar.set_ticks([0, mean_cap / 2, mean_cap])
    monthly_bar.set_ticklabels(["0", f"{mean_cap / 2:.0f}", f"{mean_cap:.0f}+"])
    monthly_bar.ax.tick_params(labelsize=5.5, colors=INK, length=2, pad=1)
    monthly_bar.outline.set_edgecolor("#a7865f")

    figure.text(0.56, 0.061, "ANNUAL SPOTLESS-DAY SHARE", color=CORE, fontsize=6.5, fontweight="bold")
    quiet_bar_ax = figure.add_axes([0.56, 0.041, 0.34, 0.012])
    quiet_bar = figure.colorbar(annual_ring, cax=quiet_bar_ax, orientation="horizontal")
    quiet_bar.set_ticks([0, quiet_cap])
    quiet_bar.set_ticklabels(["0%", f"{quiet_cap:.0%}"])
    quiet_bar.ax.tick_params(labelsize=5.5, colors=INK, length=2, pad=1)
    quiet_bar.outline.set_edgecolor("#a7865f")

    OUT.mkdir(exist_ok=True)
    output = OUT / PICTURE
    figure.savefig(output, dpi=160, facecolor=PAPER)
    page = write_interactive_site(
        years, monthly, quiet_share, annual_mean, cycle_peaks,
        mean_cap, spread_cap, annual_cap, quiet_cap,
    )
    print(f"{FILE}: {len(years)} usable years ({years[0]}–{years[-1]})")
    print(f"{len(seed_means)} monthly seeds; {len(cycle_peaks)} smoothed solar peak petals")
    print(f"Excluded {all_years - len(years)} years with fewer than {MIN_YEAR_DAYS} valid days")
    print(f"saved {output.relative_to(HERE)}")
    print(f"saved {page.relative_to(HERE)}")
    plt.show()
    plt.close(figure)


if __name__ == "__main__":
    main()