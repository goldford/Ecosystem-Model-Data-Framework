# Taylor plot generator script
# Assumes stats for two models are already computed and provided in a structured form

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
import ecospace_eval_config as cfg

# Taylor plot generator script using axisartist FloatingSubplot
import numpy as np
import glob
import re
import os
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import floating_axes as FA
from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter
from mpl_toolkits.axisartist import Subplot
from mpl_toolkits.axisartist import angle_helper
from mpl_toolkits.axisartist import floating_axes
from mpl_toolkits.axisartist import grid_finder as GF
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from mpl_toolkits.axisartist.grid_finder import MaxNLocator
from mpl_toolkits.axisartist import SubplotHost
from matplotlib.projections.polar import PolarAxes # for older matplotlib
# from matplotlib.projections import PolarAxes # for newer matplotlib
from mpl_toolkits.axisartist import GridHelperCurveLinear
import ecospace_eval_config as cfg

# ------------------------------
# User Parameters
# ------------------------------
# Provide the path to where stats CSVs are located
STATS_PATH = Path(cfg.EVALOUT_P)
PLOT_OUT_PATH = Path(cfg.FIGS_P)
FORMAT = "ecospace"   # options: "ecospace" or "ecosim"

# Define patterns & extractors per format
FILE_PATTERNS = {
    "ecospace": cfg.TY_STATS_FPAT_SPC_FMT,
    "ecosim":   cfg.TY_STATS_FPAT_SIM_FMT,
}
SCENARIO_RE = {
    "ecospace": re.compile(cfg.TY_STATS_F_SPC_FMT),
    "ecosim":   re.compile(cfg.TY_STATS_F_SIM_FMT),
}

pattern = FILE_PATTERNS[FORMAT]
rx = SCENARIO_RE[FORMAT]

# Find all matching files
stats_files = glob.glob(os.path.join(STATS_PATH, pattern))

# Extract scenario names from filenames
SCENARIOS = []
for f in stats_files:
    m = rx.match(Path(f).name)
    if m:
        SCENARIOS.append(m.group(1))

# Auto-build default MODEL_CONFIGS (with optional overrides)
MODEL_CONFIGS = {s: {"label": s, "marker": "o", "edgecolor": "black"} for s in SCENARIOS}

# Only specify exceptions you care about
OVERRIDES = {
    "SC139": {"label": "139", "edgecolor": "red"},
    # "SC150": {"label": "150*", "edgecolor": "blue"},
}
for scen, o in OVERRIDES.items():
    if scen in MODEL_CONFIGS:
        MODEL_CONFIGS[scen].update(o)

# Main loop now reads discovered files instead of formatting names
model_stats = []
for scenario in SCENARIOS:
    if FORMAT == "ecospace":
        stats_file = os.path.join(STATS_PATH, f"ecospace_bloom_timing_stats_{scenario}.csv")
    else:
        stats_file = os.path.join(STATS_PATH, f"ecosim_bloom_eval_stats_{scenario}.csv")

    df = pd.read_csv(stats_file)


# Mapping from reference dataset label to plot label and marker color
REFERENCE_CONFIGS = {
    "C09": {"label": "vs. CO9", "color": "lightgreen"},
    "sat": {"label": "vs. satellite", "color": "black"},
}


# ------------------------------
# Taylor Diagram Functions
# ------------------------------

def get_taylor_diagram_axes(fig, rect, refstd, srange, contour_levs, clr_std_cntr, extend=False):
    grid_clr = '#969696'
    cntr_clr = '#858585'
    cntr_lw = 0.8
    clbl_sz = 10
    clbl_clr = '#858585'
    axs_lbl_fs = 10
    str_sz = 10
    star_offset = 0.01

    rlocs = np.array([0, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1])
    tr = PolarAxes.PolarTransform()
    tmax = np.pi if extend else np.pi / 2
    if extend:
        rlocs = np.concatenate((-rlocs[:0:-1], rlocs))
    tlocs = np.arccos(rlocs)
    gl1 = FixedLocator(tlocs)
    tf1 = DictFormatter(dict(zip(tlocs, map(str, rlocs))))

    smin = srange[0] * refstd
    smax = srange[1] * refstd
    ghelper = FA.GridHelperCurveLinear(tr, extremes=(0, tmax, smin, smax), grid_locator1=gl1, tick_formatter1=tf1)
    ax = FA.FloatingSubplot(fig, rect, grid_helper=ghelper)
    fig.add_subplot(ax)

    ax.axis["top"].set_axis_direction("bottom")
    ax.axis["top"].toggle(ticklabels=True, label=True)
    ax.axis["top"].major_ticklabels.set_axis_direction("top")
    ax.axis["top"].label.set_axis_direction("top")
    ax.axis["top"].label.set_text("Correlation")
    ax.axis["left"].set_axis_direction("bottom")
    ax.axis["left"].label.set_text("Standard deviation")
    ax.axis["right"].set_axis_direction("top")
    ax.axis["right"].toggle(ticklabels=True)
    ax.axis["right"].major_ticklabels.set_axis_direction("left")
    ax.axis["bottom"].set_visible(False)

    _ax = ax
    ax = ax.get_aux_axes(tr)
    ax.plot([0 + star_offset], refstd, 'k*', ls='', ms=str_sz, label="Reference", zorder=10)
    t = np.linspace(0, tmax)
    r = np.zeros_like(t) + refstd
    ax.plot(t, r, '--', label='_', color=clr_std_cntr)
    rs, ts = np.meshgrid(np.linspace(smin, smax), np.linspace(0, tmax))
    rms = np.sqrt(refstd**2 + rs**2 - 2 * refstd * rs * np.cos(ts))
    _ax.grid(ls='-', color=grid_clr, lw=0.8)
    contours = ax.contour(ts, rs, rms, colors=cntr_clr, levels=contour_levs,
                          linewidths=cntr_lw, linestyles='--')
    plt.clabel(contours, inline=1, fontsize=clbl_sz, fmt='%.1f', colors=clbl_clr)
    return ax


def plot_model_points(ax, refstd, model_stats):
    # Sort so that all 'CO9' (Allen) entries are drawn first, then 'satellite' (Suchy)
    model_stats = sorted(model_stats, key=lambda x: x['source_label'])
    for entry in model_stats:
        theta = np.arccos(entry['corrcoef'])
        std = entry['model_stddev'] / entry['obs_stddev']  # Normalize
        ax.plot(theta, std, marker=entry['marker'], color=entry['color'],
                markeredgecolor=entry['edgecolor'], markersize=8,
                linestyle='None',
                label=f"{entry['label']} ({entry['source_label']})")
    # ax.legend(
    #     loc='upper left',  # Anchor position
    #     bbox_to_anchor=(0.75, 1.02),  # Shift legend box right and slightly up
    #     fontsize=7,
    #     ncols=1,
    #     frameon=True,
    #     borderaxespad=0.0
    # )


# ------------------------------
# Main Execution
# ------------------------------
def main():
    fig = plt.figure(figsize=(6, 5))
    refstd = 1.0
    srange = (0, 1.5)
    contour_levs = [0.5, 1.0, 1.5, 2.0]
    ax = get_taylor_diagram_axes(fig, 111, refstd, srange, contour_levs, clr_std_cntr='grey')

    model_stats = []
    for scenario in SCENARIOS:
        if FORMAT == "ecospace":
            stats_file = os.path.join(STATS_PATH, f"ecospace_bloom_timing_stats_{scenario}.csv")
        else:
            stats_file = os.path.join(STATS_PATH, f"ecosim_bloom_eval_stats_{scenario}.csv")

        df = pd.read_csv(stats_file)
        df.columns = df.columns.str.lower()
        if "label" not in df.columns:
            # assume first row is sat, second is C09
            default_labels = ["sat", "C09"]
            df["label"] = default_labels[:len(df)]

        for _, row in df.iterrows():
            source = row['label']

            if source not in REFERENCE_CONFIGS:
                continue
            ref_config = REFERENCE_CONFIGS[source]
            print(MODEL_CONFIGS[scenario]['label'])
            print(source)
            print(row['r'])
            model_stats.append({
                'model_stddev': row['model stddev'],
                'obs_stddev': row['obs stddev'],
                'corrcoef': row['r'],
                'label': MODEL_CONFIGS[scenario]['label'],
                'color': ref_config['color'],
                'edgecolor': MODEL_CONFIGS[scenario]['edgecolor'],
                'marker': MODEL_CONFIGS[scenario]['marker'],
                'source_label': ref_config['label']
            })

    plot_model_points(ax, refstd, model_stats)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_OUT_PATH, "ecospace_taylor_diagram_BloomTiming.png"), dpi=300)
    print("fig saved")
    plt.show()

if __name__ == "__main__":
    main()
    print("done")