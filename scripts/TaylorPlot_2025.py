# Taylor plot generator script
# Assumes stats for two models are already computed and provided in a structured form

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# Taylor plot generator script using axisartist FloatingSubplot
import numpy as np
import os
import pandas as pd
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
from matplotlib.projections import PolarAxes
from mpl_toolkits.axisartist import GridHelperCurveLinear

# ------------------------------
# User Parameters
# ------------------------------
# Provide the path to where stats CSVs are located
STATS_PATH = "../data/evaluation/"
PLOT_OUT_PATH = "../figs/"

# List model scenarios to include
SCENARIOS = ["FULLKEY_SC114_1", "FULLKEY_SC115_1", "FULLKEY_SC116_2"]
# Label and color mapping for each scenario
MODEL_CONFIGS = {
    "FULLKEY_SC114_1": {"label": "key run", "marker": "o", "edgecolor": "red"},
    "FULLKEY_SC115_1": {"label": "no wind", "marker": "s", "edgecolor": "black"},
    "FULLKEY_SC116_2": {"label": "PI curve", "marker": "v", "edgecolor": "black"},
}

# Mapping from reference dataset label to plot label and marker color
REFERENCE_CONFIGS = {
    "Allen": {"label": "vs. CO9", "color": "lightgreen"},
    "Suchy": {"label": "vs. satellite", "color": "black"},
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
    ax.legend(
        loc='upper left',  # Anchor position
        bbox_to_anchor=(0.75, 1.02),  # Shift legend box right and slightly up
        fontsize=7,
        ncols=1,
        frameon=True,
        borderaxespad=0.0
    )


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
        stats_file = os.path.join(STATS_PATH, f"bloom_timing_stats_{scenario}.csv")
        df = pd.read_csv(stats_file)

        for _, row in df.iterrows():
            source = row['Label']  # e.g., "Allen" or "Suchy"
            if source not in REFERENCE_CONFIGS:
                continue
            ref_config = REFERENCE_CONFIGS[source]
            model_stats.append({
                'model_stddev': row['Model StdDev'],
                'obs_stddev': row['Obs StdDev'],
                'corrcoef': row['R'],
                'label': MODEL_CONFIGS[scenario]['label'],
                'color': ref_config['color'],
                'edgecolor': MODEL_CONFIGS[scenario]['edgecolor'],
                'marker': MODEL_CONFIGS[scenario]['marker'],
                'source_label': ref_config['label']
            })

    plot_model_points(ax, refstd, model_stats)
    plt.tight_layout()
    plt.savefig(os.path.join(PLOT_OUT_PATH, "taylor_diagram_BloomTiming_c09Suchy.png"), dpi=300)
    print("fig saved")
    plt.show()

if __name__ == "__main__":
    main()
    print("done")