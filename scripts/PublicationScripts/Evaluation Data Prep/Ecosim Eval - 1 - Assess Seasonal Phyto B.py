# ============================================================================================
# Script:   Evaluate phytoplankton seasonal b versus model
# By: G Oldford, 2025
# Purpose:
#    - Plot seasonal b by group compared to outputs from McEwan
# Input:
#   - CSV of Ecosim output with season and date added
#   - seasonal phyto from McEwan's pub 10.1016/j.ecolmodel.2023.110402
#     Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\PhytoplanktonCommunity_McEwanSupp2023_\MODIFIED
# Output:
# - boxplots
#

# ============================================================================================
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')

# ================================
# Function definitions
# ================================
def load_and_filter_data(file_path, start_date, end_date):
    """Load CSV and filter by start and end date."""
    df = pd.read_csv(file_path, skiprows=0)
    df['date'] = pd.to_datetime(df['date'])
    mask = (df['date'] >= start_date) & (df['date'] <= end_date)
    return df.loc[mask]

def plot_boxplots_by_season(df, value_column, title_prefix=""):
    """Generate a boxplot of the specified value column grouped by season."""
    plt.figure(figsize=(10,6))
    df.boxplot(column=value_column, by='season')
    plt.title(f'{title_prefix} Boxplot of {value_column} by Season')
    plt.suptitle('')  # Removes default pandas boxplot title
    plt.xlabel('Season')
    plt.ylabel(value_column)
    plt.tight_layout()
    plt.show()


# ================================
# User inputs
# ================================
SCENARIO = "SC123"
START_DATE = '2015-01-01'
END_DATE = '2018-12-31'

# Path to CSV with date and season columns
OUTPUT_EXPORT_PATH = f"C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation/ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"

# Hardcoded McEwan observations (carbon units per m²)
obs_data = {
    'Diatoms': {'Winter': 0.65, 'Spring': 8.06, 'Summer': 1.34},
    'Nano': {'Winter': 1.57, 'Spring': 1.58, 'Summer': 1.53},
    'Other': {'Winter': 0.12, 'Spring': 0.19, 'Summer': 0.74}
}

# Model columns mapping
group_map = {
    'Diatoms': 17,
    'Nano': 18,
    'Other': 16
}

# === Load data ===
df = pd.read_csv(OUTPUT_EXPORT_PATH, skiprows=0)
df['date'] = pd.to_datetime(df['date'])
df = df[(df['date'] >= START_DATE) & (df['date'] <= END_DATE)]

# Calculate seasonal means
model_means = df.groupby('season')[[str(v) for v in group_map.values()]].mean()
model_means.rename(columns={str(v): k for k,v in group_map.items()}, inplace=True)

# Prepare McEwan data
obs_df = pd.DataFrame(obs_data).T  # Transpose so seasons become columns, groups as index

# === Build plot data ===
seasons = ['Winter', 'Spring', 'Summer']
groups = ['Diatoms', 'Nano', 'Other']

# Define group colors
group_colors = {
    'Diatoms': '#1f77b4',  # blue
    'Nano': '#ff7f0e',  # orange
    'Other': '#2ca02c'  # green
}

# Positions
x_labels = []
x_pos = []
model_heights = {g: [] for g in groups}
obs_heights = {g: [] for g in groups}

# Spacing
bar_width = 0.35
group_spacing = 1.0  # space between seasons
pair_spacing = 0.4  # space between model and obs within a season

pos = 0
for season in seasons:
    # Positions for model and obs within the season
    x_pos.append(pos - pair_spacing / 2)  # model
    x_pos.append(pos + pair_spacing / 2)  # obs
    x_labels.extend([f"Model-{season}", f"Obs-{season}"])

    # Heights
    for group in groups:
        model_heights[group].append(model_means.loc[season, group] if season in model_means.index else 0)
        obs_heights[group].append(obs_df.loc[group, season] if season in obs_df.columns else 0)

    pos += group_spacing

# === Plot ===
fig, ax = plt.subplots(figsize=(10, 6))

# Initialize bottoms
bottom_model = [0] * len(seasons)
bottom_obs = [0] * len(seasons)

# Plot model stacked bars (even indices)
for group in groups:
    values = model_heights[group]
    ax.bar(x_pos[::2], values, bar_width, bottom=bottom_model, color=group_colors[group], label=f'{group} (Model)')
    bottom_model = [i + j for i, j in zip(bottom_model, values)]

# Plot obs stacked bars (odd indices)
for group in groups:
    values = obs_heights[group]
    ax.bar(x_pos[1::2], values, bar_width, bottom=bottom_obs, color=group_colors[group], hatch='//', edgecolor='grey',
           label=f'{group} (Obs)')
    bottom_obs = [i + j for i, j in zip(bottom_obs, values)]

# === Final plot formatting ===
ax.set_xticks(x_pos)
ax.set_xticklabels(x_labels, rotation=45, ha='right')
ax.set_ylabel('Biomass (g C m$^{-2}$)')
ax.set_title(f'Model vs McEwan Seasonal Phytoplankton Biomass ({SCENARIO})')

# === Custom legend with only three entries ===
from matplotlib.patches import Patch
legend_patches = [
    Patch(facecolor=group_colors['Diatoms'], label='Diatoms'),
    Patch(facecolor=group_colors['Nano'], label='Nano'),
    Patch(facecolor=group_colors['Other'], label='Other')
]
ax.legend(handles=legend_patches, loc='upper right', fontsize=8)


plt.tight_layout()
plt.show()