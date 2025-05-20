# ============================================================================================
# Script:    Ecospace Analysis - 2b Nemcek vs Ecospace SSC Phyto v2.py
# Created:   July 11, 2024
# Revised:   May, 2025
# Author:    G. Oldford
#
# Purpose:
#    Compare observed phytoplankton data from Nemcek et al. (2023) with model outputs from:
#      - Ecospace (spatial ecosystem model)
#      - SalishSeaCast (3D physical-biogeochemical model)
#
# Inputs:
#    - CSV file output from prior matching script:
#        "Nemcek_matched_to_model_out_<ecospace_code>.csv"
#      containing matched observations and corresponding model estimates.
#
# Outputs:
#    - Summary bar plots of sample count by month and season per subdomain
#    - Seasonal average and anomaly (z-score) plots by functional group and subdomain
#    - PNG figures saved for visual comparisons of functional group trends
#
# Data Sources:
#    - Observations: Nemcek, N., Hennekes, M., Sastri, A., & Perry, R. I. (2023).
#         Seasonal and spatial dynamics of the phytoplankton community in the Salish Sea, 2015–2019.
#         *Progress in Oceanography*, 217, 103108.
#         https://doi.org/10.1016/j.pocean.2023.103108
#
#    - Model: SalishSeaCast ERDDAP download (see: https://salishsea.eos.ubc.ca/erddap/)
#
# Notes:
#    - Requires prior script to generate the matched CSV file
#    - Includes a data reduction step to prevent duplicate influence of identical timestamps
# ============================================================================================


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Configuration and Input Setup
# -------------------------------

ecospace_code = "Scv89_1"  # Choose scenario code that matches pre-processed dataset
input_dir = "C://Users//Greig/Documents//github//Ecosystem-Model-Data-Framework//data//evaluation//"
input_file = f"Nemcek_matched_to_model_out_{ecospace_code}.csv"
input_path = os.path.join(input_dir, input_file)
LOG_TRANSFORM = True  # Toggle to enable/disable log transformation
EPSILON = 1e-2        # Small constant to avoid log(0)

# -------------------------------
# Load and Inspect Data
# -------------------------------

if not os.path.exists(input_path):
    raise FileNotFoundError(f"Missing input file: {input_path}")

# Load CSV and convert time field
df = pd.read_csv(input_path)
df['Date.Time'] = pd.to_datetime(df['Date.Time'])

# Print basic metadata
print("Loaded file:", input_path)
print("Columns:", df.columns.tolist())
print("Total rows:", len(df))

# Filter to valid subdomains
valid_sdomains = ['SGS', 'SGN', 'SGI']
df_filtered = df[df['sdomain'].isin(valid_sdomains)].copy()
print(f"Filtered to subdomains {valid_sdomains}: {len(df_filtered)} rows")

# -------------------------------
# Monthly and Seasonal Sample Counts
# -------------------------------

# Add month and season columns
df_filtered['Month'] = df_filtered['Date.Time'].dt.month
season_map = {
    1: 'Winter', 2: 'Winter', 3: 'Winter',
    4: 'Spring', 5: 'Spring', 6: 'Spring',
    7: 'Summer', 8: 'Summer', 9: 'Summer',
    10: 'Fall', 11: 'Fall', 12: 'Fall'
}
df_filtered['Season'] = df_filtered['Month'].map(season_map)

# Helper: plot bar chart for counts by category
def plot_sample_counts(df_grouped, category, title, filename):
    bar_width = 0.25
    categories = df_grouped[category].unique()
    sdomains = df_grouped['sdomain'].unique()
    x = np.arange(len(categories))

    fig, ax = plt.subplots(figsize=(7, 4))
    for i, sdomain in enumerate(sdomains):
        subset = df_grouped[df_grouped['sdomain'] == sdomain]
        ax.bar(x + i * bar_width, subset['Sample Count'], bar_width, label=sdomain)

    ax.set_title(title)
    ax.set_xlabel(category)
    ax.set_ylabel('Sample Count')
    ax.legend(title='Subdomain')
    ax.set_xticks(x + bar_width * (len(sdomains) - 1) / 2)
    ax.set_xticklabels(categories)
    plt.tight_layout()
    plt.show()
    fig.savefig(filename)

# Monthly sample counts
monthly_counts = df_filtered.groupby(['Month', 'sdomain']).size().reset_index(name='Sample Count')
monthly_index = pd.MultiIndex.from_product([range(1, 13), valid_sdomains], names=['Month', 'sdomain'])
monthly_counts = monthly_counts.set_index(['Month', 'sdomain']).reindex(monthly_index, fill_value=0).reset_index()
plot_sample_counts(
    monthly_counts,
    category='Month',
    title='Number of Samples by Month and Subdomain',
    filename='..//..//figs//DIATOM_B_AVG_MONTH_2015-2018.png'
)

# Seasonal sample counts
seasonal_counts = df_filtered.groupby(['Season', 'sdomain']).size().reset_index(name='Sample Count')
seasonal_index = pd.MultiIndex.from_product([
    ['Winter', 'Spring', 'Summer', 'Fall'], valid_sdomains
], names=['Season', 'sdomain'])
seasonal_counts = seasonal_counts.set_index(['Season', 'sdomain']).reindex(seasonal_index, fill_value=0).reset_index()
plot_sample_counts(
    seasonal_counts,
    category='Season',
    title='Number of Samples by Season and Subdomain',
    filename='..//..//figs//DIATOM_B_AVG_SEASONAL_2015-2018.png'
)

# Station sample counts
station_season_counts = (
    df_filtered.groupby(['Station', 'Season'])
    .size()
    .reset_index(name='Sample Count')
)
fig, ax = plt.subplots(figsize=(10, 5))
station_order = station_season_counts['Station'].unique()
x = np.arange(len(station_order))
bar_width = 0.2

for i, season in enumerate(['Winter', 'Spring', 'Summer', 'Fall']):
    data = station_season_counts[station_season_counts['Season'] == season]
    counts = [data[data['Station'] == st]['Sample Count'].values[0] if st in data['Station'].values else 0 for st in station_order]
    ax.bar(x + i * bar_width, counts, bar_width, label=season)

ax.set_title('Sample Counts by Station and Season')
ax.set_xlabel('Station')
ax.set_ylabel('Sample Count')
ax.set_xticks(x + 1.5 * bar_width)
ax.set_xticklabels(station_order, rotation=90)
ax.legend(title='Season')
plt.tight_layout()
plt.show()


# -------------------------------
# Functional Group Construction
# -------------------------------
# List of raw count fields used in Nemcek-NAN and fields_to_average
raw_fields = ['Prasinophytes', 'Cryptophytes', 'Haptophytes', 'Dictyochophytes', 'Raphidophytes',
              'total diatoms', 'Cyanobacteria', 'PP1-DIA', 'PP2-NAN', 'PP3-PIC',
              'PZ1-CIL', 'ssc-DIA', 'ssc-FLA', 'ssc-CIL']

if LOG_TRANSFORM:
    for field in raw_fields:
        if field in df_filtered.columns:
            df_filtered[field] = np.log(df_filtered[field] + EPSILON)

# Define derived nanophytoplankton group from summed taxa
df_filtered['Nemcek-NAN'] = df_filtered[['Prasinophytes', 'Cryptophytes', 'Haptophytes',
                                         'Dictyochophytes', 'Raphidophytes']].sum(axis=1)

# Define fields to include in seasonal analysis
fields_to_average = [
    'total diatoms', 'Nemcek-NAN', 'Cyanobacteria', 'PP1-DIA', 'PP2-NAN',
    'PP3-PIC', 'PZ1-CIL', 'ssc-DIA', 'ssc-FLA', 'ssc-CIL'
]

print("Preparing seasonal means and anomalies for:")
print(fields_to_average)


# -------------------------------
# Seasonal Mean and Anomaly Calculation
# -------------------------------

# Step 1: De-duplicate by averaging values from same time/location
group_cols = ['closest_ecospace_time', 'ecospace_closest_lon', 'ecospace_closest_lat']
numeric_columns = df_filtered.select_dtypes(include='number').columns
non_numeric_columns = df_filtered.select_dtypes(exclude='number').columns

# Ensure grouping columns are excluded from numeric and non-numeric columns to avoid duplication errors
numeric_columns = [col for col in numeric_columns if col not in group_cols]
non_numeric_columns = [col for col in non_numeric_columns if col not in group_cols]

# Group by time and spatial cell, compute mean of numeric columns
df_numeric = df_filtered.groupby(group_cols)[numeric_columns].mean().reset_index()

# For non-numeric fields, take most common value (mode)
df_non_numeric = df_filtered.groupby(group_cols)[non_numeric_columns].agg(
    lambda x: x.mode().iloc[0] if not x.mode().empty else None
).reset_index()

# Merge numeric and non-numeric portions on time/space keys
df_simplified = pd.merge(df_numeric, df_non_numeric, on=group_cols, how='inner')
print("Length after deduplication:", len(df_simplified))

# Step 2: Seasonal averages by subdomain
seasonal_means = df_simplified.groupby(['Season', 'sdomain'])[fields_to_average].mean().reset_index()

# Step 3: Annual means and standard deviations by subdomain
annual_stats = df_simplified.groupby('sdomain')[fields_to_average].agg(['mean', 'std']).reset_index()

# Flatten multi-index columns from aggregation
annual_stats.columns = ['sdomain'] + [f"{var}_{stat}" for var, stat in annual_stats.columns[1:]]

# Step 4: Reshape and merge seasonal and annual stats for anomaly calc
seasonal_melted = seasonal_means.melt(id_vars=['Season', 'sdomain'], var_name='Field', value_name='Seasonal Mean')
annual_melted = annual_stats.melt(id_vars=['sdomain'], var_name='Field_Statistic', value_name='Value')
annual_melted[['Field', 'Statistic']] = annual_melted['Field_Statistic'].str.rsplit('_', n=1, expand=True)
annual_pivot = annual_melted.pivot_table(index=['sdomain', 'Field'], columns='Statistic', values='Value').reset_index()

# Merge to calculate anomaly
merged_data = pd.merge(seasonal_melted, annual_pivot, on=['sdomain', 'Field'])
merged_data['Anomaly'] = (merged_data['Seasonal Mean'] - merged_data['mean']) / merged_data['std']

print("Seasonal anomalies calculated. Example:")
print(merged_data[merged_data['Field'] == 'PP1-DIA'])

# -------------------------------
# Multi-Panel Seasonal Anomaly Plot
# -------------------------------

# Define grouped field sets for visualization
groups = {
    'Diatoms': ['total diatoms', 'PP1-DIA', 'ssc-DIA'],
    'Nanophytoplankton': ['Nemcek-NAN', 'PP2-NAN', 'ssc-FLA'],
    'Picophytoplankton': ['PP3-PIC', 'Cyanobacteria']
}

# Define color scheme
colors = {
    'Ecospace': 'blue',
    'SalishSeaCast': 'orange',
    'Observations': 'green'
}

def get_color(field):
    if field.startswith('PP') or field.startswith('PZ'):
        return colors['Ecospace']
    elif field.startswith('ssc'):
        return colors['SalishSeaCast']
    else:
        return colors['Observations']

# Subplot layout setup
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8, 8), sharex=False)
sdomains = ['SGN', 'SGS', 'SGI']
group_names = list(groups.keys())
ordered_seasons = ['Winter', 'Spring', 'Summer', 'Fall']

# Compute consistent y-limits for each subdomain (row)
y_limits = []
for sdomain in sdomains:
    all_anomalies = []
    for group in group_names:
        for field in groups[group]:
            subset = merged_data[(merged_data['sdomain'] == sdomain) & (merged_data['Field'] == field)]
            all_anomalies.extend(subset['Anomaly'].dropna().values)
    y_max = max(abs(np.min(all_anomalies)), abs(np.max(all_anomalies)))
    y_limits.append((-y_max * 1.1, y_max * 1.1))

# Plot loop
subplot_labels = ['(a)', '(d)', '(g)', '(b)', '(e)', '(h)', '(c)', '(f)', '(i)']
label_index = 0

for col, group in enumerate(group_names):
    for row, sdomain in enumerate(sdomains):
        ax = axes[row, col]
        bar_width = 0.2
        x = np.arange(len(ordered_seasons))

        for i, field in enumerate(groups[group]):
            subset = merged_data[(merged_data['sdomain'] == sdomain) & (merged_data['Field'] == field)]
            subset = subset.set_index('Season').loc[ordered_seasons].reset_index()
            ax.bar(x + i * bar_width, subset['Anomaly'], bar_width, label=field, color=get_color(field))

        ax.set_title(f'{sdomain} - {group}')
        if col == 0:
            ax.set_ylabel('Anomaly (z-score)', fontsize=9)
        ax.set_xticks(x + bar_width)
        ax.set_xticklabels(ordered_seasons, fontsize=8)
        ax.set_ylim(y_limits[row])
        ax.text(0.05, 0.95, subplot_labels[label_index], transform=ax.transAxes, fontsize=9,
                verticalalignment='top', horizontalalignment='left')
        label_index += 1

# Custom legend
custom_lines = [
    plt.Line2D([0], [0], color=colors['Ecospace'], lw=4),
    plt.Line2D([0], [0], color=colors['SalishSeaCast'], lw=4),
    plt.Line2D([0], [0], color=colors['Observations'], lw=4)
]
fig.legend(custom_lines, ['Ecospace', 'SalishSeaCast', 'Observations'],
           loc='center', ncol=3, bbox_to_anchor=(0.5, 0.02), fontsize=11)

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.show()
fig.savefig('..//..//figs//Nemcek_SSC_Ecosp_seasonal_anomaly_plots.png', bbox_inches='tight')


# -------------------------------
# Model-obs stats - exploration
# -------------------------------

# Define functional group pairings: (obs_field, ecospace_field)
comparison_pairs = {
    'Diatoms': ('total diatoms', 'PP1-DIA'),
    'Nanophytoplankton': ('Nemcek-NAN', 'PP2-NAN'),
    'Picophytoplankton': ('Cyanobacteria', 'PP3-PIC')
}

# Count valid comparisons per group, season, and subdomain
pair_counts = []

for group, (obs_var, model_var) in comparison_pairs.items():
    # Keep only rows where both obs and model are not NaN
    valid_rows = df_simplified.dropna(subset=[obs_var, model_var])

    # Count by Season and sdomain
    count_df = (
        valid_rows
        .groupby(['Season', 'sdomain'])
        .size()
        .reset_index(name='Count')
    )
    count_df['Group'] = group
    pair_counts.append(count_df)

# Combine into one DataFrame
pair_counts_df = pd.concat(pair_counts).sort_values(['Group', 'Season', 'sdomain'])

# Show result
print(pair_counts_df)
# -------------------------------
# Model-Observation Performance Metrics (Normalized Anomalies)
# -------------------------------

from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr
import numpy as np

# Define model-observation pairs for evaluation
anomaly_pairs = {
    'Diatoms': ('total diatoms', 'PP1-DIA'),
    'Nanophytoplankton': ('Nemcek-NAN', 'PP2-NAN'),
    'Picophytoplankton': ('Cyanobacteria', 'PP3-PIC')
}

# Prepare normalized anomalies in df_simplified (z-scores by subdomain)
normalized_df = df_simplified.copy()
normalized_df = normalized_df[normalized_df['closest_ecospace_time'] >= '1980-01-01']
for _, (obs_col, model_col) in anomaly_pairs.items():
    for var in [obs_col, model_col]:
        for domain in normalized_df['sdomain'].unique():
            domain_mask = normalized_df['sdomain'] == domain
            domain_vals = normalized_df.loc[domain_mask, var]
            mean_val = domain_vals.mean()
            std_val = domain_vals.std()
            norm_col = f"{var}_anom"
            if norm_col not in normalized_df.columns:
                normalized_df[norm_col] = np.nan
            normalized_df.loc[domain_mask, norm_col] = (domain_vals - mean_val) / std_val

# -------------------------------
# Compute metrics on anomalies
# -------------------------------

anomaly_results = []

for group_name, (obs_col, model_col) in anomaly_pairs.items():
    obs_anom = f"{obs_col}_anom"
    model_anom = f"{model_col}_anom"
    subset = normalized_df[[obs_anom, model_anom]].dropna()
    obs = subset[obs_anom].values
    model = subset[model_anom].values

    if len(obs) >= 5:
        bias = np.mean(model - obs)
        rmse = np.sqrt(mean_squared_error(obs, model))
        mae = mean_absolute_error(obs, model)
        r, _ = pearsonr(obs, model)
        wss = 1 - np.sum((model - obs) ** 2) / np.sum((np.abs(model - np.mean(obs)) + np.abs(obs - np.mean(obs))) ** 2)
    else:
        bias = rmse = mae = r = wss = np.nan

    anomaly_results.append({
        'Group': group_name,
        'Level': 'All',
        'Stratum': 'All',
        'N': len(obs),
        'Bias': bias,
        'RMSE': rmse,
        'MAE': mae,
        'R': r,
        'WSS': wss
    })

# Stratified by Season, Subdomain, and Season+Subdomain
strata = [('Season', s) for s in normalized_df['Season'].unique()] + \
         [('sdomain', d) for d in normalized_df['sdomain'].unique()] + \
         [('sdomain_Season', f"{row['sdomain']}_{row['Season']}")
          for _, row in normalized_df[['sdomain', 'Season']].drop_duplicates().iterrows()]

# Add combined stratum to DataFrame
normalized_df['sdomain_Season'] = normalized_df['sdomain'] + '_' + normalized_df['Season']

for level, stratum in strata:
    for group_name, (obs_col, model_col) in anomaly_pairs.items():
        obs_anom = f"{obs_col}_anom"
        model_anom = f"{model_col}_anom"
        subset = normalized_df[normalized_df[level] == stratum][[obs_anom, model_anom]].dropna()
        obs = subset[obs_anom].values
        model = subset[model_anom].values

        if len(obs) >= 5:
            bias = np.mean(model - obs)
            rmse = np.sqrt(mean_squared_error(obs, model))
            mae = mean_absolute_error(obs, model)
            r, _ = pearsonr(obs, model)
            wss = 1 - np.sum((model - obs) ** 2) / np.sum((np.abs(model - np.mean(obs)) + np.abs(obs - np.mean(obs))) ** 2)
        else:
            bias = rmse = mae = r = wss = np.nan

        anomaly_results.append({
            'Group': group_name,
            'Level': level,
            'Stratum': stratum,
            'N': len(obs),
            'Bias': bias,
            'RMSE': rmse,
            'MAE': mae,
            'R': r,
            'WSS': wss
        })

anomaly_df = pd.DataFrame(anomaly_results)
print("\nModel performance using normalized anomalies:")
print(anomaly_df.round(3))
print('done')
