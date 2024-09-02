# Created: July 11, 2024, last edited Sep 1, 2024
# Author: G Oldford
# Purpose: Compare Nemcek and Ecospace, SSCast outputs
#          extracted from closest model points
# Input files:
#   - Nemcek survey data matched to SSC and Ecospace outputs
# Output files, fig:

# Source of Data:
#  - SalishSeaCast phyto data downloaded from url below using another script
#    https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.html
# - Nemcek, N., Hennekes, M., Sastri, A., & Perry, R. I. (2023). Seasonal and spatial dynamics of the phytoplankton
#      community in the Salish Sea, 2015–2019. Progress in Oceanography, 217, 103108.
#      https://doi.org/10.1016/j.pocean.2023.103108
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

# Load data
nemcek_matched_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
model_run_code = 'SC51'
nemcek_matched_f = "Nemcek_matched_to_model_out_" + model_run_code + ".csv"
nemcek_fp = os.path.join(nemcek_matched_p, nemcek_matched_f)

if os.path.exists(nemcek_fp):
    df = pd.read_csv(nemcek_fp)
    metadata = {
        "Column Names": df.columns.tolist(),
        "Number of Rows": len(df)
    }
    print(metadata)
else:
    print("Error: File not found.")

# Convert Date.Time to datetime
df['Date.Time'] = pd.to_datetime(df['Date.Time'])

# Filter out rows where 'sdomain' is empty
df_filtered = df[df['sdomain'].isin(['SGS', 'SGN', 'SGI'])]

###################################################
############### DATA DESCRIPTION ##################
#################### MONTHLY ######################
# Extract month from Date.Time
df_filtered['Month'] = df_filtered['Date.Time'].dt.month

# Group by 'Month' and 'sdomain' and count the number of samples
samples_by_month_sdomain = df_filtered.groupby(['Month', 'sdomain']).size().reset_index(name='Sample Count')

# Ensure all combinations of months and sdomain are present
all_months = np.arange(1, 13)
all_sdomains = ['SGS', 'SGN', 'SGI']
idx = pd.MultiIndex.from_product([all_months, all_sdomains], names=['Month', 'sdomain'])
samples_by_month_sdomain = samples_by_month_sdomain.set_index(['Month', 'sdomain']).reindex(idx, fill_value=0).reset_index()

# Bar Plot
fig, ax = plt.subplots(figsize=(7, 4))

# Define width of a single bar and calculate offsets
bar_width = 0.25
months = samples_by_month_sdomain['Month'].unique()
sdomain_list = samples_by_month_sdomain['sdomain'].unique()
num_sdomains = len(sdomain_list)

# Calculate x positions for each bar group
x = np.arange(len(months))

for i, sdomain in enumerate(sdomain_list):
    subset = samples_by_month_sdomain[samples_by_month_sdomain['sdomain'] == sdomain]
    ax.bar(x + i * bar_width, subset['Sample Count'], bar_width, label=sdomain)

ax.set_title('Number of Samples by Month and Province')
ax.set_xlabel('Month')
ax.set_ylabel('Sample Count')
ax.legend(title='Province')
ax.set_xticks(x + bar_width * (num_sdomains - 1) / 2)
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

plt.show()
fig.savefig('..//figs//' + 'DIATOM_B_AVG_MONTH_' + str(2015) + "-" + str(2018) + ".png")

####################################################
################# SEASONAL HISTO ###################
# Define seasons
seasons = {
    1: 'Winter', 2: 'Winter', 3: 'Winter',
    4: 'Spring', 5: 'Spring', 6: 'Spring',
    7: 'Summer', 8: 'Summer', 9: 'Summer',
    10: 'Fall', 11: 'Fall', 12: 'Fall'
}

# Map months to seasons
df_filtered['Season'] = df_filtered['Month'].map(seasons)

# Group by 'Season' and 'sdomain' and count the number of samples
samples_by_season_sdomain = df_filtered.groupby(['Season', 'sdomain']).size().reset_index(name='Sample Count')

# Ensure all combinations of seasons and sdomain are present
all_seasons = ['Winter', 'Spring', 'Summer', 'Fall']
all_sdomains = ['SGS', 'SGN', 'SGI']
idx = pd.MultiIndex.from_product([all_seasons, all_sdomains], names=['Season', 'sdomain'])
samples_by_season_sdomain = samples_by_season_sdomain.set_index(['Season', 'sdomain']).reindex(idx, fill_value=0).reset_index()

# Bar Plot
fig, ax = plt.subplots(figsize=(6, 4))

# Define width of a single bar and calculate offsets
bar_width = 0.25
seasons = samples_by_season_sdomain['Season'].unique()
sdomain_list = samples_by_season_sdomain['sdomain'].unique()
num_sdomains = len(sdomain_list)

# Calculate x positions for each bar group
x = np.arange(len(seasons))

for i, sdomain in enumerate(sdomain_list):
    subset = samples_by_season_sdomain[samples_by_season_sdomain['sdomain'] == sdomain]
    ax.bar(x + i * bar_width, subset['Sample Count'], bar_width, label=sdomain)

ax.set_title('Number of Samples by Season and Province')
ax.set_xlabel('Season')
ax.set_ylabel('Sample Count')
ax.legend(title='Province')
ax.set_xticks(x + bar_width * (num_sdomains - 1) / 2)
ax.set_xticklabels(all_seasons)

plt.show()
fig.savefig('..//figs//' + 'DIATOM_B_AVG_SEASONAL_' + str(2015) + "-" + str(2018) + ".png")


################################################
########## Visual of means by group ############

# mean by season, sdomain, model/obs, across all seasons
df_filtered['Nemcek-NAN'] = df_filtered[['Prasinophytes', 'Cryptophytes', 'Haptophytes', 'Dictyochophytes', 'Raphidophytes']].sum(axis=1)

# List of fields to calculate the average
fields_to_average = [
    'total diatoms', 'Nemcek-NAN', 'Cyanobacteria', 'PP1-DIA', 'PP2-NAN',
    'PP3-PIC', 'PZ1-CIL', 'ssc-DIA', 'ssc-FLA', 'ssc-CIL'
]

# Group by 'Season' and 'sdomain' and calculate the mean for the specified fields
averages_by_season_sdomain = df_filtered.groupby(['Season', 'sdomain'])[fields_to_average].mean().reset_index()

# # Define the groups
# groups = {
#     'Diatoms': ['total diatoms', 'PP1-DIA', 'ssc-DIA'],
#     'Ciliates': ['PZ1-CIL', 'ssc-CIL'],
#     'Nanoplankton': ['Nemcek-NAN', 'PP2-NAN', 'ssc-FLA'],
#     'Picoplankton': ['PP3-PIC', 'Cyanobacteria']
# }
#
# # Plot each group for each sdomain
# for sdomain in ['SGS', 'SGN', 'SGI']:
#     for group_name, fields in groups.items():
#         fig, ax = plt.subplots(figsize=(12, 6))
#
#         # Define width of a single bar and calculate offsets
#         bar_width = 0.25
#         seasons = averages_by_season_sdomain['Season'].unique()
#         num_fields = len(fields)
#
#         # Calculate x positions for each bar group
#         x = np.arange(len(seasons))
#
#         for i, field in enumerate(fields):
#             subset = averages_by_season_sdomain[(averages_by_season_sdomain['sdomain'] == sdomain)]
#             ax.bar(x + i * bar_width, subset[field], bar_width, label=field)
#
#         ax.set_title(f'Average {group_name} by Season for {sdomain}')
#         ax.set_xlabel('Season')
#         ax.set_ylabel('Mean Value')
#         ax.legend(title='Field')
#         ax.set_xticks(x + bar_width * (num_fields - 1) / 2)
#         ax.set_xticklabels(seasons)
#
#         plt.show()

# Calculate annual mean and standard deviation for each field within each sdomain
annual_stats = df_filtered.groupby('sdomain')[fields_to_average].agg(['mean', 'std']).reset_index()

# Merge annual stats with seasonal means
seasonal_merged = averages_by_season_sdomain.melt(id_vars=['Season', 'sdomain'], var_name='Field', value_name='Seasonal Mean')
annual_merged = annual_stats.melt(id_vars=['sdomain'], var_name=['Field', 'Statistic'], value_name='Value')
annual_merged = annual_merged.pivot_table(index=['sdomain', 'Field'], columns='Statistic', values='Value').reset_index()
merged_data = pd.merge(seasonal_merged, annual_merged, on=['sdomain', 'Field'])

# Calculate the seasonal anomaly (in terms of standard deviations)
merged_data['Anomaly'] = (merged_data['Seasonal Mean'] - merged_data['mean']) / merged_data['std']

print(merged_data[merged_data['Field']=='PP1-DIA'])

# Define the groups
groups = {
    'Diatoms': ['total diatoms', 'PP1-DIA', 'ssc-DIA'],
    #'Ciliates': ['PZ1-CIL', 'ssc-CIL'],
    'Nanoplankton': ['Nemcek-NAN', 'PP2-NAN', 'ssc-FLA'],
    'Picoplankton': ['PP3-PIC', 'Cyanobacteria']
}

# # Plot each group for each sdomain
# for sdomain in ['SGS', 'SGN', 'SGI']:
#     for group_name, fields in groups.items():
#         fig, ax = plt.subplots(figsize=(4, 3))
#
#         # Define width of a single bar and calculate offsets
#         bar_width = 0.25
#         seasons = merged_data['Season'].unique()
#         num_fields = len(fields)
#
#         # Calculate x positions for each bar group
#         x = np.arange(len(seasons))
#
#         for i, field in enumerate(fields):
#             subset = merged_data[(merged_data['sdomain'] == sdomain) & (merged_data['Field'] == field)]
#             ax.bar(x + i * bar_width, subset['Anomaly'], bar_width, label=field)
#
#         ax.set_title(f'Seasonal Anomaly ({group_name}) by Season for {sdomain}')
#         ax.set_xlabel('Season')
#         ax.set_ylabel('Anomaly (in standard deviations)')
#         ax.legend(title='Field')
#         ax.set_xticks(x + bar_width * (num_fields - 1) / 2)
#         ax.set_xticklabels(seasons)
#
#         plt.show()

# Define colors for each category
colors = {
    'Ecospace': 'blue',
    'SalishSeaCast': 'orange',
    'Observations': 'green'
}

# Function to get the color based on the field name
def get_color(field):
    if field.startswith('PP') or field.startswith('PZ'):
        return colors['Ecospace']
    elif field.startswith('ssc'):
        return colors['SalishSeaCast']
    else:
        return colors['Observations']

# Create the figure and axes for the subplots
fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(8, 8), sharex=False)

# Plot each group for each sdomain
sdomains = ['SGN', 'SGS', 'SGI']
group_names = list(groups.keys())
ordered_seasons = ['Winter', 'Spring', 'Summer', 'Fall']

# Calculate y-limits for each row
y_limits = []
for row, group_name in enumerate(group_names):
    all_anomalies = []
    for col, sdomain in enumerate(sdomains):
        for field in groups[group_name]:
            subset = merged_data[(merged_data['sdomain'] == sdomain) & (merged_data['Field'] == field)]
            all_anomalies.extend(subset['Anomaly'].values)
    y_lim_max_abs = max(abs(min(all_anomalies)),abs(max(all_anomalies))) # most extreme val, absolute
    y_limits.append(((y_lim_max_abs*-1) * 1.1, y_lim_max_abs * 1.1))

subplot_labels = [f'({chr(97 + i)})' for i in range(9)]
label_index = 0
for row, group_name in enumerate(group_names):
    for col, sdomain in enumerate(sdomains):
        ax = axes[row, col]

        # Define width of a single bar and calculate offsets
        bar_width = 0.2
        seasons = merged_data['Season'].unique()
        num_fields = len(groups[group_name])

        # Calculate x positions for each bar group
        x = np.arange(len(seasons))

        for i, field in enumerate(groups[group_name]):
            subset = merged_data[(merged_data['sdomain'] == sdomain) & (merged_data['Field'] == field)]
            subset = subset.set_index('Season').loc[ordered_seasons].reset_index()  # Ensure the order of seasons

            ax.bar(x + i * bar_width, subset['Anomaly'], bar_width, label=field, color=get_color(field))

        # Add text box label
        # ax.text(0.95, 0.95, f'{sdomain}\n{group_name}', transform=ax.transAxes, fontsize=12,
        #         verticalalignment='top', horizontalalignment='right',
        #         bbox=dict(boxstyle='round,pad=0.5', edgecolor='black', facecolor='white'))

        ax.set_title(f'{sdomain} - {group_name}')
        ax.set_xlabel('')
        if col == 0:
            ax.set_ylabel('Anomaly ($\sigma$)')
            ax.tick_params(axis='y', labelsize=9)

        ax.set_xticks(np.arange(len(ordered_seasons)) + bar_width * (num_fields - 1) / 2)
        ax.set_xticklabels(ordered_seasons)

        # ax.set_xticks(x + bar_width * (num_fields - 1) / 2)
        # ax.set_xticklabels(seasons)

        ax.tick_params(axis='x', labelsize=8)

        # Add subplot label
        ax.text(0.05, 0.95, subplot_labels[label_index], transform=ax.transAxes, fontsize=9,
                verticalalignment='top', horizontalalignment='left')
        label_index += 1

        # Set y-axis limits based on calculated limits for each row
        ax.set_ylim(y_limits[row])

# Create a custom legend
custom_lines = [
    plt.Line2D([0], [0], color=colors['Ecospace'], lw=4),
    plt.Line2D([0], [0], color=colors['SalishSeaCast'], lw=4),
    plt.Line2D([0], [0], color=colors['Observations'], lw=4)
]

fig.legend(custom_lines, ['Ecospace', 'SalishSeaCast', 'Observations'], loc='center', ncol=3, bbox_to_anchor=(0.5, 0.02), fontsize=12)

# Adjust layout and show the plot
plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.show()

# Save the figure as a PNG file
fig.savefig('../figs/seasonal_anomaly_plots.png', bbox_inches='tight')