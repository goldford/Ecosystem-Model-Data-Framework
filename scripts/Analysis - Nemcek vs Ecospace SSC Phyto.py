# Created: July 11, 2024
# Author: G Oldford
# Purpose: Compare Nemcek and Ecospace, SSCast outputs
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
nemcek_matched_f = "Nemcek_matched_to_model_out.csv"
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
#################### SEASONAL ######################
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
