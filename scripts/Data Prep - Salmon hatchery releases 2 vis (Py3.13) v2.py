# purpose: visuals for publication of hatchery releases over time.
# by: Greig Oldford
# created: Nov 2024
# purpose: previous ipynb didn't work anymore (issue with install of altair)
#           - this script produces visuals.


import pandas as pd
import matplotlib
matplotlib.use('TkAgg') # https://matplotlib.org/stable/users/explain/figure/backends.html
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# Load the data
localpath_epad = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/1. Salmon/All Species Hatchery Releases/EPADHatcherReleasesGST"
releases = pd.read_csv(localpath_epad + "/MODIFIED/actual_releases_COORDS.csv", parse_dates=['RELEASE_YEAR'])

# Map area codes to names
area_mappings = {
    'TOMM': 'Thompson & Shuswap R.',
    'TOMF': 'Thompson & Shuswap R.',
    'GSMN': 'Georgia Strait Mainland North',
    'GSMS': 'Georgia Strait Mainland South',
    'UPFR': 'Upper Fraser R.',
    'LWFR': 'Lower Fraser R.',
    'GSVI': 'Georgia Strait Vancouver Is.'
}
releases['Area'] = releases['STOCK_PROD_AREA_CODE'].map(area_mappings)

# Filter by species and aggregate by year
species_groups = ['Coho', 'Chinook', 'Chum', 'Pink', 'Sockeye']
species_data = {}
for species in species_groups:
    species_releases = releases[releases['SPECIES_NAME'] == species]
    species_agg = species_releases.groupby('RELEASE_YEAR')['TotalRelease'].sum()
    species_data[species] = species_agg

# convert to df
stacked_df = pd.DataFrame(species_data).fillna(0)  # fill NaN with 0 for any missing years
stacked_df.index = pd.to_datetime(stacked_df.index, format='%Y').year

# Plot
fig, ax = plt.subplots(figsize=(10, 6))
stacked_df.plot(kind='bar', stacked=True, ax=ax)
ax.set_title("Hatchery Releases Over Time")
ax.set_xlabel("Release Year")
ax.set_ylabel("Total Release")
# Set y-axis to regular notation
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f'{int(x):,}'))
ax.legend(title="Species")
ax.set_xticks(range(0, len(stacked_df.index), 5))
ax.set_xticklabels(stacked_df.index[::5], rotation=45)
plt.tight_layout()
plt.show()


# plot for each species
species_groups = ['Coho', 'Chinook', 'Chum', 'Pink', 'Sockeye']
for species in species_groups:
    species_data = {}
    species_releases = releases[releases['SPECIES_NAME'] == species]
    species_agg = species_releases.groupby('RELEASE_YEAR')['TotalRelease'].sum()
    species_data[species] = species_agg

    # convert to df
    stacked_df = pd.DataFrame(species_data).fillna(0)  # fill NaN with 0 for any missing years
    stacked_df.index = pd.to_datetime(stacked_df.index, format='%Y').year

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    stacked_df.plot(kind='bar', stacked=True, ax=ax)
    ax.set_title("Hatchery Releases Over Time")
    ax.set_xlabel("Release Year")
    ax.set_ylabel("Total Release")
    ax.legend(title="Species")
    ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, _: f'{int(x):,}'))
    plt.xticks(rotation=45)
    ax.set_xticks(range(0, len(stacked_df.index), 5))
    ax.set_xticklabels(stacked_df.index[::5], rotation=45)

    plt.tight_layout()
    plt.show()
