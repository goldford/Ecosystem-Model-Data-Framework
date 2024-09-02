# Created by G Oldford
# July 23, 2024
# Purpose: compare bloom timing from ecospace to observations
#
# Source:
#
# Input: regional average from Ecospace (mean)
#
# Output: plots
#
# note: look at version 2. Best not to use the ecospace-averaged regional data (use median)


import pandas as pd
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import matplotlib.table as mtable
import numpy as np

pd.set_option('display.max_columns', None)  # Show all columns
pd.set_option('display.width', None)        # Auto-detect the width of the display and adjust accordingly
pd.set_option('display.max_colwidth', None) # Show full content of each column

# Paths
# Define the file path
# SC4
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v4 - MixingxPARLimitZ - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC4 v2
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v4_2 - MixingxPARLimitZ - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC7
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v7 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC8
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v8 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC9
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v9 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC10 - DIA responses: MixingXPAR, Temp, Light (PAR), Mixing
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v10 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC11 - DIA responses: MixingXPAR, Temp
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v11 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC12 - DIA responses: mixing (limiting at high vals), MixingXPAR (linear)
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v12 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC13 - DIA responses: mixing (light limiting at high vals), MixingXPAR, light, temp
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v13 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC14 - DIA responses: mixing (light limiting at high vals) adjusted, removed MixingxPAR, dinos are out of control
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v14 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC15 - DIA responses: mixing (light limiting at high vals) adjusted, removed MixingxPAR, dinos are out of control
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v15 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC16 - DIA responses: mixing (light limiting at high vals) adjusted, removed MixingxPAR, dinos are out of control
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v16 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC17 - DIA responses:
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v17 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC18 - timestep 24 p yr, trying now log DIA for eval. based on SC16, steeper PI
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v18 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC19 - based on old PI curve, changed the PAR to a PP driver, instead of enviro
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v19 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC20 -
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v20 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC21 -
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v21 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC22 -
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v22 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC23 - went back to PARxMixing
# file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v23 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"
# SC24 - PARxMixing, tweaked mixing response (more sensitive)
file_path = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v24 - PARMixingNut90Temp - 3DAY 1yr10yr//Ecospace_Average_Region_2_Biomass.csv"




# scenario = 'SC4'
# scenario = 'SC4_2'
scenario = 'SC24'

# path_ecospace_out = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v4 - MixingxPARLimitZ - 3DAY 1yr10yr"
# path_ecospace_out = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v7 - PARMixingNut90Temp - 3DAY 1yr10yr"
# path_ecospace_file = "Ecospace_Average_Region_2_Biomass.csv"

doy_suchy = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77]
doy_gower = [81, 72, 58, 84, 48, 72, 58] # gower via allen

dates_suchy = []
for i, doy in enumerate(doy_suchy):
    year = 2003 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dates_suchy.append(date)

suchy_bloom_timing_df = pd.DataFrame({
    'Year': range(2003, 2017),
    'Bloom Date': dates_suchy,
    'Day of Year': doy_suchy,
    "Bloom Early Late": ['avg', 'avg', 'early', 'avg', 'late', 'late', 'avg', 'avg',
                         'avg', 'avg', 'avg', 'avg', 'early', 'avg']
})

dates_gower= []
for i, doy in enumerate(doy_gower):
    year = 2003 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dates_gower.append(date)

gower_bloom_timing_df = pd.DataFrame({
    'Year_Gower': range(2003, 2009 + 1),
    'Bloom Date_Gower': dates_gower,
    'Day of Year_Gower': doy_gower,
    "Bloom Early Late Gower": ['na', 'na', 'na', 'na', 'na', 'na', 'na']
})

# Load the data, specifying the header to be on the 33rd line (32 index because it's zero-based)
# also if opened in Excel and saved I think header lines change
if (scenario == 'SC4') or (scenario == 'SC4_2'):
    header = 30
elif ((scenario == 'SC19') or (scenario == 'SC20') or
      (scenario == 'SC21') or (scenario == 'SC22') or
      (scenario == 'SC23') or (scenario == 'SC24')):
    header = 30
elif scenario == 'SC16':
    header = 33
else:
    header = 31

timesteps_pr_year = 12

def reduce_by_half(df):
    # Reshape the DataFrame by stacking pairs of rows
    reshaped = df.values.reshape(-1, 2, df.shape[1])
    # Compute the mean along the new axis
    averaged = reshaped.mean(axis=1)
    reduced_df = pd.DataFrame(averaged, columns=df.columns)
    return reduced_df


ecospace_reg2_df = pd.read_csv(file_path, header=header)
print(ecospace_reg2_df.columns)
diatoms_ecospace_df = ecospace_reg2_df[['TimeStep','PP1-DIA']]
if timesteps_pr_year == 24:
    diatoms_ecospace_df = reduce_by_half(diatoms_ecospace_df)
diatoms_ecospace_df['logPP1-DIA'] = np.log(diatoms_ecospace_df['PP1-DIA'])
diatoms_ecospace_df['Date'] = datetime(2003, 1, 2)

# Define the start date
start_date = datetime(2003, 1, 2)
dates = []

# Iterate through each year until December 31, 2018
cntr = 1
year = start_date.year
for index, row in diatoms_ecospace_df.iterrows():

    time_step = row['TimeStep']

    if cntr == 1:
        set_date = datetime(year, 1, 2)
        cntr += 1
    else:
        if cntr <= 120:
            if cntr < 120:
                set_date += timedelta(days=3)
                cntr += 1
            else:  # Handle the 120th time step
                set_date = datetime(year, 12, 28)
                year = year + 1
                cntr = 1

    diatoms_ecospace_df.at[index, 'Date'] = set_date


diatoms_ecospace_df['Day of Year'] = diatoms_ecospace_df['Date'].apply(lambda x: x.strftime('%j')).astype(int)
print(diatoms_ecospace_df.head())
# Display the first few rows to verify

# put copy in helpers 2024-07-30
def find_diatom_blooms(df, df_field):
    # Create lists to hold the results
    bloom_dates = []
    bloom_days_of_year = []
    bloom_earlylate = []

    # match spreadsheet method:
    # method: median est set fixed to 1.05x 1.65 based on med of all yrs
    # they don't quite match (2008 and 2011 are the issues)
    # if median value is set to 1.65 and is based on median across all years and threshold is high
    # secondary criteria is average of one of two following weeks must be 0.95 threshold.
    # then all years except 3 match: 2006, 2008, 2011. 2006 is probably never going to work,
    # 2008 is on cusp, 2011 is issue because either early or late but should be avg
    # if you use a median or threshold val of aroun 1.65 and a higher sub-threshold than suchy then it works
    # results overall indicate the model is generally slightly early
    #
    # match suchy method (median from within year, threshold 70% median
    # results: model tool early consistently

    # Group by year
    df['Year'] = df['Date'].dt.year
    # median_value = df[df_field].median()
    # mean_value = df[df_field].mean()
    # median_value = 1.65
    # median_value = 1.9
    grouped = df.groupby('Year')

    # Iterate through each year
    for year, group in grouped:

        median_value = group[df_field].median()
        threshold = 1.05 * median_value
        # mean_value = group[df_field].mean()
        # threshold = 1.05 * mean_value

        print("year " + str(year))
        print("threshold " + str(threshold))
        sub_threshold = 0.7 * threshold
        # sub_threshold = 0.95 * threshold

        bloom_found = False

        # Iterate through the time steps of the year
        for i in range(len(group) - 2):  # Ensure we have at least two subsequent steps
            if group.iloc[i][df_field] >= threshold:
                # print("possible bloom found")
                # print(group.iloc[i + 1]['PP1-DIA'])
                # print(group.iloc[i + 2]['PP1-DIA'])
                # print(group.iloc[i + 3]['PP1-DIA'])
                if (
                        (group.iloc[i + 1][df_field] + group.iloc[i + 2][df_field]) / 2
                        > sub_threshold
                ) or (
                        (group.iloc[i + 3][df_field] + group.iloc[i + 4][df_field]) / 2
                        > sub_threshold
                ):
                    # print("bloom passes second crit")
                    bloom_dates.append(group.iloc[i]['Date'])
                    bloom_doy = group.iloc[i]['Date'].timetuple().tm_yday
                    bloom_days_of_year.append(bloom_doy)
                    bloom_found = True
                    print(group.iloc[i + 1][df_field])
                    if bloom_doy < 68:
                        bloom_earlylate.append("early")
                    elif bloom_doy <= 107:
                        bloom_earlylate.append("avg")
                    elif bloom_doy >= 108:
                        bloom_earlylate.append("late")
                    else:
                        print("problem finding bloom timing category")
                    break

        # If no bloom is found for the year, append None
        if not bloom_found:
            bloom_dates.append(None)
            bloom_days_of_year.append(None)
            bloom_earlylate.append(None)

    return bloom_dates, bloom_days_of_year, bloom_earlylate

# specify whether log transform
log_DIA = True

# Assuming diatoms_ecospace_df is already defined and populated
if log_DIA:
    df_field = 'logPP1-DIA'
else:
    df_field = 'PP1-DIA'
bloom_dates, bloom_days_of_year, bloom_earlylate = find_diatom_blooms(diatoms_ecospace_df, df_field)

print(bloom_dates)
print(bloom_days_of_year)

# Create a DataFrame to display the results
ecospace_bloom_timing_df = pd.DataFrame({
    'Year': range(2003, 2019),
    'Bloom Date': bloom_dates,
    'Day of Year': bloom_days_of_year,
    "Bloom Early Late": bloom_earlylate

})

bloom_timing = ecospace_bloom_timing_df
bloom_timing_df = bloom_timing.merge(suchy_bloom_timing_df, left_on='Year', right_on='Year')
bloom_timing_df = bloom_timing_df.merge(gower_bloom_timing_df, left_on='Year', right_on='Year_Gower', how='left')


bloom_timing_df['bloom_match'] = bloom_timing_df.apply(lambda row: row['Bloom Early Late_x'] == row['Bloom Early Late_y'], axis=1)
print(bloom_timing_df)
# If you need to save it to a file, you can use the following line
# df_dates.to_csv("TimeSteps_Dates.csv", index=False)

df = bloom_timing_df

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

# Plotting the Ecospace model data
ax.errorbar(df['Year'], df['Day of Year_x'], yerr=3, fmt='s', color='blue', label='Ecospace Model', capsize=5)
# Plotting the Suchy data
ax.errorbar(df['Year'], df['Day of Year_y'], yerr=7, fmt='s', color='red', label='Suchy Data', capsize=5)
ax.errorbar(df['Year'], df['Day of Year_Gower'], yerr=7, fmt='s', color='orange', label='Gower Data', capsize=5)

# Adding horizontal dashed lines
ax.axhline(y=68, color='black', linestyle='--')
ax.axhline(y=108, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([x_min, x_max], 68, 108, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([2002.5, 2016.5])
ax.set_xlabel('Year')
ax.set_ylabel('Day of Year')
# ax.set_title('Diatom Bloom Timing Comparison')
ax.legend()

if log_DIA:
    trnsfrm = "logDIA"
else:
    trnsfrm = ""

# plt.savefig('..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_altmethod_' + scenario + '_' + trnsfrm + '.png')
plt.savefig('..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_suchymethod_' + '_' + trnsfrm + '.png')
plt.show()

# Create a new DataFrame for the table
table_df = df[['Year', 'Bloom Early Late_x', 'Bloom Early Late_y']]
table_df.columns = ['Year', 'Ecospace', 'Suchy']

# Define a function to apply the row colors
def color_rows(row):
    if row['Ecospace'] == row['Suchy']:
        return ['background-color: lightgreen'] * len(row)
    else:
        return ['background-color: lightcoral'] * len(row)

# Apply the row colors
styled_table = table_df.style.apply(color_rows, axis=1)

# Display the styled table
styled_table.to_html('..//figs//' + 'bloom_timing_comparison_' + scenario + '.html')


###################################
fig, ax = plt.subplots(figsize=(3, 4))  # Adjust the figure size as needed

# Hide axes
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)

# Create the table
table = ax.table(
    cellText=table_df.values,
    colLabels=table_df.columns,
    cellLoc='center',
    loc='center',
    cellColours=[
        ['lightgreen' if row['Ecospace'] == row['Suchy'] else 'lightcoral']*len(row)
        for _, row in table_df.iterrows()
    ]
)

# Set font size
table.set_fontsize(9)
table.scale(1.2, 1.2)

plt.savefig('..//figs//' + 'bloom_timing_comparison.pdf', bbox_inches='tight')
plt.savefig('..//figs//' + 'bloom_timing_comparison.png', bbox_inches='tight', dpi=300)

# Display the table
plt.show()

#####################################
# Create the plot
fig, ax = plt.subplots(figsize=(3, 4))  # Adjust the figure size as needed

# Hide axes
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_frame_on(False)

# Create the table
table = mtable.Table(ax, bbox=[0, 0, 1, 1])

# Function to color the cells
def cell_colors(val1, val2):
    if val1 == val2:
        return 'lightgreen'
    else:
        return 'lightcoral'

# Add table cells
nrows, ncols = table_df.shape
for i in range(nrows + 1):
    for j in range(ncols):
        if i == 0:
            # Header row
            text = table_df.columns[j]
            weight = 'bold'
            fontsize = 14
            facecolor = 'white'
            edgecolor = 'black'
            linewidth = 2
        else:
            # Data rows
            text = table_df.iloc[i - 1, j]
            weight = 'normal'
            fontsize = 9
            facecolor = cell_colors(table_df.iloc[i - 1, 1], table_df.iloc[i - 1, 2])
            edgecolor = 'black'
            linewidth = 1

        table.add_cell(i, j, width=1 / ncols, height=1 / (nrows + 1), text=text,
                       loc='center', facecolor=facecolor, edgecolor=edgecolor)

# Add the table to the axes
ax.add_table(table)

plt.savefig('..//figs//' + 'bloom_timing_comparison_' + scenario + '.pdf', bbox_inches='tight')
plt.savefig('..//figs//' + 'bloom_timing_comparison_' + scenario + '.png', bbox_inches='tight', dpi=300)

# Display the table
plt.show()

#####################################################
############## VISUALISE AS TIME SERIES #############

# visualise as time series
if log_DIA:
    var = 'logPP1-DIA'
    outfigname = '..//figs//' + 'region2_LOGpp1-dia_' + scenario + '.png'
else:
    var = 'PP1-DIA'
    outfigname = '..//figs//' + 'region2_pp1-dia_' + scenario + '.png'

fig, axes = plt.subplots(1, 1, figsize=(10, 3), sharex=False)
axes.plot(diatoms_ecospace_df['Date'], diatoms_ecospace_df[var], label=var, linewidth=0.6)
axes.set_ylabel(var)
plt.savefig(outfigname, bbox_inches='tight', dpi=300)

axes.legend()
plt.show()

