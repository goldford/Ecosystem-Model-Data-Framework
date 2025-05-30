# Created by G Oldford
# July 29, 2024
# Purpose: compare bloom timing from ecospace to observations
#
# Source:
#
# Input:
#   1/ NC files of Ecospace out processed from ASC (other script)
#   2/ YAML outlining 'region 2' corresponding to Suchy et al 2021 region (satellite)
#   3/ bloom timing dates from central SoG from Suchy and Gower etc
#
# Output:
#   1/ figs of estimated bloom timing from Ecospace compared to observations
#   2/
# note:
import os
import xarray as xr
from helpers import read_sdomains, make_map, adjust_map, buildSortableString, find_bloom_doy
from matplotlib.path import Path as Path
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from datetime import datetime, timedelta

import netCDF4 as nc
import numpy as np
import pandas as pd


path_ecospace_out = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
# scenario = 'SC4_2'
# ecospace_code = "Scv4_2-MixingxPARLimitZ"
# scenario = 'SC4'
# ecospace_code = "Scv4-MultxPARlimZ"
# scenario = 'SC24'
# ecospace_code = "Scv24-PAR_PI_Temp" # this one was the 'hard reset' given new analysis method instead of using regional avgs
# scenario = 'SC26' #
# ecospace_code = "Scv26-PAR_PI_Temp_mixing"
# scenario = 'SC27' # introducig hab capacity driven by MixingXPAR
# ecospace_code = "Scv27-PAR_PI_Temp_mixing_habcap"
# scenario = 'SC28' #
# ecospace_code = "Scv28-PAR_PI_Temp_mixing"
# scenario = 'SC29' # very shallow PIC curve (no effect!)
# ecospace_code = "Scv29-PAR_shallowPI_Temp_mixing"
# scenario = 'SC30' # very shallow PIC curve, but with PAR as Enviro (SC29 PAR drove PP)
# ecospace_code = "Scv30-PARenv_shallowPI_Temp_mixing"
# scenario = 'SC31' # very shallow PI curve, but with PAR as Enviro, PARxMixing as habitat cap (looks bad for pico, nano)
# ecospace_code = "Scv31-PARenv_shallowPI_Temp_MixingXPARHab"
# scenario = 'SC32' # Shallow PI, new mixing curve - good result
# ecospace_code = "Scv32-PARenv_shallowPI_MixingSteep"
# scenario = 'SC33' # Less Shallow PI, new mixing curve
# ecospace_code = "Scv33-PARenv_lessshallowPI_MixingSteep"
# KEEP IN MIND MANY SCENS ABOVE ARE INFLUENCED BY INADVERTENT PP DRIVER using MIXINGXPAR
# - NO RESPONSE HOOKED UP BUT NONETHELESS THE EXTERNAL DRIVER WAS ACTIVE
# scenario = 'SC34' # Less Shallow PI, new mixing curve, w/ temp now
# ecospace_code = "Scv34-PARenv_lessshallowPI_MixingSteep"
# scenario = 'SC35' # Less Shallow PI, new mixing curve, w/ temp now -- this one can't trust - some enviro stuck off
# ecospace_code = "Scv35-PARenv_lessshallowPI_MixingSteep"
# scenario = 'SC36' # Less Shallow PI, new mixing curve, w/ temp now - - this one can't trust - some enviro stuck off
# ecospace_code = "Scv36-PARenv_lessshallowPI_MixingSteep"
# scenario = 'SC37' # Less Shallow PI, new mixing curve, w/ temp now - this one can't trust - some enviro stuck off
# ecospace_code = "Scv37-PARenv_lessshallowPI_MixingSteep_Temp"
# scenario = 'SC38' # Has temp, wind, PAR
# ecospace_code = "Scv38-PARenv_PI_Temp_Wind"
# scenario = 'SC39' # like 38 but adjusted wind repsonse
# ecospace_code = "Scv39-PARenv_PI_Temp_Wind"
# scenario = 'SC40' # like 39 but now w/ mixing to help with habitat cap
# ecospace_code = "Scv40-PARenv_PI_Temp_Wind"
# scenario = 'SC41' # like 39 or 40 (not sure which) but different results?! bug? Seems changes made in 40 did not take effect? RE runnin 40 gave similar result to 41
# ecospace_code = "Scv41-PARenv_PI_Temp_Wind"
# scenario = 'SC42' # is 41 but with the same PI as 39, to double check
# ecospace_code = "Scv42-PARenv_PI_Temp_Wind_Mixing"
scenario = 'SC43' # is 41 but with the same PI as 39, to double check
ecospace_code = "Scv43-All_Groups_Temp"

log_transform = True # following suchy, make true
# log transforming follows Suchy (though theirs is chl)
mean_or_median = "median" # use mean or median for def of bloom?
average_from = "annual"   # 'annual' to follow suchy, set bloom def based on annual avg (alt: "all")
thrshld_fctr = 1.05 # threshold above avg, suchy used 1.05
sub_thrshold_fctr = 0.7 # the secondary threshold, if first met, the test of 1/2 following 2 weeks (does it stay in bloom)
min_y_tick = 0 # 38 in previous plots
display_gower = False # true in most previous plots

v_f = {#"NK1-COH": path_ecospace_out + "EcospaceMapBiomass-NK1-COH-{}.asc",
       # "NK2-CHI": path_ecospace_out + "EcospaceMapBiomass-NK2-CHI-{}.asc",
       # "NK3-FOR": path_ecospace_out + "EcospaceMapBiomass-NK3-FOR-{}.asc",
       # "ZF1-ICT": path_ecospace_out + "EcospaceMapBiomass-ZF1-ICT-{}.asc",
       # "ZC1-EUP": path_ecospace_out + "EcospaceMapBiomass-ZC1-EUP-{}.asc",
       # "ZC2-AMP": path_ecospace_out + "EcospaceMapBiomass-ZC2-AMP-{}.asc",
       # "ZC3-DEC": path_ecospace_out + "EcospaceMapBiomass-ZC3-DEC-{}.asc",
       # "ZC4-CLG": path_ecospace_out + "EcospaceMapBiomass-ZC4-CLG-{}.asc",
       # "ZC5-CSM": path_ecospace_out + "EcospaceMapBiomass-ZC5-CSM-{}.asc",
       # "ZS1-JEL": path_ecospace_out + "EcospaceMapBiomass-ZS1-JEL-{}.asc",
       # "ZS2-CTH": path_ecospace_out + "EcospaceMapBiomass-ZS2-CTH-{}.asc",
       # "ZS3-CHA": path_ecospace_out + "EcospaceMapBiomass-ZS3-CHA-{}.asc",
       # "ZS4-LAR": path_ecospace_out + "EcospaceMapBiomass-ZS4-LAR-{}.asc",
       # "PZ1-CIL": path_ecospace_out + "EcospaceMapBiomass-PZ1-CIL-{}.asc",
       # "PZ2-DIN": path_ecospace_out + "EcospaceMapBiomass-PZ2-DIN-{}.asc",
       # "PZ3-HNF": path_ecospace_out + "EcospaceMapBiomass-PZ3-HNF-{}.asc",
       "PP1-DIA": path_ecospace_out + "EcospaceMapBiomass-PP1-DIA-{}.asc"#,
       # "PP2-NAN": path_ecospace_out + "EcospaceMapBiomass-PP2-NAN-{}.asc",
       # "PP3-PIC": path_ecospace_out + "EcospaceMapBiomass-PP3-PIC-{}.asc", #time step format eg: 00620,
      }

# read the data into an xarray object which is versatile
yr_strt = 2003
yr_end = 2018
mo_strt = 1
da_strt = 2
mo_end = 12
da_end = 30 #unused?

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

######### ECOSPACE OUT NC FILE ############
ecospace_file = ecospace_code + "_" + str(yr_strt) + "-" + str(yr_end) + ".nc"
ecospace_nc = os.path.join(path_ecospace_out, ecospace_file)
ds = xr.open_dataset(ecospace_nc)
print(ds)

############ REGION DEF FILE ###############
domain_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
# domain_f = "analysis_domains_jarnikova.yml"
domain_f = "analysis_domains_suchy.yml"
domain_fullp = os.path.join(domain_p, domain_f)
sdomains = read_sdomains(domain_fullp)

create_suchy_mask = False
if create_suchy_mask:
    # Extract the latitude and longitude arrays
    lat = ds['lat'].values
    lon = ds['lon'].values
    dep = ds['depth'].values

    # Initialize an empty array for the mask
    combined_mask = np.zeros(lat.shape, dtype=bool)

    # Loop over the polygon domains
    for sd, coords in sdomains.items():
        print(sd)
        if sd == 'SGC2': # modified to include broader area, so there's more mixing variability
            polygon_path = Path(coords)
            points = np.vstack((lat.flatten(), lon.flatten())).T
            mask = polygon_path.contains_points(points).reshape(lat.shape)
            combined_mask |= mask

    # Add the depth condition to the mask
    depth_mask = dep > 0
    combined_mask &= depth_mask

    # Create a new xarray Dataset with the mask
    mask_ds = xr.Dataset(
        {
            'mask': (('row', 'col'), combined_mask)
        },
        coords={
            'lat': (('row', 'col'), lat),
            'lon': (('row', 'col'), lon)
        }
    )

    # Save the mask to a new NetCDF file
    mask_ds.to_netcdf('..//data//evaluation//suchy_ecospace_mask.nc')

mask_ds = xr.open_dataset('..//data//evaluation//suchy_ecospace_mask.nc')
mask = mask_ds['mask'].values
# masked_data = ds.where(mask)
# masked_pp1_dia = masked_data['PP1-DIA']
# print(np.nanmean(masked_pp1_dia))

create_maps = False

# Loop over the years and select data based on the 'time' coordinate
annual_means = []
annual_medians = []
ts_reg_means = []
ts_reg_medians = []
timestamps = []
v = list(v_f)[0] # hard coding to DIA
for year in range(yr_strt, yr_end + 1):
    print(year)
    year_data = ds.sel(time=str(year))
    masked_ds = year_data.where(mask) # apply 'regional' mask
    if log_transform:
        masked_ds[v] = np.log(masked_ds[v]+0.01)
    masked_var = masked_ds[v]
    mean_value = np.nanmean(masked_var)
    annual_means.append(mean_value)
    median_value = np.nanmedian(masked_var)
    annual_medians.append(median_value)
    # print(f"Year {year}: Mean PP1-DIA = {mean_value}")
    # print(f"Year {year}: Median PP1-DIA = {median_value}")

    # select time period within each year (one per timestep)
    # yr_plt_strt = 2005
    # yr_plt_end = 2005
    mo_strt = 1
    mo_end = 12
    da_srt = 1
    da_end = 30

    # Select time period to plot figs (one per timestep)
    start_date = pd.Timestamp(year=year, month=mo_strt, day=da_srt)
    end_date = pd.Timestamp(year=year, month=mo_end, day=da_end)

    # Get the timestep start, end from dates
    time_coords = masked_ds.time.to_index()
    # ts_strt = time_coords.get_loc(start_date, method='bfill') # deprecated
    # ts_end = time_coords.get_loc(end_date, method='ffill') + 1
    ts_strt = time_coords.get_indexer([start_date], method='bfill')[0]
    ts_end = time_coords.get_indexer([end_date], method='ffill')[0] + 1

    map_num = 0
    for ts in range(ts_strt, ts_end):
    # for ts in range(ts_strt - 1, ts_end - 1):

        v1 = masked_var # to further develop, to loop over v's
        d1 = v1[ts]
        ts_mean = np.nanmean(d1)
        ts_median = np.nanmedian(d1)
        ts_reg_means.append(ts_mean)
        ts_reg_medians.append(ts_median)

        timestamp = pd.Timestamp(d1.time.values)
        timestamps.append(timestamp)
        year1 = timestamp.year
        month1 = timestamp.month
        day1 = timestamp.day

        if create_maps and map_num == 0: # remove map_num if you want map per time step
            ################### FIGS PER TIME STEP #######################
            # get the ecospace output indices for start and end from desired date range to plot
            # test map of ECOSPACE domain by row /col
            # note the bottom left should be fixed (nans)
            fig, axs = plt.subplots(figsize=(5, 8))
            # custom colormap to emphasise lower values
            colors = [(1, 1, 1), (0.8, 0.8, 0.8), (0.6, 0.6, 0.6), (0.4, 0.4, 0.4), (0, 0, 1), (0, 1, 0), (1, 1, 0),
                      (1, 0, 0)]
            # sc1
            # nodes = [0.0, 0.1, 0.2, 0.35, 0.6, 0.8, 0.9, 1.0] # Adjust these values to emphasize certain ranges
            # sc2
            nodes = [0.0, 0.15, 0.25, 0.3, 0.5, 0.65, 0.75, 1.0]
            custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))

            # specific month coloring experiments
            nodes = [0.0, 0.15, 0.25, 0.3, 0.69, 0.77, 0.86, 1.0]
            custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))

            # PPT colors, custom
            colors = [(1, 1, 1), (0.8, 0.8, 0.8), (0.6, 0.6, 0.6), (0.4, 0.4, 0.4), (0, 0, 1), (0, 1, 0), (1, 1, 0)]
            nodes = [0.0, 0.15, 0.25, 0.3, 0.45, 0.85, 1.0]
            custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))

            # sc 1
            vmin = 0
            vmax = 3

            w = plt.pcolormesh(d1.data, cmap=custom_cmap, vmin=vmin, vmax=vmax)
            title = ecospace_code + "-" + str(year) + "-" + str(month1) + "-" + str(day1) + "-" + str(ts)
            title = str(year) + "-" + str(month1) + "-" + str(day1) + "-" + str(ts)
            axs.set_title('{v} - {d}'.format(v=v, d=title))
            plt.colorbar(w, ax=axs)
            axs.invert_yaxis()
            # plt.show()

            # Save the figure in different formats
            month01 = buildSortableString(month1, 2)
            day01 = buildSortableString(day1, 2)
            # out_fig_name = "..//figs//{}_{}-{}{}{}-{}.{}"
            out_fig_name = "..//figs//TEMP_{}_{}-{}{}{}-{}.{}"
            # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "PDF"), format='pdf')
            out_fig_name = out_fig_name.format(ecospace_code, v, year, month01, day01, ts, "PNG")
            fig.savefig(out_fig_name, format='png')
            # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "EPS"), format='eps')
            print("Saved fig " + out_fig_name)
            map_num += 1


#####################################################
############## VISUALISE AS TIME SERIES #############
# visualise as time series
# Plot the regional average values in a TS
# using horizontal lines to show annual avg and vertical lines to show satellite derived bloom timing
v = list(v_f)[0]
if log_transform:
    outfigname = '..//figs//' + 'Ecospace_out_suchy_region_LOG' + v + '_' + scenario + '.png'
else:
    outfigname = '..//figs//' + 'Ecospace_out_suchy_region_' + v + '_' + scenario + '.png'

fig, axes = plt.subplots(1, 1, figsize=(10, 3), sharex=False)
axes.plot(timestamps, ts_reg_means, label=v, linewidth=0.6)
axes.plot(timestamps, ts_reg_medians, label="", linewidth=0.4)

axes.set_ylabel(v)

# Add vertical dashed lines at each date in dates_suchy
for date in dates_suchy:
    axes.axvline(date, color='black', linestyle='--', linewidth=1)

for date in dates_suchy:
    axes.axvline(datetime(date.year, 1, 1), color='blue', linestyle='--', linewidth=0.6)

# note I am using the suchy threshold (median * 1.05)

i = 0
for year in range(yr_strt, yr_end + 1):
    start_of_year = pd.Timestamp(f'{year}-01-01')
    end_of_year = pd.Timestamp(f'{year}-12-31')
    axes.hlines(annual_means[i] * 1.05, start_of_year, end_of_year, colors='red', linestyle='-', linewidth=1)
    axes.hlines(annual_medians[i] * 1.05, start_of_year, end_of_year, colors='orange', linestyle='-', linewidth=1)
    i += 1

plt.savefig(outfigname, bbox_inches='tight', dpi=300)
axes.legend()
plt.show()


################################################
############### DOY PLOT #######################
print("DoY bloom comparison plot")
# Create a DataFrame to display the results


if mean_or_median == "median":
    var_avg = ts_reg_medians
else:
    var_avg = ts_reg_means

# create a simple df to pass to function
ecospace_df = pd.DataFrame({
    'Year': pd.to_datetime(timestamps).year,
    'Date': pd.to_datetime(timestamps),
    v: var_avg
})

bloom_dates, bloom_days_of_year, bloom_earlylate = find_bloom_doy(ecospace_df, v,
                                                                  thrshld_fctr=thrshld_fctr,
                                                                  sub_thrshld_fctr=sub_thrshold_fctr,
                                                                  average_from=average_from,
                                                                  mean_or_median=mean_or_median
                                                                  )
ecospace_bloom_timing_df = pd.DataFrame({
    'Year': range(yr_strt, yr_end + 1),
    'Bloom Date': bloom_dates,
    'Day of Year': bloom_days_of_year,
    "Bloom Early Late": bloom_earlylate
})

# join the model and obs bloom info
bloom_timing = ecospace_bloom_timing_df
bloom_timing_df = bloom_timing.merge(suchy_bloom_timing_df, left_on='Year', right_on='Year')
bloom_timing_df = bloom_timing_df.merge(gower_bloom_timing_df, left_on='Year', right_on='Year_Gower', how='left')
bloom_timing_df['bloom_match'] = bloom_timing_df.apply(lambda row: row['Bloom Early Late_x'] == row['Bloom Early Late_y'], axis=1)
print(bloom_timing_df.columns)
df = bloom_timing_df

# Plotting
fig, ax = plt.subplots(figsize=(10, 6))

# Plotting the Ecospace model data
ax.errorbar(df['Year'], df['Day of Year_x'], yerr=3, fmt='s', color='blue', label='Ecospace Model', capsize=5)
# Plotting the Suchy data
ax.errorbar(df['Year'], df['Day of Year_y'], yerr=8, fmt='s', color='red', label='Suchy Data', capsize=5)

if display_gower:
    ax.errorbar(df['Year'], df['Day of Year_Gower'], yerr=7, fmt='s', color='orange', label='Gower Data', capsize=5)

# Adding horizontal dashed lines
ax.axhline(y=68, color='black', linestyle='--')
ax.axhline(y=108, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([x_min, x_max], 68, 108, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([2002.5, 2016.5])
y_min, y_max = ax.get_ylim()
ax.set_ylim([min_y_tick, y_max])
ax.set_xlabel('Year')
ax.set_ylabel('Day of Year')
# ax.set_title('Diatom Bloom Timing Comparison')
ax.legend()

if log_transform:
    trnsfrm = "logDIA"
else:
    trnsfrm = ""

plt.title(scenario + " " + trnsfrm)
# plt.savefig('..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_altmethod_' + scenario + '_' + trnsfrm + '.png')
plt.savefig('..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_suchymethod_' + scenario + '_' + trnsfrm + '_rev2.png')
plt.show()


################################################
############## STOP LIGHT PLOT #################
# Create a new DataFrame for the table
print(df)
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

plt.savefig('..//figs//' + scenario + '_' + trnsfrm + '_bloom_timing_comparison.pdf', bbox_inches='tight')
plt.savefig('..//figs//' + scenario + '_' + trnsfrm + '_bloom_timing_comparison.png', bbox_inches='tight', dpi=300)

# Display the table
plt.show()