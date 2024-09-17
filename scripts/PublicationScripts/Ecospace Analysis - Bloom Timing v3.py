# Created by G Oldford
# Aug 29, 2024
# Purpose: compare bloom timing from ecospace to observations
#
# Input:
#   1/ NC files of Ecospace out processed from ASC by another script
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
from helpers import find_nearest_point_from1D, find_nearest_point
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
# scenario = 'SC43' # is 41 but with the same PI as 39, to double check
# ecospace_code = "Scv43-All_Groups_Temp"
# scenario = 'SC45' # is 41 but with the same PI as 39, to double check
# ecospace_code = "Scv45-PARPI_Temp_Wind"
# SCENARIOS 45 vs 46 tell us that temp responses applied to other PP groups let DIA get more nutrients and grow faster
#                    thus we get earlier bloom.
# scenario = 'SC46' # is 41 but with the same PI as 39, to double check
# ecospace_code = "Scv46-PARPI_AllTemp_Wind"
# Scenario 47 is done to tune wind such that we get a little later bloom timing again, after adding the temp responses
#               for the other PP groups in 46 led to bloom timing that was generally a bit too early
# scenario = 'SC47' # is 46 but with slope changed of wind for DIA to try to compensate for effects seen in 47
# ecospace_code = "Scv47-PARPI_AllTemp_Wind"
# scenario = 'SC48' # is 47 but with ecopath biomass instead of habitat adjusted used
# ecospace_code = "Scv48-PARPI_AllTemp_Wind"
# scenario = 'SC49' # is 48 but with wind response back to what is used in 46 and prior (keeping the Ecopath base B at init)
# ecospace_code = "Scv49-PARPI_AllTemp_Wind"
# scenario = 'SC50' # is 49 but with 'habitat adjusted biomasses' and right-shoulder PI PAR env resp functions
# ecospace_code = "Scv50-RSPI_AllTemp_Wind"
# scenario = 'SC51' # is 50 but with forcings for 2000-2002
# ecospace_code = "Scv51-RSPI_AllTemp_Wind"
# scenario = 'SC52' # is 51 but with ecopath initial B's
# ecospace_code = "Scv52-RSPI_AllTemp_Wind"
# scenario = 'SC53' # is 52 but now with proxy shallow temp PROD response in lieu of nutrients, applied to DIA
# ecospace_code = "Scv53-RSPI_AllTemp_Wind"
# scenario = 'SC54' # is 53 but now with proxy shallow temp MORT response in lieu of nutrients, applied to DIA
# ecospace_code = "Scv54-RSPI_AllTemp_Wind"
# scenario = 'SC56' # is 54, 55 but with Hab Cap as dynamic driver directly by light, keeping also the responses w/
#                   # like proxy shallow temp MORT response in lieu of nutrients, applied to DIA
#                   # The hab cap driver seems to bump the variability of growth rates up
#                   # whereas the enviro response method seems only to penalise
# ecospace_code = "Scv56-RSPI_AllTemp_Wind"
# scenario = 'SC58' # is 54, 55 but with Hab Cap as dynamic driver directly by light, keeping also the responses w/
#                   # like proxy shallow temp MORT response in lieu of nutrients, applied to DIA
#                   # The hab cap driver seems to bump the variability of growth rates up
#                   # so removed temp penalty driver on DIA to see how much variability is reduced
# ecospace_code = "Scv58-RSPI_AllTemp_Wind"
# scenario = 'SC59' # is 58 but with no hab capacity as driver and with spinup data (day 200 repeated)
#                   # and B initialised is Ecospace's habitat adjusted B
# ecospace_code = "Scv59-RSPI_AllTemp_Wind"
# scenario = 'SC64' # is 58 but with no hab capacity as driver and with spinup data (day 200 repeated)
#                   # and B initialised is Ecospace's habitat adjusted B
# ecospace_code = "Scv64-RSPI_AllTemp_Wind"
# scenario = 'SC56_2' # 56, trying to recreate it, ecopath adj b's, pbmax dia 20
# ecospace_code = "Scv56_2-RSPI_AllTemp_Wind"
# scenario = 'SC56_3' # 56, trying to recreate it, ecopath adj b's, pbmax 1000
# ecospace_code = "Scv56_3-RSPI_AllTemp_Wind"
# scenario = 'SC56_4' # 56, trying to recreate it, ecopath hab adj b's, pbmax 1000
# ecospace_code = "Scv56_4-RSPI_AllTemp_Wind"
# scenario = 'SC56_5' # 56, trying to recreate it, ecopath hab adj b's, pbmax 20 all
# ecospace_code = "Scv56_5-RSPI_AllTemp_Wind"
# scenario = 'SC51_2' # 56, trying to recreate it, ecopath hab adj b's, pbmax 20 all
# ecospace_code = "Scv51_2-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC51_3' # 56, trying to recreate it, ecopath hab adj b's, pbmax 20 all
# ecospace_code = "Scv51_3-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC50_2' # 50 attempt 2, to recreate...  ecosim nut_f 90, PBMax 1000, ecopath hab adjst B
# ecospace_code = "Scv50_2-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC50_3' # 50 attempt to recreate...  ecosim nut_f 95, PBMax 1000, ecopath hab adjst B
# ecospace_code = "Scv50_3-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC51_4' #
# ecospace_code = "Scv51_4-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC51_4_2' #
# ecospace_code = "Scv51_4_2-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC70' # same as 4_2, debug mode, but with 4 threads
# ecospace_code = "Scv70-PAR_PI_AllPPTemp_Wind"
# scenario = 'SC71' # same as 70, with all v's = 2 and with 4 threads
# ecospace_code = "Sc71-PAR_PI_AllPPTemp_Wind"
scenario = 'FULLKEY_SC51_5'
ecospace_code = "FULLKEY_Scv51_5-PAR_PI_AllPPTemp_Wind"

filenm_yr_strt = 1978

log_transform = True # following suchy, make true
# log transforming follows Suchy (though theirs is chl)
mean_or_median = "median" # use mean or median for def of bloom?
average_from = "annual"   # 'annual' to follow suchy, set bloom def based on annual avg (alt: "all")
thrshld_fctr = 1.05 # threshold above avg, suchy used 1.05
sub_thrshold_fctr = 0.7 # the secondary threshold, if first met, the test of 1/2 following 2 weeks (does it stay in bloom)
min_y_tick = 30 # 38 in previous plots
display_gower = False # true in most previous plots
display_allen = True

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
yr_strt = 1980
yr_end = 2018
mo_strt = 1
da_strt = 2
mo_end = 12
da_end = 30 #unused?

# first indices below are for year 2003
doy_suchy = [100, 68, 50, 83, 115, 115, 100, 100, 92, 88, 92, 92, 55, 77]
doy_gower = [81, 72, 58, 84, 48, 72, 58] # gower via allen

# first indices below are 1980 (+/- 4 days should be used for error bars)
doy_allen = [94, 78, 81, 82, 87, 76, 75, 87, 99, 76,
             81, 78, 77, 55, 86, 86, 67, 87, 77, 104,
             66, 70, 92, 86, 81, 62, 88, 100 ,90, 97, 104] # from 1D hindcast, peak bloom dates

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

dates_allen = []
for i, doy in enumerate(doy_allen):
    year = 1980 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dates_allen.append(date)

allen_bloom_timing_df = pd.DataFrame({
    'Year_Allen': range(1980, 2010 + 1),
    'Bloom Date_Allen': dates_allen,
    'Day of Year_Allen': doy_allen,
    "Bloom Early Late Allen": np.full(len(range(1980,2010+1)), 'na')
})


######### ECOSPACE OUT NC FILE ############
ecospace_file = ecospace_code + "_" + str(filenm_yr_strt) + "-" + str(yr_end) + ".nc"
ecospace_nc = os.path.join(path_ecospace_out, ecospace_file)
ds = xr.open_dataset(ecospace_nc)
print(ds)

############################################
############ MASKS - REGIONS ###############
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
    mask_ds.to_netcdf('..//..//data//evaluation//suchy_ecospace_mask.nc')

mask_ds = xr.open_dataset('..//..//data//evaluation//suchy_ecospace_mask.nc')
mask = mask_ds['mask'].values
# masked_data = ds.where(mask)
# masked_pp1_dia = masked_data['PP1-DIA']
# print(np.nanmean(masked_pp1_dia))

#############################################
############## MASK ALLEN 1D ################

create_allen_mask = False
if create_allen_mask:

    lat_allen1D = 49.1887166
    lon_allen1D = -123.4966697
    lats_mod = ds['lat'].values
    lons_mod = ds['lon'].values
    deps_mod = ds['depth'].values
    mask_ecospace = deps_mod != 0
    dx = 0.01  # 0.01 deg lat = ~1.1 km
    allen1d_idx_and_depth = find_nearest_point(lon_allen1D, lat_allen1D, lons_mod, lats_mod, mask_ecospace, dx)
    nearest_idx = allen1d_idx_and_depth[0]
    if not nearest_idx is np.nan:
        print('found it')
        print(allen1d_idx_and_depth)
    else:
        print('did not find it')

allen1d_idx_and_depth = [52, 100, 370.68523301071156] # row, col, depth of allen 1d model location


#############################################
# get the data for a single model point
annual_means_allen = [] # allen refers to location of allen 1D model
annual_medians_allen = []
ts_reg_means_allen = []
ts_reg_medians_allen = []
timestamps_allen = []

v = list(v_f)[0] # hard coding to DIA
print('ecospace model data at allen 1D model location')
for year in range(yr_strt, yr_end + 1):
    ############################################
    # pull model results for location corresponding to allen's 1D model

    print(year)
    year_data = ds.sel(time=str(year))
    masked_ds = ds.sel(row=allen1d_idx_and_depth[1],
                       col=allen1d_idx_and_depth[0])
    if log_transform:
        masked_ds[v] = np.log(masked_ds[v]+0.01)
    masked_var = masked_ds[v]
    mean_value = np.nanmean(masked_var)
    annual_means_allen.append(mean_value)
    median_value = np.nanmedian(masked_var)
    annual_medians_allen.append(median_value)

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

    for ts in range(ts_strt, ts_end):
        # for ts in range(ts_strt - 1, ts_end - 1):

        v1 = masked_var  # to further develop, to loop over v's
        d1 = v1[ts]
        ts_mean = np.nanmean(d1)
        ts_median = np.nanmedian(d1)
        ts_reg_means_allen.append(ts_mean)
        ts_reg_medians_allen.append(ts_median)

        timestamp = pd.Timestamp(d1.time.values)
        timestamps_allen.append(timestamp)

print(len(ts_reg_medians_allen))



################################################
############### DOY PLOT 1 - long ###############
### Shows Allen's 1D versus Ecospace output ###
### Bloom timing based on Allen 1D location ###
print("DoY bloom comparison plot 1")
# Create a DataFrame to display the results

xlim_min = 1979.5
xlim_max = 2016.5

if mean_or_median == "median":
    var_avg = ts_reg_medians_allen
else:
    var_avg = ts_reg_means_allen

# create a simple df to pass to function
ecospace_df = pd.DataFrame({
    'Year': pd.to_datetime(timestamps_allen).year,
    'Date': pd.to_datetime(timestamps_allen),
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
bloom_timing_df = ecospace_bloom_timing_df.merge(allen_bloom_timing_df, left_on='Year', right_on='Year_Allen', how='left')
print(bloom_timing_df.columns)
bloom_timing_df['bloom_match'] = bloom_timing_df.apply(lambda row: row['Bloom Early Late'] == row['Bloom Early Late Allen'], axis=1)
print(bloom_timing_df.columns)
df = bloom_timing_df

# Plotting
fig_width = 12
fig_height = 6
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plotting the Ecospace model data
ax.errorbar(df['Year'], df['Day of Year'], yerr=3, fmt='s', color='blue', label='Ecospace Model', capsize=3)
ax.errorbar(df['Year'], df['Day of Year_Allen'], yerr=8, fmt='s', color='green', label='C09 1D Model', capsize=3)

# Adding horizontal dashed lines
ax.axhline(y=68, color='black', linestyle='--')
ax.axhline(y=108, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([xlim_min, x_max], 68, 108, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([xlim_min, xlim_max])
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
plt.savefig('..//..//figs//' + 'smBloomTiming_Ecospace_vs_Allen_' + str(yr_strt) + '_' + str(yr_end) + 'alttosuchymethod_' + scenario + '_' + trnsfrm + '_rev2.png')
plt.show()

exit()


#############################################
create_maps = False

# Loop over the years and select data based on the 'time' coordinate
annual_means = []
annual_medians = []
ts_reg_means = []
ts_reg_medians = []
timestamps = []

v = list(v_f)[0] # hard coding to DIA
for year in range(yr_strt, yr_end + 1):
    ############################################
    # pull model results for 'region' (mask corresponding to suchy)
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
            out_fig_name = "..//..//figs//TEMP_{}_{}-{}{}{}-{}.{}"
            # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "PDF"), format='pdf')
            out_fig_name = out_fig_name.format(ecospace_code, v, year, month01, day01, ts, "PNG")
            fig.savefig(out_fig_name, format='png')
            # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "EPS"), format='eps')
            print("Saved fig " + out_fig_name)
            map_num += 1



################################################
############### DOY PLOT 2 #####################
### Shows Suchy et al and optionally Allen's ###
### With bloom timing based on suchy region ####

print("DoY bloom comparison plot")
# Create a DataFrame to display the results

xlim_min = 1979.5
xlim_max = 2016.5

if display_allen: # rescale x
    xlim_min = 1979.5

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
bloom_timing_df = bloom_timing.merge(suchy_bloom_timing_df, left_on='Year', right_on='Year', how='left')
bloom_timing_df = bloom_timing_df.merge(gower_bloom_timing_df, left_on='Year', right_on='Year_Gower', how='left')
if display_allen:
    bloom_timing_df = bloom_timing_df.merge(allen_bloom_timing_df, left_on='Year', right_on='Year_Allen', how='left')
bloom_timing_df['bloom_match'] = bloom_timing_df.apply(lambda row: row['Bloom Early Late_x'] == row['Bloom Early Late_y'], axis=1)
print(bloom_timing_df.columns)
df = bloom_timing_df

# Plotting
if display_allen:
    # fig_width = 17
    # fig_height = 7
    fig_width = 10 # 2024-09
    fig_height = 4
else:
    fig_width = 10
    fig_height = 6
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plotting the Ecospace model data
ax.errorbar(df['Year'], df['Day of Year_x'], yerr=3, fmt='s', color='blue', label='Ecospace Model', capsize=5)
# Plotting the Suchy data
ax.errorbar(df['Year'], df['Day of Year_y'], yerr=8, fmt='s', color='red', label='Satellite Data', capsize=5)

if display_gower:
    ax.errorbar(df['Year'], df['Day of Year_Gower'], yerr=7, fmt='s', color='orange', label='Gower Data', capsize=5)

if display_allen:
    ax.errorbar(df['Year'], df['Day of Year_Allen'], yerr=8, fmt='s', color='green', label='C09 1D Model', capsize=5)

# Adding horizontal dashed lines
ax.axhline(y=68, color='black', linestyle='--')
ax.axhline(y=108, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([xlim_min, x_max], 68, 108, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([xlim_min, xlim_max])
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
# plt.savefig('..//..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_altmethod_' + scenario + '_' + trnsfrm + '.png')
plt.savefig('..//..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_' + str(yr_strt) + '_2016_suchymethod_' + scenario + '_' + trnsfrm + '_rev2.png')
plt.show()


exit()


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
# styled_table.to_html('..//..//figs//' + 'bloom_timing_comparison_' + scenario + '.html')


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

# plt.savefig('..//..//figs//' + scenario + '_' + trnsfrm + '_bloom_timing_comparison.pdf', bbox_inches='tight')
plt.savefig('..//..//figs//' + scenario + '_' + trnsfrm + '_bloom_timing_comparison.png', bbox_inches='tight', dpi=300)

# Display the table
plt.show()




#####################################################
############## VISUALISE AS TIME SERIES #############
# visualise as time series
# Plot the regional average values in a TS
# using horizontal lines to show annual avg and vertical lines to show satellite derived bloom timing
v = list(v_f)[0]
if log_transform:
    outfigname = '..//..//figs//' + 'Ecospace_out_suchy_region_LOG' + v + '_' + scenario + '.png'
else:
    outfigname = '..//..//figs//' + 'Ecospace_out_suchy_region_' + v + '_' + scenario + '.png'

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
