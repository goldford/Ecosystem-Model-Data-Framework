# Created by G Oldford
# Aug 29, 2024
# Purpose: compare bloom timing from ecospace to observations
#
# Input:
#   1/ NC files of Ecospace out processed from ASC by script #1
#   2/ YAML outlining 'region 2' corresponding to Suchy et al 2021 region (satellite)
#   3/ bloom timing dates from central SoG from empirical obs and models
#
# Output:
#   1/ saved bloom timing df's so it does not have to be computed each time (csv's)
#   2/ figs of estimated bloom timing from Ecospace compared to observations
# note:
#  2025-05-05 - added exclude_decjan to helpers function to make annual average more comparable to satellite method
#  2025-05-07 - be careful with this script! It does not update the bloom timing df's unless you flip switches
#  2025-05-08 - this script is currently a total mess

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
# scenario = 'FULLKEY_SC51_5'
# ecospace_code = "FULLKEY_Scv51_5-PAR_PI_AllPPTemp_Wind"
# scenario = 'FULLKEY_SC80_1'
# ecospace_code = "Scv80_1-All_Groups_20250501"
# scenario = 'FULLKEY_SC82_1' # depth responses for eup and dec removed
# ecospace_code = "Scv82_1-All_Groups_20250506"
# scenario = 'FULLKEY_SC83_1' # now with many env responses for zoop
# ecospace_code = "Scv83_1-All_Groups_20250506"
# scenario = 'FULLKEY_SC84_1' # now with many env responses for zoop AND nutrients back up from 90 to 96
# ecospace_code = "Scv84_1-All_Groups_20250506"
scenario = 'FULLKEY_SC85_1' # now with many env responses for zoop AND nutrients back up from 90 to 96
ecospace_code = "Scv85_1-All_Groups_20250506"

filenm_yr_strt = 1978

log_transform = True # following suchy, make true
# log transforming follows Suchy (though theirs is chl)
mean_or_median = "median" # use mean or median for def of bloom?
average_from = "annual"   # 'annual' to follow suchy, set bloom def based on annual avg (alt: "all")
# TECHNICALLY not annual in suchy
thrshld_fctr = 1.05 # threshold above avg, suchy used 1.05
sub_thrshold_fctr = 0.7 # the secondary threshold, if first met, the test of 1/2 following 2 weeks (does it stay in bloom)
min_y_tick = 38 # 38 in previous plots
display_gower = False # true in most previous plots
display_allen = False # display allen 1D model results comparison in plot?
display_ssc_s3 = False

# compute ecospace bloom date at CO9 1D location (no need to run again if already done)
compute_ecospace_1d = True
# compute bloom date for Suchy et al 2022 SSoG region (unless already done)
compute_ecospace_SSoG = True
display_allen = False
exclude_decjan = False

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

# Suchy et al 2022; first indices below are for year 2003
doy_suchy = [100, 68, 50, 83, 115, 115, 100, #2003 - 2009
             100, 92, 88, 92, 92, 55, 77]

bloom_early_suchy = 68 # day is inclusive
bloom_late_suchy = 108 # day is inclusive

# Gower & King, 2018; > 5 mg Chla; 2000 - 2016
doy_gower = [83, 68, 76, 85, 71, 56, 84, 49, 71, 58, #2000 -2009
             90, 97, 91, 69, 94, 66, 86] #2010 - 2016

# C09 Collins et al 2009, Allen & Wolfe, 2013; Allen & Latournell, 2024 (SOPO)
# from 1D hindcast, peak bloom dates, station s3
# first indices below are 1980 (+/- 4 days should be used for error bars)
doy_allen = [94, 78, 81, 82, 87, 76, 75, 87, 99, 76, #1980 - 1989
             81, 78, 77, 55, 86, 86, 67, 87, 77, 104, #1990-1999
             66, 70, 92, 86, 81, 62, 88, 100 ,90, 97, #2000-2009
             104, #1980-2010 (from 2013 pub)
             103, 98, 88, 86, 77, 88, 104, 77, #2011-2018 # from SOPO 2024 rep.
             ]

std_dev_allen = np.std(doy_allen)
mean_dev_allen = np.mean(doy_allen)
# bloom_early_C09 = 71 # from pub
# bloom_late_C09 = 96 # from pub
bloom_early_C09 = mean_dev_allen - std_dev_allen # day is inclusive for early
bloom_late_C09 = mean_dev_allen + std_dev_allen # day is inclusive for late

# Salish Sea Cast 3D 'S3' bloom timing 2007 - 2020
# https://salishsea-meopar-docs.readthedocs.io/en/latest/bloom_timing/bloom_notebooks/201905EnvironmentalDrivers_S3.html
# roughly digitised by GO 2025-03-26; 201905 model v
doy_allen_3D_metric3 = [105, 98, 101, 109, 113, 105,
                        89, 88, 74, 92, 105, 82, 81, 87] # metric 3


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
    'Year_Gower': range(2000, 2016 + 1),
    'Bloom Date_Gower': dates_gower,
    'Day of Year_Gower': doy_gower,
    "Bloom Early Late Gower": ['na', 'na', 'na', 'na', 'na', 'na', 'na', 'na', 'na', 'na',
                               'na', 'na', 'na', 'na', 'na', 'na', 'na']
})

dates_allen = []
for i, doy in enumerate(doy_allen):
    year = 1980 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dates_allen.append(date)

allen_bloom_timing_df = pd.DataFrame({
    'Year_Allen': range(1980, 2018 + 1),
    'Bloom Date_Allen': dates_allen,
    'Day of Year_Allen': doy_allen,
    "Bloom Early Late Allen": np.full(len(range(1980,2018+1)), 'na')
})

def classify_bloom(doy):
    if doy + 4 <= bloom_early_C09:
        return "early"
    elif doy - 4 >= bloom_late_C09:
        return "late"
    elif (doy + 4 <= bloom_late_C09 - 1) and (doy - 4 >= bloom_early_C09 + 1):
        return "avg"
    else:
        return "cusp"
allen_bloom_timing_df["Bloom Early Late Allen"] = allen_bloom_timing_df["Day of Year_Allen"].apply(classify_bloom)

# added by GO - 2025-03-26
dates_allen_3D = []
for i, doy in enumerate(doy_allen_3D_metric3):
    year = 2007 + i
    date = datetime(year, 1, 1) + timedelta(days=doy - 1)
    dates_allen_3D.append(date)

allen3D_bloom_timing_df = pd.DataFrame({
    'Year_Allen3D': range(2007, 2020 + 1),
    'Bloom Date_Allen3D': dates_allen_3D,
    'Day of Year_Allen3D': doy_allen_3D_metric3,
    "Bloom Early Late Allen3D": np.full(len(range(2007, 2020 + 1)), 'na')
})


######## OPEN ECOSPACE OUT NC FILE ##########
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
mask_satel = mask_ds['mask'].values
# masked_data = ds.where(mask)
# masked_pp1_dia = masked_data['PP1-DIA']
# print(np.nanmean(masked_pp1_dia))

#############################################
############## MASK ALLEN 1D ################

create_allen_mask = False
if create_allen_mask:

    # station S3
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

# hardcoded the results from above
allen1d_idx_and_depth = [52, 100, 370.68523301071156] # row, col, depth of allen 1d model location
allen3d_idx_and_depth = [52, 100, 370.68523301071156] # row, col, depth of allen 3d 'S3' timing extract - GO 2025-03-26


#############################################
# get the ecospace data for a single model point
annual_means_allen = [] # allen refers to location of allen 1D model
annual_medians_allen = []
ts_reg_means_allen = []
ts_reg_medians_allen = []
timestamps_allen = []


exclude_juntojan = False # keep it in to try to be more in line with C09 rather than satellite data
use_satellite_mask = False
if use_satellite_mask:
    satmask_suffx = "_satmask_"
else:
    satmask_suffx = "_C09loc_"
if compute_ecospace_1d:
    v = list(v_f)[0] # hard coding to DIA
    print('ecospace model data at allen 1D model location')
    for year in range(yr_strt, yr_end + 1):
        ############################################
        # pull ecospace results for comparison to C09

        print(year)
        year_data = ds.sel(time=str(year))

        if use_satellite_mask:
            masked_ds = year_data.where(mask_satel)
            # masked_pp1_dia = masked_data['PP1-DIA']
        else:
            masked_ds = year_data.sel(row=allen1d_idx_and_depth[1],
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

        # Select time period
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
                                                                      mean_or_median=mean_or_median,
                                                                      exclude_juntojan=exclude_juntojan,
                                                                      bloom_early=bloom_early_C09,
                                                                      bloom_late=bloom_late_C09
                                                                      )
    ecospace_bloom_timing_df = pd.DataFrame({
        'Year': range(yr_strt, yr_end + 1),
        'Bloom Date': bloom_dates,
        'Day of Year': bloom_days_of_year,
        "Bloom Early Late": bloom_earlylate
    })

    # Export
    ecospace_bloom_timing_df.to_csv('..//..//data//evaluation//ecospace_bloom_timing_CSoG_1D_' +
                                    satmask_suffx +
                                    scenario +
                                    '.csv', index=False)

else:
    bloom_p = "..//..//data//evaluation//"
    ecospace_bloom_timing_df = pd.read_csv(os.path.join(bloom_p, 'ecospace_bloom_timing_CSoG_1D_' +
                                                        satmask_suffx +
                                                        scenario +
                                                        '.csv'))




################################################
############### DOY PLOT 1 - long ###############
### Shows Allen's 1D versus Ecospace output ###
### Bloom timing based on Allen 1D location ###
print("DoY bloom comparison plot 1")
# Create a DataFrame to display the results

xlim_min = 1979.5
xlim_max = 2018.5

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
jitter_x = 0.0
ax.errorbar(df['Year']-jitter_x, df['Day of Year'], yerr=1.5, fmt='s',
            color='blue', label='Ecospace',
            capsize=0, marker='o', elinewidth=1,
            markersize=4)

ax.errorbar(df['Year']+jitter_x, df['Day of Year_Allen'], yerr=4, fmt='s',
            color='green', label='C09 Model',
            capsize=0, marker='s', elinewidth=1,
            markersize=4)

# Adding horizontal dashed lines
ax.axhline(y=bloom_early_C09, color='black', linestyle='--')
ax.axhline(y=bloom_late_C09, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([xlim_min, x_max], bloom_early_C09, bloom_late_C09, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([xlim_min, xlim_max])
y_min, y_max = ax.get_ylim()
ax.set_ylim([min_y_tick, y_max])
ax.set_xlabel('Year')
ax.set_ylabel('Day of Year')
# ax.set_title('Diatom Bloom Timing Comparison')

highlight_years = False
# for satellite data years, see if ecospace and sat agree (overrides above)
years_to_check = [i for i in range(1980, 2018)]
if highlight_years:

    df_agreement = df[['Year', 'Day of Year', 'Day of Year_Allen']]
    df_agreement.columns = ['Year', 'Ecospace', 'Allen']
    # do allen and ecospace agree?
    df_agreement['agree_tf'] = (
            (df_agreement['Ecospace'] + 1.5 >= df_agreement['Allen'] - 4) &
            (df_agreement['Ecospace'] - 1.5 <= df_agreement['Allen'] + 4)
    )

    df_agreement.loc[df_agreement['Year'].isin(years_to_check), 'agree_tf'] = (
            (df_agreement['Ecospace'] + 1.5 >= df_agreement['Allen'] - 4) &
            (df_agreement['Ecospace'] - 1.5 <= df_agreement['Allen'] + 4)
    )

    # Iterate over the years where agree_tf is False
    for _, row in df_agreement[df_agreement['agree_tf'] == False].iterrows():
        year = row['Year']
        # Shade between year - 0.2 and year + 0.2
        ax.fill_betweenx(ax.get_ylim(), year - 0.2, year + 0.2, color='red', alpha=0.1)


ax.legend()

if log_transform:
    trnsfrm = "logDIA"
else:
    trnsfrm = ""

# plt.title(scenario + " " + trnsfrm)
plt.title("")

outfig = '..//..//figs//' + 'smBloomTiming_ecospace_vs_C09' + str(yr_strt) + '_' + \
         str(yr_end) + 'suchymethod_' + scenario + '_' + trnsfrm + '_rev4.png'
# outfig = '..//..//figs//' + 'smBloomTiming_justAllen1S_' + str(yr_strt) + '_' + \
#          str(yr_end) + 'alttosuchymethod_' + scenario + '_' + trnsfrm + '_rev4.png'
plt.savefig(outfig)
plt.show()

print("done")
print(outfig)

############## STOP LIGHT PLOT #################
# Create a new DataFrame for the table
print(df)
table_df = df[['Year', 'Bloom Early Late', 'Bloom Early Late Allen']]
table_df.columns = ['Year', 'Ecospace', 'Allen']
table_df = table_df.loc[table_df['Year'].isin(years_to_check)]

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


fig, ax = plt.subplots(figsize=(3, 10))  # Adjust the figure size as needed

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
        ['lightgreen' if row['Ecospace'] == row['Allen'] else 'lightcoral']*len(row)
        for _, row in table_df.iterrows()
    ]
)

# Set font size
table.set_fontsize(9)
table.scale(1.2, 1.2)

# plt.savefig('..//..//figs//' + scenario + '_' + trnsfrm + '_bloom_timing_comparison.pdf', bbox_inches='tight')
plt.savefig('..//..//figs//' + scenario + '_' + trnsfrm + '_C09_bloom_timing_stoplight.png', bbox_inches='tight', dpi=300)

# Display the table
plt.show()

#################################################
############### DOY PLOT 2 ###############
###### Shows SSC 3D versus Ecospace output ######
##### Bloom timing based on SSC S3 location #####
print("DoY bloom comparison plot 2")
# Create a DataFrame to display the results

xlim_min = 2006.5
xlim_max = 2018.5

# join the model and obs bloom info
bloom_timing_df = ecospace_bloom_timing_df.merge(allen3D_bloom_timing_df, left_on='Year', right_on='Year_Allen3D', how='left')
print(bloom_timing_df.columns)
bloom_timing_df['bloom_match'] = bloom_timing_df.apply(lambda row: row['Bloom Early Late'] == row['Bloom Early Late Allen3D'], axis=1)
print(bloom_timing_df.columns)
df = bloom_timing_df

# Plotting
fig_width = 6
fig_height = 6
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plotting the Ecospace model data
ax.errorbar(df['Year']-0.1, df['Day of Year'], yerr=1.5, fmt='s',
            color='blue', label='Ecospace Model',
            capsize=0, marker='o', elinewidth=1,
            markersize=6)

ax.errorbar(df['Year']+0.1, df['Day of Year_Allen3D'], yerr=1, fmt='s',
            color='green', label='SSC v201905 at S3 (metric 3)',
            capsize=0, marker='s', elinewidth=1,
            markersize=6)

# Adding horizontal dashed lines
ax.axhline(y=bloom_early_suchy, color='black', linestyle='--')
ax.axhline(y=bloom_late_suchy, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([xlim_min, x_max], bloom_early_suchy, bloom_late_suchy, color='lightgrey', alpha=0.5, zorder=0)
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
outfig = '..//..//figs//' + 'smBloomTiming_Ecospace_vs_SSC3D_S3_' + str(yr_strt) + '_' + \
         str(yr_end) + 'suchymetric3_' + scenario + '_' + trnsfrm + '_rev1.png'
plt.savefig(outfig)
plt.show()

print("done")
print(outfig)



##############################################
########### SUCHY REGIONAL ANALYSIS ##########
########### plot 3 ##########################

if exclude_decjan:
    include_months = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11]# months to include during averaging (Suchy only had data Feb - Nov)
else:
    include_months = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]# months to include during averaging (Suchy only had data Feb - Nov)

create_maps = False

# Loop over the years and select data based on the 'time' coordinate
annual_means = []
annual_medians = []
ts_reg_means = []
ts_reg_medians = []
timestamps = []

# TO DO: the map gen code is embedded below - should be in it's own section
if compute_ecospace_SSoG:
    v = list(v_f)[0] # hard coding to DIA

    for year in range(yr_strt, yr_end + 1):
        ############################################
        # pull model results for 'region' (mask corresponding to suchy)
        print(year)
        #year_data = ds.sel(time=str(year)) # why does this work?? - GO 2025
        year_data = ds.sel(time=str(year)).where(ds.time.dt.month.isin(include_months), drop=True)

        masked_ds = year_data.where(mask_satel) # apply 'regional' mask
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
                                                                      mean_or_median=mean_or_median,
                                                                      exclude_juntojan=False,
                                                                      bloom_early=bloom_early_suchy,
                                                                      bloom_late=bloom_late_suchy
                                                                      )
    ecospace_bloom_timing_df = pd.DataFrame({
        'Year': range(yr_strt, yr_end + 1),
        'Bloom Date': bloom_dates,
        'Day of Year': bloom_days_of_year,
        "Bloom Early Late": bloom_earlylate
    })

    # save to file
    if exclude_decjan:
        decjan = "_noDecJan_"
    else:
        decjan = "_"
    ecospace_bloom_timing_df.to_csv('..//..//data//evaluation//ecospace_bloom_timing_SSoG_' +
                                    decjan +
                                    scenario +
                                    '.csv', index=False)

else:
    if exclude_decjan:
        decjan = "_noDecJan_"
    else:
        decjan = "_"
    bloom_p = "..//..//data//evaluation//"
    ecospace_bloom_timing_df = pd.read_csv(os.path.join(bloom_p, 'ecospace_bloom_timing_SSoG_' +
                                                        decjan +
                                                        scenario +
                                                        '.csv')) #alt ecospace_bloom_timing_SSoG_usingFebtoNov.csv



################################################
############### DOY PLOT 3 #####################
### Shows Suchy et al and optionally Allen's ###
### With bloom timing based on suchy region ####

print("DoY bloom comparison plot")
# Create a DataFrame to display the results
min_y_tick = 30 # 38 in previous plots
xlim_min = 2002.5
xlim_max = 2016.5
jitter_x = 0.0

# display comparison with allen 1D C09 model at S3?
if display_allen: # rescale x
    xlim_min = 1979.5
    xlim_max = 2018.5

if display_ssc_s3:
    xlim_max = 2016.5
    xlim_max = 2018.5


# join the model and obs bloom info
bloom_timing = ecospace_bloom_timing_df
# suchy satellite
bloom_timing_df = bloom_timing.merge(suchy_bloom_timing_df, left_on='Year', right_on='Year', how='left')

# gower satellite
bloom_timing_df = bloom_timing_df.merge(gower_bloom_timing_df, left_on='Year', right_on='Year_Gower', how='left')

# 1D C09 model @ S3
if display_allen:
    bloom_timing_df = bloom_timing_df.merge(allen_bloom_timing_df, left_on='Year', right_on='Year_Allen', how='left')

# 3D SSC model @ S3
if display_ssc_s3:
    bloom_timing_df = bloom_timing_df.merge(allen3D_bloom_timing_df, left_on='Year', right_on='Year_Allen3D', how='left')


bloom_timing_df['bloom_match'] = bloom_timing_df.apply(lambda row: row['Bloom Early Late_x'] == row['Bloom Early Late_y'], axis=1)
print(bloom_timing_df.columns)

df = bloom_timing_df

# Plotting
if display_allen | display_ssc_s3:
    # fig_width = 17
    # fig_height = 7
    fig_width = 13 # 2024-09
    fig_height = 5
else:
    fig_width = 7
    fig_height = 5
fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# Plotting the Ecospace model data
# note that x can be offset by adding or subtracting from year
ax.errorbar(df['Year'], df['Day of Year_x'], yerr=1.5, fmt='s',
            color='black', label='Ecospace Model', #capsize=5,
            capsize=0, marker='o', elinewidth=1,
            markersize=4
            )

# ax.errorbar(df['Year'], df['Day of Year_x'], yerr=1.5, fmt='s',
#             color='black', label='Ecospace Model (' + scenario + ')', #capsize=5,
#             capsize=0, marker='o', elinewidth=1,
#             markersize=4
#             )

# Plotting the Suchy data
ax.errorbar(df['Year']+jitter_x, df['Day of Year_y'], yerr=4, fmt='s',
            color='blue', label='Satellite Data',
            capsize=0, marker='s', elinewidth=1,
            markersize=4)

# if display_gower:
#     ax.errorbar(df['Year'], df['Day of Year_Gower'], yerr=3.5, fmt='s',
#                 color='orange', label='Gower Data',
#                 capsize=0, marker='^', elinewidth=1,
#                 markersize=5)
#
if display_allen:
    ax.errorbar(df['Year']-jitter_x, df['Day of Year_Allen'], yerr=4, fmt='s',
                color='green', label='C09',
                capsize=0, marker='^', elinewidth=1,
                markersize=5)
#
if display_ssc_s3:
    ax.errorbar(df['Year']-jitter_x, df['Day of Year_Allen3D'], yerr=1.5, fmt='s',
                color='red', label='SSCv201905 - S3',
                capsize=0, marker='^', elinewidth=1,
                markersize=5)


# Adding horizontal dashed lines
ax.axhline(y=bloom_early_suchy, color='black', linestyle='--')
ax.axhline(y=bloom_late_suchy, color='black', linestyle='--')

# Setting the x-limits to ensure the grey rectangle fills the entire plot width
x_min, x_max = ax.get_xlim()

# Adding a rectangle filled with light grey that spans the entire width of the plot
ax.fill_between([xlim_min, x_max], bloom_early_suchy, bloom_late_suchy, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([xlim_min, xlim_max])
y_min, y_max = ax.get_ylim()
ax.set_ylim([min_y_tick, y_max])
ax.set_xlabel('Year')
ax.set_ylabel('Day of Year')
# ax.set_title('Diatom Bloom Timing Comparison')

# highlight years with disagreement?
highlight_years = True
if highlight_years:

    if display_allen:
        df_agreement = df[['Year', 'Day of Year_x', 'Day of Year_y', 'Day of Year_Allen']]
        df_agreement.columns = ['Year', 'Ecospace', 'Suchy', 'Allen']
        # do allen and ecospace agree?
        df_agreement['agree_tf'] = (
            (df_agreement['Ecospace'] + 1.5 >= df_agreement['Allen'] - 4) &
            (df_agreement['Ecospace'] - 1.5 <= df_agreement['Allen'] + 4)
            )
    else:
        df_agreement = df[['Year', 'Day of Year_x', 'Day of Year_y']]
        df_agreement.columns = ['Year', 'Ecospace', 'Suchy']


    # for satellite data years, see if ecospace and sat agree (overrides above)
    years_to_check = [2003, 2004, 2005, 2006, 2007, 2008,
                      2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016]
    df_agreement.loc[df_agreement['Year'].isin(years_to_check), 'agree_tf'] = (
            (df_agreement['Ecospace'] + 1.5 >= df_agreement['Suchy'] - 4) &
            (df_agreement['Ecospace'] - 1.5 <= df_agreement['Suchy'] + 4)
    )

    # Iterate over the years where agree_tf is False
    for _, row in df_agreement[df_agreement['agree_tf'] == False].iterrows():
        year = row['Year']
        # Shade between year - 0.2 and year + 0.2
        ax.fill_betweenx(ax.get_ylim(), year - 0.2, year + 0.2, color='red', alpha=0.1)

ax.legend(loc="lower left")

if log_transform:
    trnsfrm = "logDIA"
else:
    trnsfrm = ""

# plt.title(scenario + " " + trnsfrm)
plt.title("")

# plt.savefig('..//..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_altmethod_' + scenario + '_' + trnsfrm + '.png')
if display_allen:
    plot_nm_seg_allen1d = "_vs_C09_"
else:
    plot_nm_seg_allen1d = ""

if display_ssc_s3:
    plot_nm_seg_ssc = "_vs_SSC_"
else:
    plot_nm_seg_ssc = ""

plot_nm_seg = 'BloomTiming_Ecospace_vs_Suchy' + plot_nm_seg_allen1d + plot_nm_seg_ssc

figname= '..//..//figs//' + \
          plot_nm_seg + \
         str(yr_strt) + \
         '_2016_suchymethod_' \
         + scenario + '_' + trnsfrm + '_rev2.png'

plt.savefig(figname)
plt.show()

print("done")
print(figname)


#####################################################
############## VISUALISE AS TIME SERIES #############
# at present only works if compute SSoG above is true
# visualise as time series
# Plot the regional average values in a TS
# using horizontal lines to show annual avg and vertical lines to show satellite derived bloom timing
if compute_ecospace_SSoG:
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

    i = 0
    for year in range(yr_strt, yr_end + 1):
        print(year)
        start_of_year = pd.Timestamp(f'{year}-01-01')
        end_of_year = pd.Timestamp(f'{year}-12-31')
        axes.hlines(annual_means[i] * 1.05, start_of_year, end_of_year, colors='red', linestyle='-', linewidth=1)
        axes.hlines(annual_medians[i] * 1.05, start_of_year, end_of_year, colors='orange', linestyle='-', linewidth=1)
        i += 1

    plt.savefig(outfigname, bbox_inches='tight', dpi=300)
    axes.legend()
    plt.show()


############## STOP LIGHT PLOT #################
# Create a new DataFrame for the table
print(df)
table_df = df[['Year', 'Bloom Early Late_x', 'Bloom Early Late_y']]
table_df.columns = ['Year', 'Ecospace', 'Suchy']
table_df = table_df.loc[table_df['Year'].isin(years_to_check)]

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
exit()



