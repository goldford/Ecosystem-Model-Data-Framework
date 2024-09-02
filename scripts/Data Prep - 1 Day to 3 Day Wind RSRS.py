# Created by G Oldford
# Aug 7-8, 2024
# Purpose: process daily RDRS v2.1 wind files (NC) to 3-day fields for testing of forcing Ecospace
#
# Input:
#   1/ NC files of RDRS daily processed on server to daily averages
#
# Output:
#   1/ Ecospace ASC in 3 day blocked MAXIMUM daily means of wind speed (sqrt(u^2 + v^2)) and stress (u^2 + v^2)
#   2/ Ecosim TS, same as above, but averaged Max (this is probably a bit convoluted
#
# note:
#   - not running this on the graham server because it seems to be down.
#   - 2024-08-08 - since maximum winds (sustained for some threshold period) are what is important for disrupting spring
#     PP blooms, I am using the maximum 3-day wind speeds and stress, but these are computed from daily mean
#     wind speeds. I can't do anything about this because the original data are stuck on graham and must be processed
#     there, practically speaking. I also don't know the threshold, so daily is probably okay.

import os
import numpy as np
import netCDF4 as nc
import regex as re
from helpers import is_leap_year, saveASCFile, getDataFrame, buildSortableString
import csv

wind_daily_p = "C:/Users/Greig/Downloads/RDRS_Wind/"
wind_daily_f = "RDRS21_NEMOgrid_wind_daily_{}.nc"
ecospace_out_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecospace"
ecospace_out_speed_subdir = "speed/"
ecospace_out_stress_subdir = "stress"
ecosim_out_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/LTL_MODELS/RDRS forcings/Wind_RDRS/Ecosim"

plume_dir = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap/"
ecospacegrid_dir = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap/"

scaling = "sqrt"
startyear = 1980
endyear = 2002
startMonth = 1
endMonth = 12
out_code = "var"
timestep = "3day"

# ################# ASC Params #################
saveTemplate = "{}_{}_{}.asc" # new: ewename, year, dayofyear
ASCheader = "ncols        93\nnrows        151 \nxllcorner    -123.80739\nyllcorner    48.29471 \ncellsize     0.013497347 \nNODATA_value  0.0"

# ################# MAP CLIP #################
# array offsets for clipping EwE grid out of NEMO grid
# Calculated in "Basemap Converter (Py3).ipynb"
upperleft_row_ewe = 253
upperleft_col_ewe = 39
bottomleft_row_ewe = 102
bottomleft_col_ewe = 39

# ################# PLUME MASK #################
# Create the plume mask
PlumeRegionFilename = os.path.join(plume_dir, "Fraser_Plume_Region.asc")
dfPlumeRegion = getDataFrame(PlumeRegionFilename, "-9999.00000000")
dfPlumeMask = dfPlumeRegion == 1  # Plume is region one

# ################# LAND MASK #################
ecospacegrid_f = os.path.join(ecospacegrid_dir, "ecospacedepthgrid.asc")
df_dep = getDataFrame(ecospacegrid_f, "-9999.00000000")
dfLandMask = df_dep == 0  # mask is where elev =0


################## MAIN LOOP ##################
speed_ecosim_all = []
stress_ecosim_all = []
save_ecospace = True
save_ecosim = True
sigdigfmt_spd = '%0.2f'
sigdigfmt_strs = '%0.1f'
lump_final = True # lump final 121th and 122th with 120th block? (results in 6 hrs = 1 mo)
var_mult = 1

i = 0
for iyear in range(startyear, endyear + 1):

    print(iyear)

    wind_full_path = os.path.join(wind_daily_p, wind_daily_f.format(iyear))
    ds_var = nc.Dataset(wind_full_path)

    u_comp = ds_var['u10m'][:]
    v_comp = ds_var['v10m'][:]
    speed = np.sqrt(u_comp ** 2 + v_comp ** 2)

    # ############### LOOP 3DAY BLOCKS #################
    leapTF = is_leap_year(iyear)
    num_days = 365
    if leapTF:
        num_days = 366

    num_days = num_days // 3  # 3 day 'years'
    if not leapTF:
        num_days += 1  # this helps catch day 364 and 365

    j = 0  # for file naming - the hacked 3day year
    for iday in range(1, num_days + 1):

        # default block of 3 days
        day_strt = (iday - 1) + (j * 2)
        day_end = day_strt + 2
        middle_day = day_strt + 2

        # exceptions for last blocks of days of year -->
        # lump final 121th and 122th with 120th block?
        if (lump_final) and (iday == 120):

            if not leapTF:
                # day_end = day_strt+5
                middle_day = day_strt + 3
            else:
                # day_end = day_strt+6
                middle_day = day_strt + 4
            varday1 = speed[day_strt:, :, :]  # (365, 299, 132)

        else:
            # catch if last block of days (364,365) is only 2 long
            if not leapTF:
                if iday == num_days + 1:
                    day_end = day_strt + 1
                    middle_day = day_strt + 1
            varday1 = speed[day_strt:day_end, :, :]

        speed_3day = np.ma.max(varday1, axis=0)
        speed_3day = speed_3day * var_mult

        # Mask out invalid and infinite values
        masked_speed_3day = np.ma.masked_invalid(speed_3day)
        masked_speed_3day = np.ma.masked_where(np.isinf(masked_speed_3day), masked_speed_3day)

        # Clip the values to avoid overflow before squaring
        clipped_speed_3day = np.clip(masked_speed_3day, a_min=None,
                                     a_max=np.sqrt(np.finfo(masked_speed_3day.dtype).max))
        stress_3day = np.square(clipped_speed_3day)
        stress_3day = np.ma.masked_where(np.isinf(stress_3day), stress_3day)

        # save to ECOSPACE forcing file
        if save_ecospace:

            # wind speed
            savepath = os.path.join(ecospace_out_p, ecospace_out_speed_subdir)
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)

            savefn = saveTemplate.format("RDRS_windspeed10m", iyear, buildSortableString(middle_day, 2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, speed_3day, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe, sigdigfmt_spd,
                        ASCheader, dfPlumeMask, dfLandMask)

            # wind stress
            savepath = os.path.join(ecospace_out_p, ecospace_out_stress_subdir)
            if (os.path.isdir(savepath) != True):
                os.mkdir(savepath)

            savefn = saveTemplate.format("RDRS_windstress10m", iyear, buildSortableString(middle_day, 2))
            savefullpath = os.path.join(savepath, savefn)
            print("Saving file " + savefullpath)
            saveASCFile(savefullpath, stress_3day, bottomleft_row_ewe, upperleft_row_ewe, upperleft_col_ewe, sigdigfmt_strs,
                        ASCheader, dfPlumeMask, dfLandMask)

        if save_ecosim:
            # ############# ECOSIM CSV Forcing files #############
            # ecosim forcing grid requires monthly TS for long-term forcing
            # so here we avg across map and then expand, inserting

            var_ecosim = speed_3day.flatten()
            var_ecosim = var_ecosim * var_mult
            var_ecosim = var_ecosim[var_ecosim != 0]  # drop zeros for averaging
            var_ecosim = np.mean(var_ecosim)
            # Extract the rounding decimal places from format string
            match = re.search(r'\.(\d+)f', sigdigfmt_spd)
            if match:
                precision = int(match.group(1))
                var_ecosim = np.round(var_ecosim, precision)
            else:
                print("Issue with finding sig digits to round for ecosim out.")
            speed_ecosim_all.append([iyear, middle_day, var_ecosim])

            var_ecosim = stress_3day.flatten()
            var_ecosim = var_ecosim * var_mult
            var_ecosim = var_ecosim[var_ecosim != 0]  # drop zeros for averaging
            var_ecosim = np.mean(var_ecosim)
            # Extract the rounding decimal places from format string
            match = re.search(r'\.(\d+)f', sigdigfmt_spd)
            if match:
                precision = int(match.group(1))
                var_ecosim = np.round(var_ecosim, precision)
            else:
                print("Issue with finding sig digits to round for ecosim out.")
            stress_ecosim_all.append([iyear, middle_day, var_ecosim])

        if lump_final and (iday == 120):
            break

        j += 1
    i += 1
if save_ecosim:

    # wind speed
    # add index corresponding to hacked 3day year
    var_ecosim_all_idx = [[index + 1] + sublist for index, sublist in enumerate(speed_ecosim_all)]
    # expand the array since ecosim wants values per month not year in enviro forcing grid
    var_ecosim_all_expnd = [sublist for sublist in var_ecosim_all_idx for _ in range(12)]
    # index by time step (fake 3day 'months' = 6hr)
    var_ecosim_all_expnd_idx = [[i + 1] + row for i, row in enumerate(var_ecosim_all_expnd)]
    column_names = ['threeday_yrmo', 'threeday_yr', 'year', 'dayofyear', "RDRS_windspeed10m"]

    savepath = os.path.join(ecosim_out_p)
    if (os.path.isdir(savepath) != True):
        os.mkdir(savepath)

    savefn = "RDRS_windspeed10m" + "_" + timestep + "_" + str(startyear) + "-" + str(endyear) + ".csv"
    savefullpath = os.path.join(savepath, savefn)
    print("Saving Ecosim file " + savefullpath)

    with open(savefullpath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(column_names)
        writer.writerows(var_ecosim_all_expnd_idx)

    # do wind stress
    var_ecosim_all_idx = [[index + 1] + sublist for index, sublist in enumerate(stress_ecosim_all)]
    var_ecosim_all_expnd = [sublist for sublist in var_ecosim_all_idx for _ in range(12)]
    var_ecosim_all_expnd_idx = [[i + 1] + row for i, row in enumerate(var_ecosim_all_expnd)]
    column_names = ['threeday_yrmo', 'threeday_yr', 'year', 'dayofyear', "RDRS_windstress10m"]

    savepath = os.path.join(ecosim_out_p)
    if (os.path.isdir(savepath) != True):
        os.mkdir(savepath)

    savefn = "RDRS_windstress10m" + "_" + timestep + "_" + str(startyear) + "-" + str(endyear) + ".csv"
    savefullpath = os.path.join(savepath, savefn)
    print("Saving Ecosim file " + savefullpath)

    with open(savefullpath, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(column_names)
        writer.writerows(var_ecosim_all_expnd_idx)
