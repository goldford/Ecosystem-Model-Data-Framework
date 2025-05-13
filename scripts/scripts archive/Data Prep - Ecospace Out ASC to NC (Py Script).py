# Create 2024-31
# G Oldford
# Purpose: convert Ecospace ASC out to NetCDF
#          For use with the 3Day Model

import numpy as np
import xarray as xr
import pandas as pd
import os

import cartopy as cp
import matplotlib.pyplot as plt
import cartopy
from cartopy import crs, feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.patches import Rectangle
import cmocean as cm
import cartopy.crs as ccrs
from helpers import buildSortableString, is_leap_year
from datetime import datetime, timedelta

# paths
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//ECOSPACE 216 2024 - 10yr1yr 2003-2018//asc//"
#ecospace_code = "Scv1-NoMultMix"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v1 - NoMultMixing - 3DAY 1yr10yr//asc//"
#ecospace_code = "Scv2-MultMix"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v2 - MultMixing - 3DAY 1yr10yr//asc//"
#ecospace_code = "Scv3-MultxPAR"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v3 - MixingxPAR - 3DAY 1yr10yr//asc//"
#ecospace_code = "Scv5-PARMixingNut90"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v5 - PARMixingNut90 - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv7-PARMixingNut90Temp"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v7 - PARMixingNut90Temp - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv4_2-MixingxPARLimitZ"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v4 - MixingxPARLimitZ - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv24-PAR_PI_Temp"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v25 - PAR_PI_Temp - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv26-PAR_PI_Temp_mixing"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v26 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv27-PAR_PI_Temp_mixing_habcap"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v27 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv28-PAR_PI_Temp_mixing" # trying non-linear mixing
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v28 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv29-PAR_shallowPI_Temp_mixing" # very shallow PI slope
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v29 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv30-PARenv_shallowPI_Temp_mixing" # very shallow PI slope, PAR is now an enviro driver again (was PP)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v30 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv31-PARenv_shallowPI_Temp_MixingXPARHab" # very shallow PI slope, PAR is now an enviro driver again, MixingxPAR driving hab cap (crazy response from PIC and NAN)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v31 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv32-PARenv_shallowPI_MixingSteep" # very shallow PI slope, PAR is now an enviro driver again, MixingxPAR driving hab cap (crazy response from PIC and NAN)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v32 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv33-PARenv_lessshallowPI_MixingSteep" # very shallow PI slope, PAR is now an enviro driver again, MixingxPAR driving hab cap (crazy response from PIC and NAN)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v33 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv34-PARenv_lessshallowPI_MixingSteep" # very shallow PI slope, PAR enviro driver again, MixingxPAR driving hab cap (crazy response from PIC and NAN)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v34 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv35-PARenv_lessshallowPI_MixingSteep" # this one might be messed up - enviro drivers not working
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v35 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv36-PARenv_lessshallowPI_MixingSteep" # # this one might be messed up - enviro drivers not working
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v36 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv37-PARenv_lessshallowPI_MixingSteep_Temp" # this one might be messed up - enviro drivers not working
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v37 - PAR_PI_Temp_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv38-PARenv_PI_Temp_Wind" # introduces wind, drivers working, 2003 poor but other yrs generally better
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v38 - PAR_PI_Temp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv39-PARenv_PI_Temp_Wind" # like 38 but sensitivity test of slope of wind response
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v39 - PAR_PI_Temp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv40-PARenv_PI_Temp_Wind" # like 39 but with mixing response for simulating habitat volume
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v40 - PAR_PI_Temp_Wind_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv41-PARenv_PI_Temp_Wind" # should be identical to 39 or 40 (check which) but getting different results
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v41 - PAR_PI_Temp_Wind_mixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv42-PARenv_PI_Temp_Wind_Mixing" # should be identical to 39 or 40 (check which) but getting different results
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v42 - PAR_PI_Temp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv43-All_Groups_Temp" # should be identical to 39 or 40 (check which) but getting different results
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v43 - All Groups w Temp - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv45-All_Groups_Temp" # should be identical to 39
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v45 - All Groups w Temp - 3DAY 1yr10yr//asc//"
ecospace_code = "Scv80_1-All_Groups_20250501" # should be identical to 39
path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_v12_ewe6_7_19295_SC51_04 - DEBUG5//Sc216 v80_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"

path_out = "..//data//ecospace_out//"
nemo_ewe_csv = "..//data//basemap//Ecospace_grid_20210208_rowscols.csv"

# mxng_p = "NEMO_prepped_as_ASC/{var}/"
# tmp_p = "NEMO_prepped_as_ASC/{var}/"
# li_p = "ECOSPACE_in_PAR3_Sal10m/{var}/"
# k_p = "ECOSPACE_in_PAR3_Sal10m/RUN203_{var}/"
# sal_p = "NEMO_prepped_as_ASC/{var}/"
# path_data2 = "../data/forcing/"
# li_p2 = "RDRS_light_monthly_ASC/{var}/"
# wi_p = "RDRS_wind_monthly_ASC/{var}/"

# asc to nc
def asc_to_nc_3day(v_f, outfilename, nemo_ewe_csv,
                   rows=151, cols=93, skiprows=6,
                   yr_strt=2003, yr_end=2003,
                   mo_strt=1, da_strt=2,
                   mo_end=12, da_end=30):

    df = pd.read_csv(nemo_ewe_csv)

    # get the day of year in 3DAY BLOCKS
    # while dumping the remainder (blocks 121,122) b/c this is how 3D model driver data prepped
    # ie 1 year = '10 years' thus 120 time steps per real year, 5-6 days remainder
    days_of_year = range(2, 360, 3)
    date_list = []
    months = list(range(1, 12))
    pd_timestamps = []
    time_step_model = []
    i = 0
    for yr in range(yr_strt, yr_end+1):
        for doy in days_of_year:
            date1 = datetime(yr, 1, 1) + timedelta(days=doy - 1)
            month1 = date1.month
            if yr == yr_end and (month1 == mo_end + 1):
                break
            day_of_month = date1.day
            date_list.append([yr, month1, day_of_month, doy])
            pd_timestamps.append(pd.Timestamp(yr, month1, day_of_month))
            time_step_model.append(buildSortableString(i+1, 5))
            i += 1
        if yr == yr_end:
            break

    # https://pandas.pydata.org/docs/user_guide/timeseries.html#timeseries-offset-aliases
    # time = pd.date_range(start='{yr_strt}-{mo_strt}-{da_strt}'.format(yr_strt=yr_strt, mo_strt=mo_strt, da_strt=da_strt),
    #                               end='{yr_end}-{mo_end}-{da_end}'.format(yr_end=yr_end, mo_end=mo_end, da_end=da_end),
    #                               freq='3D')
    # for t in range(1,len(pd_timestamps)+1):
    #     time_step_model.append(buildSortableString(t,5))

    # Create an empty xarray dataset
    ds = xr.Dataset(
        coords={
            'time': pd_timestamps,
            'row': range(1, rows + 1),
            'col': range(1, cols + 1),
            'lat': (('row', 'col'), np.full((rows, cols), np.nan)),
            'lon': (('row', 'col'), np.full((rows, cols), np.nan)),
            'depth': (('row', 'col'), np.full((rows, cols), np.nan)),
            'EWE_col': (('row', 'col'), np.full((rows, cols), np.nan)),
            'EWE_row': (('row', 'col'), np.full((rows, cols), np.nan)),
            'NEMO_col': (('row', 'col'), np.full((rows, cols), np.nan)),
            'NEMO_row': (('row', 'col'), np.full((rows, cols), np.nan)),
        },
        attrs={'description': 'dataset of monthly ASC files'}
    )

    # Populate the dataset with lat, lon, depth, EWE_col, EWE_row from the CSV file
    for _, row in df.iterrows():
        ewe_row = int(row['EWE_row']) - 1
        ewe_col = int(row['EWE_col']) - 1
        ds['lat'][ewe_row, ewe_col] = row['lat']
        ds['lon'][ewe_row, ewe_col] = row['lon']
        ds['depth'][ewe_row, ewe_col] = row['depth']
        ds['EWE_col'][ewe_row, ewe_col] = row['EWE_col']
        ds['EWE_row'][ewe_row, ewe_col] = row['EWE_row']
        ds['NEMO_col'][ewe_row, ewe_col] = row['NEMO_col']
        ds['NEMO_row'][ewe_row, ewe_col] = row['NEMO_row']

    # create empty variable with correct shape
    for v in v_f:
        ds[v] = xr.DataArray(
            np.nan * np.zeros((len(pd_timestamps), rows, cols)),
            dims=('time', 'row', 'col'),
            attrs={'description': f'{v} data'}
        )

    # load these ASCs into a nice xarray dataset
    for v in v_f:
        attribute = v_f[v]
        print(attribute)
        for t in range(0, len(time_step_model)):
            f_n = v_f[v].format(time_step_model[t])
            ti = pd_timestamps[t]
            yr = ti.year
            mo = ti.month
            da = ti.day

            with open(f_n) as f:
                data = np.loadtxt(f, skiprows=skiprows)

                # homogenize what nans are
                data[data == -9999.0] = ['nan']
                data[data == 0.0] = ['nan']

                # fix issue with bottom left area in map
                data[140:, :15] = ['nan']

                ds[f'{v}'.format(var=v)].loc[{'time': f'{yr}-{mo}-{da}'.format(year=yr, month=mo, day=da)}] = xr.DataArray(
                    data,
                    dims=('row', 'col'),
                    attrs={'description': f'{v} data for year {yr} month {mo} day {da}'.format(var=v, year=yr, month=mo, day=da)}
                )

    ######## SAVE TO NETCDF ########
    # Create an encoding dictionary for compression
    encoding = {var: {'zlib': True, 'complevel': 5} for var in ds.data_vars}

    # Write the xarray dataset to a compressed NetCDF file
    output_file = outfilename
    ds.to_netcdf(output_file, format='NETCDF4', encoding=encoding)

    # Check if the original and loaded datasets are identical
    ds_loaded = xr.open_dataset(output_file)
    datasets_identical = ds.equals(ds_loaded)
    print("Datasets are identical:", datasets_identical)

    return ds, datasets_identical


v_f = {"NK1-COH": path_ecospace_out + "EcospaceMapBiomass-NK1-COH-{}.asc",
       "NK2-CHI": path_ecospace_out + "EcospaceMapBiomass-NK2-CHI-{}.asc",
       "NK3-FOR": path_ecospace_out + "EcospaceMapBiomass-NK3-FOR-{}.asc",
       "ZF1-ICT": path_ecospace_out + "EcospaceMapBiomass-ZF1-ICT-{}.asc",
       "ZC1-EUP": path_ecospace_out + "EcospaceMapBiomass-ZC1-EUP-{}.asc",
       "ZC2-AMP": path_ecospace_out + "EcospaceMapBiomass-ZC2-AMP-{}.asc",
       "ZC3-DEC": path_ecospace_out + "EcospaceMapBiomass-ZC3-DEC-{}.asc",
       "ZC4-CLG": path_ecospace_out + "EcospaceMapBiomass-ZC4-CLG-{}.asc",
       "ZC5-CSM": path_ecospace_out + "EcospaceMapBiomass-ZC5-CSM-{}.asc",
       "ZS1-JEL": path_ecospace_out + "EcospaceMapBiomass-ZS1-JEL-{}.asc",
       "ZS2-CTH": path_ecospace_out + "EcospaceMapBiomass-ZS2-CTH-{}.asc",
       "ZS3-CHA": path_ecospace_out + "EcospaceMapBiomass-ZS3-CHA-{}.asc",
       "ZS4-LAR": path_ecospace_out + "EcospaceMapBiomass-ZS4-LAR-{}.asc",
       "PZ1-CIL": path_ecospace_out + "EcospaceMapBiomass-PZ1-CIL-{}.asc",
       "PZ2-DIN": path_ecospace_out + "EcospaceMapBiomass-PZ2-DIN-{}.asc",
       "PZ3-HNF": path_ecospace_out + "EcospaceMapBiomass-PZ3-HNF-{}.asc",
       "PP1-DIA": path_ecospace_out + "EcospaceMapBiomass-PP1-DIA-{}.asc",
       "PP2-NAN": path_ecospace_out + "EcospaceMapBiomass-PP2-NAN-{}.asc",
       "PP3-PIC": path_ecospace_out + "EcospaceMapBiomass-PP3-PIC-{}.asc",
       "BA1-BAC": path_ecospace_out + "EcospaceMapBiomass-BA1-BAC-{}.asc" #time step format eg: 00620,
      }

rows = 151
cols = 93
skiprows = 6 # header

# corresponding lats and lons to centres of grid cells
# read the data into an xarray object which is versatile
yr_strt = 1978
yr_end = 2018
mo_strt = 1
da_strt = 2
mo_end = 12
da_end = 30

out_filename = ecospace_code + "_" + str(yr_strt) + "-" + str(yr_end) + ".nc"
out_filename = os.path.join(path_out, out_filename)


original_ds, are_identical = asc_to_nc_3day(v_f, out_filename, nemo_ewe_csv,
                                            rows, cols, skiprows,
                                            yr_strt, yr_end, mo_strt, da_strt,
                                            mo_end, da_end)
print("Crosscheck ASC in = NC data out:", are_identical)