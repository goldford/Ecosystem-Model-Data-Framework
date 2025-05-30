# Created 2024-08
# By: G Oldford
# Purpose: convert Ecospace ASC out to NetCDF
#          For use with the 3Day Model
# Input: Biomass maps as ASC files at each time step from Ecospace
#
# Output: A single NetCDF file that combines and compresses the inputs
#         Given size issues and github, these are moved manually before subsequent analysis
#         to a sync folder (see shortcuts)
#

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
# ecospace_code = "Scv45-PARPI_Temp_Wind" # should be identical to 39
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v45 - PAR_PI_Temp_Wind - 3DAY 1yr10yr//asc//"

# ecospace_code = "Scv46-PARPI_AllTemp_Wind" # same as 45, 39 but now with temp responses for PP groups other than DIA
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v46 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"

# ecospace_code = "Scv47-PARPI_AllTemp_Wind" # same as 46 but attempted to punish blooms with wind more to compensate for inclusion of temp resp for other PP causing dIA to bloom early in 46
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v47 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv48-PARPI_AllTemp_Wind" # same as 47 but changed to Ecopath biomasses instead of 'habitat adjusted' for init
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v48 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv49-PARPI_AllTemp_Wind" # same as 48 but changed back the wind stress response of DIA to original one used in 46 and prior
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v49 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"

# ecospace_code = "Scv50-RSPI_AllTemp_Wind" # same as 49 but reverted to habitat adjusted biomasses and changing the PI to right shoulder
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v50 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"

# ecospace_code = "Scv51-RSPI_AllTemp_Wind" # same as 50 but introducing data from 2000 - 2002 as 'spin up'
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v51 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv52-RSPI_AllTemp_Wind" # same as 51 but with inital Ecopath B's instead of hab adjusted
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v52 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv53-RSPI_AllTemp_Wind" # same as 52 but with inital Ecopath B's and w/ shallow temp response for DIA (proxy for nutrients)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v53 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv54-RSPI_AllTemp_Wind" # same as 53 but with inital Ecopath B's and w/ shallow temp response for DIA (proxy for nutrients)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v54 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv56-RSPI_AllTemp_Wind" # same as 54 and 55 w/ inital Ecopath B's and w/ shallow temp response for DIA (proxy for nutrients) AND Hab cap driven by light
#                # The hab cap driver seems to bump the variability of growth rates up
#                 # whereas the enviro response method seems only to penalise
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v56 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv58-RSPI_AllTemp_Wind" # same as 57 w/ inital Ecopath B's and w/ shallow temp response for DIA (proxy for nutrients) AND Hab cap driven by light
#                # The hab cap driver seems to bump the variability of growth rates up, removed temp response.
#                # hab cap drivers make it more variable but hope this dampens it
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v58 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv59-RSPI_AllTemp_Wind" # same as 58 w/ inital Ecopath B's and w/ shallow temp response for DIA (proxy for nutrients) AND Hab cap driven by light
#                # The hab cap driver seems to bump the variability of growth rates up, removed temp response.
#                # hab cap drivers make it more variable but hope this dampens it, PB max is 3 or 4 (1000 in others) and nutrients at 95
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v59 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# learned from 59 that PBmax matters, lowering it makes results more variable but possibly to downside
# ecospace_code = "Scv64-RSPI_AllTemp_Wind" # major revamp, ran with the debug version to fix issues
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v15 - DEBUG2//Sc216 v64- PAR_PI_AllPPTemp_Wind - 3DAY//asc//"
# ecospace_code = "Scv56_2-RSPI_AllTemp_Wind" # re run of 56 to cross check results can be repeated
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v56_2 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv56_3-RSPI_AllTemp_Wind" #
# path_ecospace_out = "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v56_3- PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv56_4-RSPI_AllTemp_Wind" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v56_4- PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv56_5-RSPI_AllTemp_Wind" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v16//Sc216 v56_5 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv51_2-PAR_PI_AllPPTemp_Wind" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11 - DEBUG3//Sc216 v51_2 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv51_3-PAR_PI_AllPPTemp_Wind" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11 - DEBUG3//Sc216 v51_3 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv51_4-PAR_PI_AllPPTemp_Wind" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11 - DEBUG3//Sc216 v51_4 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv50_2-PAR_PI_AllPPTemp_Wind" # trying to recreate a good run = note this is starting at 2000 not 2003 as originally
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11 SC50_DEBUG//Sc216 v50_2 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv50_3-PAR_PI_AllPPTemp_Wind" # as above with upped nutrients from 90 to 95 in ecosim
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11 SC50_DEBUG//Sc216 v50_3 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Sc216 v70-PAR_PI_AllPPTemp_Wind" # same as 51_4 (good one), run in debug mode, but now w/ 3 threads
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v17 - DEBUG3//Sc216 v70 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# result of above same as 4_2
# ecospace_code = "Sc216 v71-PAR_PI_AllPPTemp_Wind" # same as 51_4 (good one), run in debug mode, but now w/ 3 threads
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v18 - DEBUG3//Sc216 v71 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"

# ecospace_code = "FULLKEY_Scv51_5-PAR_PI_AllPPTemp_Wind" # same as 51_4 (good one), run in debug mode, but now w/ 3 threads
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2024_Carb_3day_v11 SC51_04 - DEBUG3//Sc216 v51_5 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv80_1-All_Groups_20250501"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_v12_ewe6_7_19295_SC51_04 - DEBUG5//Sc216 v80_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv82_1-All_Groups_20250506" # removed depth response of eup and dec to see if anything changes vs 80_1
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v82_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv83_1-All_Groups_20250506" # now w/ zoop enviro responses
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v83_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv84_1-All_Groups_20250506" # w/ env resp and nutrients up to 96
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v84_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv85_1-All_Groups_20250506" # same as 84 but with DIA temp response adjustment
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v85_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv86_1-All_Groups_20250506" # same as 85 but with DIA PAR response adjustment (_1 adjst 1, _2 adjst 2, _3 adjst again
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v86_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv87_1-All_Groups_20250506" # same as 86 but with a trapezoidal temperature response from diatoms
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v87_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv88_1-All_Groups_20250506" # same as 85 but NO WIND
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v88_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv89_1-All_Groups_20250506" # lost track here - is this where V's up??
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v89_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv90_1-All_Groups_20250506" # same as 89 but now with wind (sc39)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v90_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv91_1-All_Groups_20250506" # same as 90 but sc47 wind
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v91_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv92_1-All_Groups_20250506" # same as 91 but reverted to older temp response for DIA (SC7)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v92_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv93_1-All_Groups_20250506" # same as 92 but reverted to lower V's for groups >13, nutrients from 90 to 96 (previous runs above mislabelleD)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v93_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv94_1-All_Groups_20250506" # same as 93 but with new intermediate wind SC94
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v94_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv95_1-All_Groups_20250506" # same as 92 but PBMax to 100 from 1000
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v95_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc///"
# ecospace_code = "Scv88_2-All_Groups_20250506" # trying to get good QU39 fit to phyto as well as good b timing - based on v88_1 which I assume was low V's, WIND NOW IN, tweaked temp to match 92, nutr 96 instead of 90
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v88_2 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv96_1-All_Groups_20250506" # similar to 92 except low v's and steeper dia temp, sc85 (shares this with v88)
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v96_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv97_1-All_Groups_20250506" # same as 96 with trapezoidal DIA temp response
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v97_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv98_1-All_Groups_20250506" # tweaked all dia responses to keep winter not as low and penalise summer conditions
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v98_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv99_1-All_Groups_20250506" # this is a test of effect of NO wind. Also, adjusted temp response and reverted PAR back
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//Sc216 v99_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv102_1-All_Groups_20250523" # wind back, but just spring-winter wind
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v102_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv103_1-All_Groups_20250523" # same as 102 except tweak to PAR for dia
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v103_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv104_1-All_Groups_20250523" # same as 103 except tweak to PAR and light for dia
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v104_1 - PAR_PI_AllPPTemp_Wind - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv105_1-All_Groups_20250523" # same as 88_2 except with a trapzdl temp rspnse for dia.
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v105_2 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv106_1-All_Groups_20250523" # same as 105 except added wind 'floor' of 0.2 env resp for dia
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v106_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv107_1-All_Groups_20250523" # same as 106 but increased steepness from SC47 env resp for DIA PAR and added higher floor than SC107 of 0.25
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v107_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv108_1-All_Groups_20250523" # same as 107 but increased steepness from SC47 env resp for DIA PAR SC107 of 0.25; LP 18 RP 30, w/ 0.25 floor
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v108_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv109_1-All_Groups_20250523" # same as 108 but now with shallower PAR curve
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v109_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv110_1-All_Groups_20250523" # same as 109 but now with steep dia temp trapz resp w/ 0.2 floor
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v110_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv111_1-All_Groups_20250523" # same as 110 but essentially fixed glitch in declining part of dia temp trap resp
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v111_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv112_1-All_Groups_20250523" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v112_1 - PAR_PI_AllPPTemp_Wind//asc//"
# ecospace_code = "Scv113_1-All_Groups_20250523" #
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v113_1 - PAR_PI_AllPPTemp_Wind//asc//"
ecospace_code = "Scv114_1-All_Groups_20250523" #
path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v12 - DEBUG//v114_1 - PAR_PI_AllPPTemp_Wind//asc//"



#
yr_strt = 1978
yr_end = 2018
mo_strt = 1
da_strt = 2
mo_end = 12
da_end = 30

path_out = "..//..//data//ecospace_out//"
nemo_ewe_csv = "..//..//data//basemap//Ecospace_grid_20210208_rowscols.csv"

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
       "PP3-PIC": path_ecospace_out + "EcospaceMapBiomass-PP3-PIC-{}.asc", #time step format eg: 00620,

       "BA1-BAC": path_ecospace_out + "EcospaceMapBiomass-BA1-BAC-{}.asc",
      }

rows = 151
cols = 93
skiprows = 6 # header

# corresponding lats and lons to centres of grid cells
# read the data into an xarray object which is versatile


out_filename = ecospace_code + "_" + str(yr_strt) + "-" + str(yr_end) + ".nc"
out_filename = os.path.join(path_out, out_filename)


original_ds, are_identical = asc_to_nc_3day(v_f, out_filename, nemo_ewe_csv,
                                            rows, cols, skiprows,
                                            yr_strt, yr_end, mo_strt, da_strt,
                                            mo_end, da_end)
print("Crosscheck ASC in = NC data out:", are_identical)
print("done.")