# Created by G Oldford
# Aug 7, 2024
# Purpose: process daily RDRS v2.1 wind files (NC) to 3-day fields for testing of forcing Ecospace
#
# Source:
#
# Input:
#   1/ NC files of RDRS daily processed on server
#
# Output:
#   1/ figs of estimated bloom timing from Ecospace compared to observations
#   2/
#
# note:
#   - not running this on the graham server because it seems to be down.

import os
import numpy as np

wind_daily_p = "C://Users//Greig//Downloads//RDRS_Wind"

year_st = 1980
year_en = 1981
file example = "RDRS21_NEMOgrid_wind_daily_2018.nc"