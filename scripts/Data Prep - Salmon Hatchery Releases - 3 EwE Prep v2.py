
##3) Prep Release Data for EwE
#Sep 2020 By: G Oldford (last mod: Jan 2025)

#Data In:
#Actual_releases_COORDS.csv - has coordinates and other metadata - includes all species
# note differences between 'STOCK_...' and 'REL_...' fields are 'stock' and 'release' (stock
# location may not be release location)

#Original Data (see previous notebook):
#EPAD data from Carl (DFO / SEP) - has no coordinates - 'actual_releases.csv'
#RMIS release location data (edits by SOGDC) - has coordinates - 'rmis_locations_2019.csv'
#PSF release and release location data - has coordinates - 'PSEReleasesAndLocations2019.csv'
#Coordinate table by Greig - has a few problem coordinates georeferenced with best-guess

#Purpose:
#Export a series of CSV (Biomass, Row, Col), one for each month, a set of series of files for
# each functional group (each species and life stage)
#Step 1: extract start and end month of release (take average)
#Step 2: Aggregate releases by month, species, and "group" (as defined in PhD work)
#Step 3: set a row / col corresponding to the custom, rotated EwE map

#Notes:
# this was ported from old notebooks in Jan 2025
# much of notebook was lost in this repo - had to retrieve from the failed CEM-SOG repo Jan 2025
#EPAD data from Carl Walters and RMIS locations data from SOGDC
#rmis_smolt_releases dataset from 'rmis_releases.csv' from http://sogdatacentre.ca/search-data/spatial-data/ from all_layers->rmis->rmis_smolt_releases
#prioritized effort for coordinate matching on just coho and Chinook


### 1. READ FILES, INSPECT ###
import pandas as pd
import numpy as np
import datetime, calendar
import math

# locations table from the SSMSP SOGDC (may have more lats / lons added than source at RMIS)
path = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\All Species Hatchery Releases\\EPADHatcherReleasesGST\\MODIFIED\\"
releases_df = pd.read_csv(path + "actual_releases_COORDS.csv")

# table above was edited in ArcMap to find closest model row / cols
# read in result:
fromarcmap_df = pd.read_csv(path + "allspecies_coordseditedinarcmap.csv")
fromarcmap_df['UniqueID'] = fromarcmap_df['Field1']

print("Quick check - chinook in 1988 should be around 33 mil, with no other years higher")
print(releases_df.loc[(releases_df['SPECIES_NAME']=='Chinook')].groupby(['BROOD_YEAR'])['TotalRelease'].sum())
print(releases_df.columns)
print(releases_df[['UniqueID','SPECIES_NAME','REL_CU_NAME','BIOMASS_MT','START_YR_REL',
             'START_MO_REL','START_DAY_REL','START_DATE','END_YR_REL','END_MO_REL',
             'END_DAY_REL','RELEASE_YEAR']])


#### 2. Fix issues with dates ####
# fix NaN's in days (put day = 15 in this case)
releases_df['START_DAY_REL'] = releases_df['START_DAY_REL'].fillna(15)
releases_df['END_DAY_REL'] = releases_df['END_DAY_REL'].fillna(15)
releases_df['START_MO_REL'] = releases_df['START_MO_REL'].fillna(1)
releases_df['END_MO_REL'] = releases_df['END_MO_REL'].fillna(1)

# convert to integers
releases_df['START_DAY_REL'] = releases_df['START_DAY_REL'].astype(int)
releases_df['END_DAY_REL'] = releases_df['END_DAY_REL'].astype(int)
releases_df['START_MO_REL'] = releases_df['START_MO_REL'].astype(int)
releases_df['END_MO_REL'] = releases_df['END_MO_REL'].astype(int)

#print(releases_df[['START_MO_REL','START_DAY_REL','END_MO_REL','END_DAY_REL']])
# helpful datetime functions
def makedate(year_col, month_col, day_col):
    return (datetime.date(year_col, month_col, day_col))


def getavgrel_month(year_start, month_start, day_start, year_end, month_end, day_end):
    date_start = makedate(year_start, month_start, day_start)
    date_end = makedate(year_end, month_end, day_end)

    date_avg = date_start + ((date_end - date_start) / 2)
    return (date_avg.month)


def getavgrel_year(year_start, month_start, day_start, year_end, month_end, day_end):
    date_start = makedate(year_start, month_start, day_start)
    date_end = makedate(year_end, month_end, day_end)

    date_avg = date_start + ((date_end - date_start) / 2)
    return (date_avg.year)


def getavgrel_datetime(year_start, month_start, day_start, year_end, month_end, day_end):
    date_start = makedate(year_start, month_start, day_start)
    date_end = makedate(year_end, month_end, day_end)

    date_avg = date_start + ((date_end - date_start) / 2)
    return (date_avg)

############################################
#### 3. Calculate avg month of release  ####
# using average of start and end day
# find avg date of releases
releases_df['release_avg_month'] = releases_df.apply(lambda x: getavgrel_month(x.START_YR_REL,
                                                                               x.START_MO_REL,
                                                                               x.START_DAY_REL,
                                                                               x.END_YR_REL,
                                                                               x.END_MO_REL,
                                                                               x.END_DAY_REL), axis=1)
releases_df['release_avg_year'] = releases_df.apply(lambda x: getavgrel_year(x.START_YR_REL,
                                                                             x.START_MO_REL,
                                                                             x.START_DAY_REL,
                                                                             x.END_YR_REL,
                                                                             x.END_MO_REL,
                                                                             x.END_DAY_REL), axis=1)

#REL_CU_NAME, SPECIES_NAME, release_avg_month, release_avg_year
print(releases_df['release_avg_month'].unique())

print(releases_df.loc[(releases_df['SPECIES_NAME']=='Chinook')].groupby(['REL_CU_NAME',
                                                                         'release_avg_year',
                                                                         'release_avg_month'])['BIOMASS_MT'].sum().reset_index())

# create a datetime column (python format) for average release date
releases_df['release_avg_date'] = releases_df.apply(lambda x: getavgrel_datetime(x.START_YR_REL,
                                                                                 x.START_MO_REL,
                                                                                 x.START_DAY_REL,
                                                                                 x.END_YR_REL,
                                                                                 x.END_MO_REL,
                                                                                 x.END_DAY_REL), axis=1)

print(releases_df.head())
outpath = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\All Species Hatchery Releases\\EPADHatcherReleasesGST\\MODIFIED\\"
releases_df.to_csv(outpath + 'temp_2025.csv')



#########################################
#### 4. Assign model map row / col  #####
#########################################
# - have most lats / lons but many are off the marine portion of the map -
# requires either a script or manual adjustments
# - decided to do it in ArcMap, all manually, including assigning row / cols
# - I began investigating how to do it with Release CU Name (see below) but too many
# NULLS so did the moving of release points and then the assigning of rows and cols manually.

# this is just a hack first pass (see above)
# row / col for each location
#fraser riv (IFR, LFR): 131, 45
#cowichan (COW): 110, 7
#upper georgia strait (UGS): 34, 12
#georgia strait mainland (UGS): 50, 29
#howe sound (?): 104, 55
print(releases_df['REL_CU_NAME'].unique())

# there are issues with 31% of the records - empty / NaN / null values for REL_CU_NAME
releases_df['problemswithCUNAME'] = releases_df['REL_CU_NAME'].isnull()
nullvalues = len(releases_df.loc[(releases_df['problemswithCUNAME']==True)])
notnullvalues = len(releases_df.loc[(releases_df['problemswithCUNAME']==False)])
print(nullvalues / (notnullvalues + nullvalues))



################################################################################
#############################################################
#Summary of matching EWE map row / col to release lat / lon
# - I took all hatchery releases and found a corresponding Row / Col using ArcMap.
# For most I had to manually find the most reasonable nearest map model row / col.
# - next step is to re-export as a CSV time series for EwE
# - need to know the EwE functional group code, though
# utility code - not used (did changes in ArcMap)
#df.loc[df.my_channel > 20000, 'my_channel'] = 0
#releases_df.loc[(releases_df['REL_CU_NAME']=='LOWER FRASER RIVER_FA_0.3'),'EWE_ROW'] = 131
#releases_df.loc[(releases_df['REL_CU_NAME']=='LOWER FRASER RIVER_FA_0.3'),'EWE_COL'] = 45
#releases_df.loc[(releases_df['REL_CU_NAME']=='LOWER FRASER'),'EWE_ROW'] = 131
#releases_df.loc[(releases_df['REL_CU_NAME']=='LOWER FRASER'),'EWE_COL'] = 45
#releases_df.columns

#releases_df['UniqueID'] = releases_df[.index]
releases_df.head()

# join
releases_ewerowscols = pd.merge(releases_df, fromarcmap_df, on=['UniqueID'], how='left')
# these columns should match except I've inserted 15 when start_day_rel = 0
releases_ewerowscols[['START_DAY_REL','START_DAY_']]
print(releases_ewerowscols.columns)

#1/3 of records have no release conservation unit
temp = releases_ewerowscols[['REL_CU_NAME','REL_CU_INDEX']]
#temp.drop_duplicates()
print(temp[temp['REL_CU_NAME'].notnull()])

# how many have no CU for after 1978?
temp = releases_ewerowscols[['REL_CU_NAME','REL_CU_INDEX','START_YR_REL']].loc[releases_ewerowscols['START_YR_REL']>1978]
print(temp[temp['REL_CU_NAME'].notnull()])
# roughly same

# conclusion: the release conservation unit name is empty 1/3 of the time.
print("the release conservation unit name is empty 1/3 of the time.")

#note that 1/3 of records have no release conservation unit
temp = releases_ewerowscols[['START_YR_REL','GAZETTED_NAME','REL_CU_NAME','REL_CU_INDEX','ROW_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]

# 'STOCK_PROD_AREA_CODE' will have to be used for records without conservation unit and stock_cu index
temp[(temp['STOCK_CU_INDEX'].isnull())&(temp['REL_CU_INDEX'].isnull())&
     (temp['SPECIES_NAME']!='Steelhead')&(temp['SPECIES_NAME']!='Cutthroat')]

# temp[(temp['STOCK_CU_INDEX'].isnull())&(temp['REL_CU_INDEX'].isnull())&
#      (temp['SPECIES_NAME']!='Steelhead')&(temp['SPECIES_NAME']!='Cutthroat')]


################ 5) Match releases to the EwE Functional Group  ################
# EwE functional groups ID and text name as stored in data model
# (actual Ecospace or Ecopath ID will differ)
# Groups:
# Chinook-H-IFR-2 - Interior Fraser Spring 4_2 / 1.2 (CK-16, CK-17), Fraser Spring 5_2 / 1.3 (CK-10, CK-12, CK-14, CK-18), Fraser Summer 5_2  / 1.3 (CK-09, CK-11, CK-19), Fraser Summer 4_1    / 0.3 (CK-07, CK-13, CK-15); CWT Indicator Stocks: Nicola (NIC), Shuswap (SHU), Middle Shuswap (MSH)
# Chinook-H-LFR-2 - Management Unit and Conservation Units: Fraser Fall 4_1 / 0.3 (CK-2, CK-3, CK-4, CK-5, CK-6, CK-20,CK-9006, CK-9007,CK-9008); CWT Indicator Stocks: Harrison (HAR), Chilliwack (CHI)
# Chinook-H-COW-2 - Management Unit and Conservation Units: Lower Georgia Strait Nanaimo to Cowichan (CK-21, CK-22, CK-25); CWT Indicator Stocks: Cowichan (COW)
# Chinook-H-UGS-2 - Management Unit and Conservation Units: Upper Georgia Strait (CK-28, CK-29, CK-27, CK-83); CWT Indicator Stocks: Quinsam (QUI), Phillips (PHI), Puntledge (PPS), Big Qualicum (BQR)
# Coho-H-IFR-2 - Management Unit and Conservation Units: Interior Fraser (CO-4, CO-5, CO-6, CO-7, CO-8, CO-9, CO-48). CWT Indicator stocks (uncertain - see CW / JK spreadsheet; coldwater, deadman, spius, dunn, louis, lemieux, eagle)
# Coho-H-LFR-2 - Management Unit and Conservation Units: Lower Fraser & Boundary Bay (CO-1, CO-2, CO-3, CO-10, CO-47). CWT Indicator stocks: Inch Cr, Louis (Dunn Creek H), Chilliwack
# Coho-H-UGS-2 - (lat < 49, CO-13, CO-11) Management Unit and Conservation Units: East Coast Vancouver Island + Georgia Strait (CO-11,CO-13); CWT indicator stocks: Quinsam, Big Qualicum, Black, Puntledge, Goldstream
# Coho-H-COW-2 - (lat >= 49, CO-13, CO-11) Management Unit and Conservation Units: East Coast Vancouver Island + Georgia Strait (CO-11, CO-13); CWT indicator stocks: none

# Rearing type codes (REARING_TYPE_CODE)
# H - hatchery, seapen, lakepen, rearing channel.
# W - wild unfed
# F - wild fed (held short term in pens in river prior to release)
# U - unknown

# stock type codes (STOCK_TYPE_CODE)
# H - Hatchery
# W - Wild
# M - Mixed (hatchery and wild)
# U - Unknown

#Souce: Serbic, G. 1992. The Finclip Recovery Database & Reporting System

releases_ewerowscols[(releases_ewerowscols["REARING_TYPE_CODE"]=='W')&(releases_ewerowscols["SPECIES_NAME"]=='Chinook')]
releases_ewerowscols["REARING_TYPE_CODE"].unique()

################ 5a) Round 1 matching ##################
# Round 1 matching: use REL_CU_INDEX (release conservation unit ID) to match releases to EWE model group
# Chinook
releases_ewerowscols["EWE_GROUP_CODE"] = "x"
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CK-16") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-17") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-10") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-12") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-14") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-18") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-09") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-11") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-19") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-7") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-13") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-15")), "Chinook-H-IFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CK-2") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-3")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-4") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-5") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-6") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-9006")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-9007")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-9008")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-20")), "Chinook-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CK-21") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-22") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-25")), "Chinook-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CK-27") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-28") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-29")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CK-83")), "Chinook-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])
# Coho
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CO-4") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-5") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-6") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-7") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-8") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-9")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-48")), "Coho-H-IFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CO-1") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-2") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-3") |
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-10")|
                                                   (releases_ewerowscols["REL_CU_INDEX"] == "CO-47")), "Coho-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CO-13") |
                                                  (releases_ewerowscols["REL_CU_INDEX"] == "CO-11")) &
                                                  (releases_ewerowscols["FINAL_LAT"] >= 49), "Coho-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["REL_CU_INDEX"] == "CO-13") &
                                                  (releases_ewerowscols["FINAL_LAT"] < 49), "Coho-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])

temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME','REL_CU_NAME','REL_CU_INDEX','ROW_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]
print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))


############## 5b) Round 2 matching: use stock_cu_name #################
# Round 1 matching: use STOCK_CU_INDEX (release conservation unit ID) to match releases to EWE model group
# Chinook
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CK-16") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-17") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-10") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-12") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-14") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-18") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-09") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-11") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-19") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-7") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-13") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-15")), "Chinook-H-IFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CK-2") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-3") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-4") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-5") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-6") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-9006") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-9007") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-9008")|
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-20")), "Chinook-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CK-21") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-22")|
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-25")), "Chinook-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CK-27") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-28") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-29")|
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CK-83")), "Chinook-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])
# Coho
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CO-4") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-5") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-6") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-7") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-8") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-9") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-48")), "Coho-H-IFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CO-1") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-2") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-3") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-10") |
                                                   (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-47")), "Coho-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CO-13") |
                                                  (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-11")) &
                                                  (releases_ewerowscols["FINAL_LAT"] >= 49), "Coho-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"] == "x") &
                                                  (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-13") &
                                                  (releases_ewerowscols["FINAL_LAT"] < 49), "Coho-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])

temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME','REL_CU_NAME','REL_CU_INDEX','ROW_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]
print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))


temp2 = temp[(temp['EWE_GROUP_CODE']=='x') &
             (temp['REL_CU_NAME'].notnull()) &
             (temp['REARING_TYPE_CODE'] == 'H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]
print(temp2['REL_CU_INDEX'].unique())

# some of the remaining unmatched COW and UGS coho releases have no 'final_lat', so use the 'ewe_row' instead
# note rows begin index number in north (northern row = 0)
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CO-13")|
                                                  (releases_ewerowscols["REL_CU_INDEX"] == "CO-11")) &
                                                  (releases_ewerowscols["ROW_EWE"] < 100), "Coho-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  ((releases_ewerowscols["REL_CU_INDEX"] == "CO-13") |
                                                  (releases_ewerowscols["REL_CU_INDEX"] == "CO-11") ) &
                                                  (releases_ewerowscols["ROW_EWE"] >= 100), "Coho-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CO-13")|
                                                  (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-11")) &
                                                  (releases_ewerowscols["ROW_EWE"] < 100), "Coho-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  ((releases_ewerowscols["STOCK_CU_INDEX"] == "CO-13") |
                                                  (releases_ewerowscols["STOCK_CU_INDEX"] == "CO-11") ) &
                                                  (releases_ewerowscols["ROW_EWE"] >= 100), "Coho-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])


temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME','REL_CU_NAME','REL_CU_INDEX','ROW_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]
print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))



####################### 5c) Round 3 matching: use EWE row / col and species #########################
# Rows and cols corresponding to the model were matched earlier to lats and lons of releases
# the small number of unmatched releases can be matched to EWE functional groups using ewe release rows / cols,
# and species name

# if row <= 90 then release likely can be associated with 'upper strait of georgia' groups

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["ROW_EWE"] <= 90), "Coho-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["ROW_EWE"] <= 90), "Chinook-H-UGS-2", releases_ewerowscols["EWE_GROUP_CODE"])

temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME','REL_CU_NAME','REL_CU_INDEX','ROW_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]
print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))

####################### 5d) final matching - ad hoc ####################
# what remains is small (1-2% of records) and unusual releases
temp2 = temp[(temp['EWE_GROUP_CODE']=='x') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]

# Some of the records are Chinook Nitinat river stock released in Esquimalt harbour
#print(temp2[temp2['STOCK_CU_INDEX']=="CK-31"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["STOCK_CU_INDEX"] =="CK-31"), "Chinook-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])

temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME','REL_CU_NAME','REL_CU_INDEX','ROW_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]
print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))

# use the 'stock_prod_area_code' for the remainder
# - did not use this prior because the stock production area is different from the release area
#   and long distance transplants are common
temp2 = temp[(temp['EWE_GROUP_CODE']=='x') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]

# chinook
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["STOCK_PROD_AREA_CODE"] =="LWFR"), "Chinook-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  ((releases_ewerowscols["STOCK_PROD_AREA_CODE"] =="UPFR") |
                                                  (releases_ewerowscols["STOCK_PROD_AREA_CODE"] =="TOMF")), "Chinook-H-IFR-2", releases_ewerowscols["EWE_GROUP_CODE"])

# coho
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["STOCK_PROD_AREA_CODE"] =="LWFR"), "Coho-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  ((releases_ewerowscols["STOCK_PROD_AREA_CODE"] =="UPFR") |
                                                  (releases_ewerowscols["STOCK_PROD_AREA_CODE"] =="TOMF")), "Coho-H-IFR-2", releases_ewerowscols["EWE_GROUP_CODE"])

temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME',
                             'REL_CU_NAME','REL_CU_INDEX','ROW_EWE','COL_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]
print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))

# the remaining unmatched appear distributed in the south.
# assign the remaining ones to groups purely based on release row / col
# temp2 = temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
#      ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]
# temp2['COL_EWE'].unique()

# chinook
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["ROW_EWE"] > 90) &
                                                  (releases_ewerowscols["COL_EWE"] <= 20), "Chinook-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Chinook") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["ROW_EWE"] > 90) &
                                                  (releases_ewerowscols["COL_EWE"] > 20), "Chinook-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])

# coho
releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["ROW_EWE"] > 90) &
                                                  (releases_ewerowscols["COL_EWE"] <= 20), "Coho-H-COW-2", releases_ewerowscols["EWE_GROUP_CODE"])

releases_ewerowscols["EWE_GROUP_CODE"] = np.where((releases_ewerowscols["SPECIES_NAME"] == "Coho") &
                                                  (releases_ewerowscols["REARING_TYPE_CODE"] == "H") &
                                                  (releases_ewerowscols["EWE_GROUP_CODE"]=="x") &
                                                  (releases_ewerowscols["ROW_EWE"] > 90) &
                                                  (releases_ewerowscols["COL_EWE"] > 20), "Coho-H-LFR-2", releases_ewerowscols["EWE_GROUP_CODE"])


temp = releases_ewerowscols[['EWE_GROUP_CODE','REARING_TYPE_CODE','START_YR_REL','GAZETTED_NAME',
                             'REL_CU_NAME','REL_CU_INDEX','ROW_EWE','COL_EWE',
                             'RELEASE_SITE_NAME','STOCK_NAME','STOCK_CU_INDEX','SPECIES_NAME','STOCK_PROD_AREA_CODE']]

print("the number of release records:")
print(print(len(temp[((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook')) &
            (temp['REARING_TYPE_CODE']=='H')])))

print("the number of unmatched release records:")
print(len(temp[(temp['EWE_GROUP_CODE']=='x') & (temp['REARING_TYPE_CODE']=='H') &
     ((temp['SPECIES_NAME']=='Coho')|(temp['SPECIES_NAME']=='Chinook'))]))

print(releases_ewerowscols.columns)



temp3 = releases_ewerowscols[['EWE_GROUP_CODE','SPECIES_NAME',
                              'AVE_WEIGHT',
                              'TotalRelease','BIOMASS_MT',
                              'release_avg_date',
                              'FINAL_LAT','FINAL_LON',
                             'ROW_EWE','COL_EWE']]
# temp3[temp3['EWE_GROUP_CODE']=="Chinook-H-LFR-2"]
# temp3[temp3['EWE_GROUP_CODE']=="Chinook-H-LFR-2"].to_csv(outpath + "temp3.csv")
temp3.to_csv(outpath + "temp3.csv")
print(releases_ewerowscols.columns)


################### WRITE TO FILE ######################
outpath = "C:\\Users\\Greig\\Sync\\6. SSMSP Model\\Model Greig\\Data\\1. Salmon\\All Species Hatchery Releases\\EPADHatcherReleasesGST\\MODIFIED\\"

temp3['TOTRELEASE_NO'] = temp3[['TotalRelease']]
out_df = temp3.drop(['TotalRelease'], axis=1)
# ID from source_meta table in data framework
out_df['SOURCE_ID'] = 2
out_df.to_csv(outpath + 'HatcheryRel_TS_ForNextstep_2025.csv')