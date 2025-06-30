# ============================================================================================
# Script:   Process Ecosim Biomass output to CSV that has dates and seasons
# By: G Oldford, 2025
# Purpose:
#    Create 3-day averaged time series from Nemcek et al. observations,
#    with custom subdomain groupings
# ============================================================================================

# Import pandas
import pandas as pd
from datetime import datetime, timedelta

# Define the path to the CSV file
ECOPATH_NM = "ECOSPACE_KEYRUN_LTL_2025_Carb_3day_ewe6_7_19295_v13_BigPC"
SCENARIO = "SC123"
OUTPUT_F = "biomass_monthly.csv"
FILE_PATH = f"C:/Users/Greig/Documents/EwE output/{ECOPATH_NM}/ecosim_Ecosim_{SCENARIO}/{OUTPUT_F}"
OUTPUT_EXPORT_PATH = f"C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation/ecosim_{SCENARIO}_onerun_B_dates_seasons.csv"
TIMESTEP_DAYS = 3
HEADER_N = 14
MAX_TIMESTEPS = 120 # per year (will lump >120 with 120)
# the model timestep output is not assumed to be correct in the outputs
# because of using the 3-day time step etc.
YEAR_START = 1978
YEAR_END = 2018

# === ADD SEASON COLUMN ===
def get_season(date):
    m = date.month
    if m in [11, 12, 1, 2]:
        return 'Winter'
    elif m in [3, 4, 5]:
        return 'Spring'
    elif m in [6, 7, 8]:
        return 'Summer'
    else:
        return 'Fall'

# Read the CSV file, skipping the first 14 lines (0-based)
df = pd.read_csv(FILE_PATH, skiprows=HEADER_N)

# Show the first few rows to confirm it loaded correctly
print(df.head())

# === CALCULATE DATES ===
timestep_col = df.columns[0]  # identifies 'timestep\\group'
start_date = datetime(YEAR_START, 1, 1)

# Function to calculate date with year rollover
def calc_date(timestep):
    year_offset = (timestep - 1) // MAX_TIMESTEPS
    timestep_in_year = ((timestep - 1) % MAX_TIMESTEPS) + 1
    year = YEAR_START + year_offset
    year_start_date = datetime(year, 1, 1)
    # Middle day of timestep
    date = year_start_date + timedelta(days=(timestep_in_year - 1) * TIMESTEP_DAYS + TIMESTEP_DAYS // 2)
    return date

df['date'] = df[timestep_col].apply(calc_date)
df['season'] = df['date'].apply(get_season)

# === OUTPUT PREVIEW ===
# df_date = df[['date', 'season']]
# print(df[[timestep_col, 'date']].head())

# === EXPORT TO CSV ===
df.to_csv(OUTPUT_EXPORT_PATH, index=False)
print(f"Exported updated dataframe to {OUTPUT_EXPORT_PATH}")
