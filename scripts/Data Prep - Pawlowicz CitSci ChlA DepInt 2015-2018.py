# Created May 22 2024
# By G Oldford
# Purpose: open the Chl A SoG survey data from R. Pawlowicz
#          output summary stats and ChlA by subregion in SoG
#          for use evaluating the ECOSPACE MODEL (and ECOSIM model)
# To Do:
#       map of sampling locations by year
#       histogram of sampling counts by year and month
#       summary statistics table of sampling by year, month, subregion
#       output summary time series by year,month,day and subregion
# Notes:
#      see archive for some sample ipynb

import pandas as pd
from matplotlib import pyplot as plt, dates
import matplotlib.gridspec as grid_spec
from datetime import datetime, timedelta
from calendar import month_name

# #laptop path
# path_chla = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Depth Int Chlorophyl Cit Sci Pawlowicz 2015-2019//ORIGINAL//"
# #desktop path
# #path_chla = "D://Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Depth Int Chlorophyl Cit Sci Pawlowicz 2015-2019//ORIGINAL//"
# file_chla = "CitSci_dIchl_20220117.csv"
#
# skipheaderlines = 3
# chla_df = pd.read_csv(path_chla + file_chla, skiprows=skipheaderlines, parse_dates=True)
#
# # skip units line
# chla_df = chla_df[1:]
# chla_df["datetime_pd"] = pd.to_datetime(chla_df["datetime"])
# # https://www.dataindependent.com/pandas/pandas-to-datetime/
#
# chla_df = chla_df[chla_df["chl"] != "     NaN"].reset_index()
# chla_df["chl"] = pd.to_numeric(chla_df["chl"])
# #chla_df = chla_df[["chl", "datetime_pd"]]
# chla_df["year_"] = chla_df.datetime_pd.dt.year
# chla_df["month_"] = chla_df.datetime_pd.dt.month
# chla_df["day_"] = chla_df.datetime_pd.dt.day
#
import csv
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
from helpers import find_nearest_point_from1D
import os
import pytz

# Paths
path_chla = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Depth Int Chlorophyl Cit Sci Pawlowicz 2015-2019//ORIGINAL//"
file_chla = "CitSci_dIchl_20220117.csv"
path_modpts = "C://Users//Greig//Documents//github//Ecosystem-Model-Data-Framework//data//basemap//"
file_modpts = "Ecospace_grid_20210208_rowscols.csv"
pathout = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Depth Int Chlorophyl Cit Sci Pawlowicz 2015-2019//MODIFIED//"
out_chla = "CitSci_dIchl_20220117_PST_wEcospaceModPts.csv"

# Initialize lists to hold data
chl_datetimes = []; chl_values = []
chl_lats = []; chl_lons = []

# Read the chl csv file
with open(path_chla + file_chla, mode='r', encoding='utf-8') as file:
    reader = csv.reader(file)
    # Skip the first three header lines
    for _ in range(4):
        if reader.line_num == 3:
            header = next(reader)
            print("Chl-a column names:", header)
        else:
            next(reader)
    # Skip the units line
    next(reader)

    # Process the remaining rows
    for row in reader:
        try:
            # Parse the datetime string to a datetime object
            date = datetime.strptime(row[0], '%d/%m/%Y %H:%M:%S')
            # Append the datetime object to the list
            chl_datetimes.append(date)
            # Append the chlorophyll value to the list (convert to float)
            chl_lats.append(float(row[1].strip()))
            chl_lons.append(float(row[2].strip()))
            chl_values.append(float(row[3].strip()))
        except ValueError:
            # Skip rows with invalid data
            continue

# Convert lists to numpy arrays for easier processing

# Convert UTC to PST
utc_zone = pytz.utc
pst_zone = pytz.timezone('US/Pacific')
chl_datetimes_pst = np.array([utc_zone.localize(date).astimezone(pst_zone) for date in chl_datetimes])

#chl_datetimes = np.array(chl_datetimes)
chl_lats = np.array(chl_lats)
chl_lons = np.array(chl_lons)
chl_values = np.array(chl_values)

# Extract year, month, and day
chl_years = np.array([date.year for date in chl_datetimes_pst])
chl_months = np.array([date.month for date in chl_datetimes_pst])
chl_days = np.array([date.day for date in chl_datetimes_pst])

# Filter data for years 2015-2018
valid_years = (chl_years >= 2015) & (chl_years <= 2018)
chl_years = chl_years[valid_years]
chl_months = chl_months[valid_years]
chl_days = chl_days[valid_years]
chl_values = chl_values[valid_years]


####################### FIND CLOSEST MODEL POINT #########################
# Initialize lists to hold model point data
model_lats = []; model_lons = []
model_deps = []; model_rows = []; model_cols = []

# Read the map file
with open(path_modpts + file_modpts, mode='r', encoding='utf-8') as file:
    reader = csv.reader(file)

    # Read the header line
    header = next(reader)
    print("Model points column names:", header)

    # Process the remaining rows
    for row in reader:
        try:
            # Append the latitude and longitude to the lists (convert to float)
            model_lats.append(float(row[1].strip()))
            model_lons.append(float(row[2].strip()))
            model_deps.append(float(row[3].strip()))
            model_cols.append(float(row[4].strip()))
            model_rows.append(float(row[5].strip()))
        except ValueError:
            # Skip rows with invalid data
            continue

# Convert lists to numpy arrays for easier processing
model_lats = np.array(model_lats)
model_lons = np.array(model_lons)
model_deps = np.array(model_deps)
model_rows = np.array(model_rows)
model_cols = np.array(model_cols)

print(model_lons.shape, model_lats.shape, model_deps.shape, model_cols.shape, model_rows.shape)

#0.01 deg lat = ~1.1 km
dx = 0.01
mask = model_deps != 0

mod_lats_nearest = []; mod_lons_nearest = [];
mod_deps_nearest = []; mod_rows_nearest = []; mod_cols_nearest = []
for lat, lon in zip(chl_lats, chl_lons):
    p = find_nearest_point_from1D(lon, lat, model_lons, model_lats, mask, dx)
    nearest_idx = p[0]
    if not nearest_idx is np.nan:
        mod_lons_nearest.append(model_lons[nearest_idx])
        mod_lats_nearest.append(model_lats[nearest_idx])
        mod_deps_nearest.append(model_deps[nearest_idx])
        mod_rows_nearest.append(model_rows[nearest_idx])
        mod_cols_nearest.append(model_cols[nearest_idx])
    else:
        mod_lons_nearest.append(-999)
        mod_lats_nearest.append(-999)
        mod_deps_nearest.append(-999)
        mod_rows_nearest.append(-999)
        mod_cols_nearest.append(-999)

mod_lons_nearest = np.array(mod_lons_nearest)
mod_lats_nearest = np.array(mod_lats_nearest)
mod_deps_nearest = np.array(mod_deps_nearest)
mod_rows_nearest = np.array(mod_rows_nearest)
mod_cols_nearest = np.array(mod_cols_nearest)
print(mod_lons_nearest.shape, mod_lats_nearest.shape, mod_deps_nearest.shape, mod_rows_nearest.shape ,mod_cols_nearest.shape)
print(chl_lons.shape, chl_lats.shape, chl_values.shape)


####################### EXPORT TO CSV ########################
output_file = os.path.join(pathout, out_chla)
with open(output_file, mode='w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)
    # Write the header
    writer.writerow(["year_pst", "month_pst", "day_pst", "datetimepst", "lat", "lon", "chla_val",
                     "nearest_mod_lon", "nearest_model_lat",
                     "nearest_mod_dep", "nearest_mod_row", "nearest_mod_col"])

    # Write the data rows
    for i in range(len(chl_years)):
        writer.writerow([chl_years[i], chl_months[i], chl_days[i], chl_datetimes_pst[i], chl_lats[i], chl_lons[i], chl_values[i],
                         mod_lons_nearest[i], mod_lats_nearest[i], mod_deps_nearest[i], mod_rows_nearest[i], mod_cols_nearest[i]])

print("wrote file ", file)

count_negative_999 = np.count_nonzero(mod_lons_nearest == -999)
print("number of chla samples outside domain")
print(count_negative_999)
print("number inside")
print(len(mod_lons_nearest) - count_negative_999)

# separate script
# anomaly method - log(mean) chl? clima and anom from from Apr - Sep; anom based on Pawlowicz (used monthly);
# subregions? - N and C for now?
# review McEwan, Mitra, Carls'

####################### VISUAL OF COUNTS BY MONTH ########################
# Group by year, month, and day and count occurrences
unique_dates, counts = np.unique([f'{year}-{month:02d}-{day:02d}' for year, month, day in zip(chl_years, chl_months, chl_days)],
                                 return_counts=True)
unique_dates = [datetime.strptime(date, '%Y-%m-%d') for date in unique_dates]

# Sort by date
sorted_indices = np.argsort(unique_dates)
unique_dates = np.array(unique_dates)[sorted_indices]
counts = counts[sorted_indices]

# Plot the bar chart
plt.figure(figsize=(7.5, 4))
plt.bar(unique_dates, counts, width=1)

# Set x-axis major locator to show every month
# plt.gca().xaxis.set_major_locator(plt.MultipleLocator(30))
#plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
# plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=3))
# plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
plt.gca().xaxis.set_major_locator(mdates.YearLocator(1, month=7, day=3))
plt.gca().xaxis.set_minor_locator(mdates.MonthLocator(bymonthday=1, interval=3))
plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
plt.gca().xaxis.set_minor_formatter(mdates.DateFormatter('%b'))

#plt.xlabel('Date (Year-Month-Day)')
plt.ylabel('Count of Measurements')
plt.title('Count of Measurements by Year, Month, and Day (2015-2018)')
plt.xticks(rotation=0, minor=False, horizontalalignment='center', y=-0.05)
plt.xticks(rotation=0, minor=True, fontsize='small')

#plt.grid(True)

# Save as EPS, PNG, PDF
plt.savefig("..//figs//count_chla_pawlowicz_2015-2018.eps", format='eps', bbox_inches='tight')
plt.savefig("..//figs//count_chla_pawlowicz_2015-2018.png", format='png', bbox_inches='tight')
plt.savefig("..//figs//count_chla_pawlowicz_2015-2018.pdf", format='pdf', bbox_inches='tight')

plt.show()

####################### VISUAL OF AVG CHLA BY DAY ########################
# Combine into a DataFrame
print(chl_years.shape, chl_months.shape, chl_days.shape, chl_values.shape)
data = {
    'year': chl_years,
    'month': chl_months,
    'day': chl_days,
    'chl_value': chl_values,
    'chl_obs_lat': chl_values
}
df = pd.DataFrame(data)

# Create a column for day of the year (1 to 365/366)
df['date'] = pd.to_datetime(df[['year', 'month', 'day']])
df['day_of_year'] = df['date'].dt.dayofyear

# Group by year and day of the year to calculate the average chl values
daily_avg = df.groupby(['year', 'day_of_year'])['chl_value'].mean().reset_index()

# Pivot the DataFrame to have days as rows and years as columns for easier plotting
daily_avg_pivot = daily_avg.pivot(index='day_of_year', columns='year', values='chl_value')

# Plotting
plt.figure(figsize=(12, 6))
for year in daily_avg['year'].unique():
    yearly_data = daily_avg[daily_avg['year'] == year]
    plt.scatter(yearly_data['day_of_year'], yearly_data['chl_value'], label=str(year), s=10)  # s is the size of the points

plt.xlabel('Day of Year')
plt.ylabel('Average Chlorophyll Value')
plt.title('Average Chlorophyll Values by Day for Each Year')
plt.legend(title='Year')
plt.grid(True)
plt.show()

print(len(daily_avg))
print(daily_avg)
print(df)