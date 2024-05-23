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

# Paths
path_chla = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Depth Int Chlorophyl Cit Sci Pawlowicz 2015-2019//ORIGINAL//"
file_chla = "CitSci_dIchl_20220117.csv"

# Initialize lists to hold data
datetimes = []
chlorophyll_values = []

# Read the CSV file
with open(path_chla + file_chla, mode='r', encoding='utf-8') as file:
    reader = csv.reader(file)
    # Skip the first three header lines
    for _ in range(3):
        next(reader)
    # Skip the units line
    next(reader)

    # Process the remaining rows
    for row in reader:
        try:
            # Parse the datetime string to a datetime object
            date = datetime.strptime(row[0], '%d/%m/%Y %H:%M:%S')
            # Append the datetime object to the list
            datetimes.append(date)
            # Append the chlorophyll value to the list (convert to float)
            chlorophyll_values.append(float(row[1].strip()))
        except ValueError:
            # Skip rows with invalid data
            continue

# Convert lists to numpy arrays for easier processing
datetimes = np.array(datetimes)
chlorophyll_values = np.array(chlorophyll_values)

# Extract year, month, and day
years = np.array([date.year for date in datetimes])
months = np.array([date.month for date in datetimes])
days = np.array([date.day for date in datetimes])

# Filter data for years 2015-2018
valid_years = (years >= 2015) & (years <= 2018)
years = years[valid_years]
months = months[valid_years]
days = days[valid_years]

# Group by year, month, and day and count occurrences
unique_dates, counts = np.unique([f'{year}-{month:02d}-{day:02d}' for year, month, day in zip(years, months, days)],
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

# Save as EPS
plt.savefig("..//figs//count_chla_pawlowicz_2015-2018.eps", format='eps', bbox_inches='tight')

# Save as PNG
plt.savefig("..//figs//count_chla_pawlowicz_2015-2018.png", format='png', bbox_inches='tight')

# Save as PDF
plt.savefig("..//figs//count_chla_pawlowicz_2015-2018.pdf", format='pdf', bbox_inches='tight')


plt.show()