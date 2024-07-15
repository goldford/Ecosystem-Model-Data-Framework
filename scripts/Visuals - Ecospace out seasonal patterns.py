import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as dates
import calendar

#ecospace_out_code = "Scv2-MultMix"
#ecospace_out_code = "Scv3-MultxPAR"
#ecospace_out_code = "Scv4-MultxPARlimZ"
#ecospace_out_code = "Scv5-PARMixingNut90"
ecospace_out_code = "Scv7-PARMixingNut90Temp"
year_start = 2003
year_end = 2018
# path_in = "..//data//ecospace_out//"
path_in = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
nc_filename = ecospace_out_code + "_" + str(year_start) + "-" + str(year_end) + ".nc"
#nc_filename = ecospace_out_code + "_" + str(year_start) + "-" + str(year_end) + "_2" + ".nc"
# Load the dataset
ds = xr.open_dataset(os.path.join(path_in, nc_filename))

# years to visualize
year_start = 2003
year_end = 2018

# Variable to plot (replace with the relevant variable name from NC file)
#variable_to_plot = "PP1-DIA"
# Groups to plot
groups = ["PP1-DIA", "PP2-NAN", "PP3-PIC", "PZ1-CIL", "PZ2-DIN"]
scaling_factors = {"PZ1-CIL": 0.125, "PZ2-DIN": 0.125} # proportion autotrophic

# Mapping of variable names to readable labels
readable_labels = {
    "PP1-DIA": "Diatoms",
    "PP2-NAN": "Nanophytoplankton",
    "PP3-PIC": "Picophytoplankton",
    "PZ1-CIL": "Ciliates",
    "PZ2-DIN": "Dinoflagellates",
    "PZ3-HNF": "Heterotrophic Nanoplankton"
}

# Initialize a dictionary to store the mean values for each group
mean_values = {}

fig, ax = plt.subplots(figsize=(6, 7))

# Plotting properties
bbox = {'boxstyle': 'round', 'facecolor': 'w', 'alpha': 0.9}
#cmap = plt.get_cmap('tab10')
#palette = [cmap(0), cmap(0.2), 'k', cmap(0.1), cmap(0.3)]

# Define a colormap and generate palette
cmap = plt.get_cmap('viridis')
num_colors = len(groups)
palette = [cmap(i / num_colors) for i in range(num_colors)]

# Loop over each group
for variable_to_plot, color in zip(groups, palette):

    # Extract the data for the specific variable and time range
    # Note: Adjust this extraction based on the actual structure of your NetCDF file
    data_var = ds[variable_to_plot].sel(time=slice(f"{year_start}-01-01", f"{year_end}-12-31"))

    # Average the data over the spatial dimensions (rows and columns)
    data_avg = data_var.mean(dim=['row', 'col'])

    # Apply scaling if the group is in scaling_factors
    if variable_to_plot in scaling_factors:
        data_avg *= scaling_factors[variable_to_plot]

    # Store the mean values for later summation
    mean_values[variable_to_plot] = data_avg.groupby('time.dayofyear').mean(dim='time')

    # Convert the time coordinate to datetime objects
    time_coords = data_avg.time.to_index().to_pydatetime()

    # Plot the data
    # for year, color in zip(range(year_start, year_end + 1), palette):
    #     index = np.logical_and(time_coords >= datetime(year, 1, 1), time_coords < datetime(year, 12, 31))
    #     x = [datetime(2007, t.month, t.day) for t in time_coords[index]]
    #     y = data_avg[index].values.flatten()  # Flatten the array for plotting
    #     ax.plot(x, y, color=color, label=year)

    # Average the data over the spatial dimensions (rows and columns)
    data_avg_annual = data_var.groupby('time.dayofyear').mean(dim='time')
    data_avg_annual = data_avg_annual.mean(dim=['row', 'col'])
    #time_coords_annual = data_avg_annual.dayofyear.to_index().to_pydatetime()
    days_in_year = [datetime(2007, 1, 2) + timedelta(days=i*3) for i in range(120)]
    # Plot the mean line
    ax.plot(days_in_year, data_avg_annual, color=color, linestyle='-', linewidth = 2, label=readable_labels[variable_to_plot])

# Sum the mean values across all groups
total_mean = sum(mean_values.values())

# Plot the total mean line
#ax.plot(days_in_year, total_mean, color='black', linestyle='-', label='Total Mean')

# Set x-axis limits and labels
ax.set_xlim([datetime(2007, 1, 1), datetime(2007, 12, 31)])
ax.set_ylabel('Biomass (C g m^2)')
ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
ax.legend()


fig.savefig('..//figs//' + 'seasonal_PP_' + str(year_start) + "-" + str(year_end) + ".png" )
plt.show()

######## PLOT JUST DIATOMS TO SHOW INTERYEAR VARIABILITY ########
group = "PP1-DIA"
print("plotting group " + group)
fig, ax = plt.subplots(figsize=(6, 7))

# Extract the data for the specific variable and time range
# Note: Adjust this extraction based on the actual structure of your NetCDF file
data_var = ds[group].sel(time=slice(f"{year_start}-01-01", f"{year_end}-12-31"))

# Average the data over the spatial dimensions (rows and columns)
data_avg = data_var.mean(dim=['row', 'col'])

variable_to_plot = group

# # Apply scaling if the group is in scaling_factors
# if variable_to_plot in scaling_factors:
#     data_avg *= scaling_factors[variable_to_plot]

# Convert the time coordinate to datetime objects
time_coords = data_avg.time.to_index().to_pydatetime()

cmap = plt.get_cmap('viridis')
num_colors = year_end - year_start + 1
palette = [cmap(i / num_colors) for i in range(num_colors)]

# Plot the data
print("plotting each year")
for year, color in zip(range(year_start, year_end + 1), palette):
    index = np.logical_and(time_coords >= datetime(year, 1, 1), time_coords < datetime(year, 12, 31))
    x = [datetime(2007, t.month, t.day) for t in time_coords[index]]
    y = data_avg[index].values.flatten()  # Flatten the array for plotting
    ax.plot(x, y, color=color, label=year)

print("done plotting each year")
ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
fig.savefig('..//figs//' + 'DIATOM_B_' + str(year_start) + "-" + str(year_end) + ".png" )

# # Average the data over the spatial dimensions (rows and columns)
# data_avg_annual = data_var.groupby('time.dayofyear').mean(dim='time')
# data_avg_annual = data_avg_annual.mean(dim=['row', 'col'])
# #time_coords_annual = data_avg_annual.dayofyear.to_index().to_pydatetime()
# days_in_year = [datetime(2007, 1, 2) + timedelta(days=i*3) for i in range(120)]
# # Plot the mean line
# ax.plot(days_in_year, data_avg_annual, color=color, linestyle='-', linewidth = 2, label=readable_labels[variable_to_plot])
#
# # Sum the mean values across all groups
# total_mean = sum(mean_values.values())

# Plot the total mean line
#ax.plot(days_in_year, total_mean, color='black', linestyle='-', label='Total Mean')

# Set x-axis limits and labels
ax.set_xlim([datetime(2007, 1, 1), datetime(2007, 12, 31)])
ax.set_ylabel('Biomass (C g m^2)')
ax.xaxis.set_major_locator(dates.MonthLocator())
ax.xaxis.set_major_formatter(dates.DateFormatter('%b'))
# ax.legend()

# fig.savefig('..//figs//' + 'seasonal_PP_' + str(year_start) + "-" + str(year_end) + ".png" )
plt.show()





# #####
# # Extract month information from datetime objects
# months = [dt.month for dt in time_coords]
#
# # Define the number of bins for the heatmap
# num_bins = 38
#
# # Create a 2D histogram of the data
# hist, xedges, yedges = np.histogram2d(months, data_avg.values.flatten(), bins=[12, num_bins])
#
# # Plotting properties
# fig, ax = plt.subplots(figsize=(12, 6))
#
# # Plot the heatmap
# pcm = ax.pcolormesh(xedges, yedges, hist.T, cmap='viridis', shading='auto')
#
# # Add colorbar
# cbar = fig.colorbar(pcm, ax=ax, label='Frequency')
#
# # Set labels and title
# ax.set_xlabel('Month')
# ax.set_ylabel('Value')
# ax.set_title('Annual Pattern')
#
# # Set x-axis ticks and labels
# ax.set_xticks(np.arange(1, 13))
# ax.set_xticklabels(calendar.month_abbr[1:])  # Assuming you imported 'calendar'
#
# # Show plot
# plt.show()



# # For the second group PP2-NAN, you can follow similar steps as above with necessary modifications.
# # Extract data for PP2-NAN
# variable_to_plot_2 = "PP2-NAN"
# data_var_2 = ds[variable_to_plot_2].sel(time=slice(f"{year_start}-01-01", f"{year_end}-12-31"))
# data_avg_2 = data_var_2.mean(dim=['row', 'col'])
#
# # Build figure layout
# fig, ax = plt.subplots(figsize=(15, 6))
#
# # Plot the data for PP2-NAN
# for year, color in zip(range(year_start, year_end + 1), palette):
#     index = np.logical_and(time_coords >= datetime(year, 1, 1), time_coords < datetime(year, 12, 31))
#     x = [datetime(year, t.month, t.day) for t in time_coords[index]]
#     y = data_avg_2[index].values.flatten()  # Flatten the array for plotting
#     ax.plot(x, y, color=color, label=f'{variable_to_plot_2} - {year}')
#
# # Plot the mean line
# mean_line_2 = data_avg_2.mean(dim='time')
# ax.plot(time_coords, mean_line_2, color='black', linestyle='--', label=f'{variable_to_plot_2} - Mean')
#
# # Set x-axis limits and labels
# ax.set_xlim([datetime(year_start, 1, 1), datetime(year_end, 12, 31)])
# ax.set_ylabel('B')
# ax.xaxis.set_major_locator(dates.YearLocator())
# ax.xaxis.set_minor_locator(dates.MonthLocator())
# ax.xaxis.set_major_formatter(dates.DateFormatter('%Y'))
# ax.legend()
#
# # Show plot
# plt.show()