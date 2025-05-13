import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Load the NetCDF file
file_path = "C:/Users/Greig/Downloads/hotssea_1980_08_temperature_mean.nc"  # Replace with your file path
ds = xr.open_dataset(file_path, decode_times=False)

# Select the variable to plot
var_name = "votemper_bottom"  # Replace with the actual variable name if different
var_name = "votemper_avg150mtoBot"
var_name = "votemper_surface"
data = ds[var_name]

# Get the latitude and longitude
lats = ds["nav_lat"]
lons = ds["nav_lon"]

# Select the first time step (if time is a dimension)
if "time_counter" in data.dims:
    data = data.isel(time_counter=8)

# Determine extent correctly
extent = [lons.min().values, lons.max().values, lats.min().values, lats.max().values]

# Plot the data
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(extent, crs=ccrs.PlateCarree())  # Ensure extent is interpreted in the correct CRS
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, linestyle=':')

# Contour plot
plt.pcolormesh(lons, lats, data, transform=ccrs.PlateCarree(), cmap='viridis')

plt.colorbar(label=var_name)
plt.title(f"{var_name} Visualization")

plt.show()

import cftime
import numpy as np
import pandas as pd

# Extract time variable
time_var = ds["time_counter"]

# Get time attributes
time_units = time_var.attrs["units"]  # e.g., "seconds since 1900-01-01 00:00:00"
calendar = time_var.attrs.get("calendar", "gregorian")

# Identify and mask invalid values (overflowing large numbers)
fill_value = 9.969210e+36
valid_times = np.where(time_var.values > 1e+10, np.nan, time_var.values)  # Mask very large values

# Convert only non-NaN values
valid_time_values = valid_times[~np.isnan(valid_times)]

# Convert using cftime
converted_dates = cftime.num2date(valid_time_values, time_units, calendar=calendar)

# Convert to pandas datetime (handling NaNs)
time_values = pd.to_datetime(converted_dates)

# Extract Year, Month, and Day
years = time_values.year
months = time_values.month
days = time_values.day

print(years, months, days)