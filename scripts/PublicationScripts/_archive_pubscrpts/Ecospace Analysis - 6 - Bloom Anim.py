# Created Sep 17 2024 by G Oldford
# Purpose: Create quality still images and animation of spring bloom
# inputs:
#  - nc file of Ecospace outputs
# outputs:
#  - png images
#  - videos (mpg?)

import numpy as np
import os
import xarray as xr
from helpers import read_sdomains, make_map, adjust_map, buildSortableString
from matplotlib.path import Path as Path
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import ListedColormap
from datetime import datetime, timedelta
import netCDF4 as nc
import numpy as np
import pandas as pd
import cmocean

path_ecospace_out = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT//"
scenario = 'FULLKEY_SC51_5'
ecospace_code = "FULLKEY_Scv51_5-PAR_PI_AllPPTemp_Wind"
filenm_yr_strt = 1978

log_transform = False

# don't need since depth is in the ds
# with nc.Dataset('..//..//data//bathymetry//bathy_salishsea_1500m_20210706.nc') as b:
#     # print(mesh.variables)
#     navlon = b.variables['nav_lon'][:]
#     navlat = b.variables['nav_lat'][:]
#     bathy = b.variables['Bathymetry'][:]
#     #     bathy=np.sum(mesh.variables['bathymetry'][0,:])#*mesh.variables['e3t_0'][0,:,:,:],0)
#
#     tmask = b.variables['Bathymetry'][:, :]

v_f = {#"NK1-COH": path_ecospace_out + "EcospaceMapBiomass-NK1-COH-{}.asc",
       # "NK2-CHI": path_ecospace_out + "EcospaceMapBiomass-NK2-CHI-{}.asc",
       # "NK3-FOR": path_ecospace_out + "EcospaceMapBiomass-NK3-FOR-{}.asc",
       # "ZF1-ICT": path_ecospace_out + "EcospaceMapBiomass-ZF1-ICT-{}.asc",
       # "ZC1-EUP": path_ecospace_out + "EcospaceMapBiomass-ZC1-EUP-{}.asc",
       # "ZC2-AMP": path_ecospace_out + "EcospaceMapBiomass-ZC2-AMP-{}.asc",
       # "ZC3-DEC": path_ecospace_out + "EcospaceMapBiomass-ZC3-DEC-{}.asc",
       # "ZC4-CLG": path_ecospace_out + "EcospaceMapBiomass-ZC4-CLG-{}.asc",
       # "ZC5-CSM": path_ecospace_out + "EcospaceMapBiomass-ZC5-CSM-{}.asc",
       # "ZS1-JEL": path_ecospace_out + "EcospaceMapBiomass-ZS1-JEL-{}.asc",
       # "ZS2-CTH": path_ecospace_out + "EcospaceMapBiomass-ZS2-CTH-{}.asc",
       # "ZS3-CHA": path_ecospace_out + "EcospaceMapBiomass-ZS3-CHA-{}.asc",
       # "ZS4-LAR": path_ecospace_out + "EcospaceMapBiomass-ZS4-LAR-{}.asc",
       # "PZ1-CIL": path_ecospace_out + "EcospaceMapBiomass-PZ1-CIL-{}.asc",
       # "PZ2-DIN": path_ecospace_out + "EcospaceMapBiomass-PZ2-DIN-{}.asc",
       # "PZ3-HNF": path_ecospace_out + "EcospaceMapBiomass-PZ3-HNF-{}.asc",
       "PP1-DIA": path_ecospace_out + "EcospaceMapBiomass-PP1-DIA-{}.asc",
       # "PP2-NAN": path_ecospace_out + "EcospaceMapBiomass-PP2-NAN-{}.asc",
       # "PP3-PIC": path_ecospace_out + "EcospaceMapBiomass-PP3-PIC-{}.asc", #time step format eg: 00620,
      }

# read the data into an xarray object which is versatile
yr_strt = 1980
yr_end = 2018
mo_strt = 1
da_strt = 2
mo_end = 12
da_end = 30 #unused?


######### ECOSPACE OUT NC FILE ############
ecospace_file = ecospace_code + "_" + str(filenm_yr_strt) + "-" + str(yr_end) + ".nc"
ecospace_nc = os.path.join(path_ecospace_out, ecospace_file)
ds = xr.open_dataset(ecospace_nc)
print(ds)


# water mask
# Extract the latitude and longitude arrays
lat = ds['lat'].values
lon = ds['lon'].values
dep = ds['depth'].values

# Initialize an empty array for the mask (see bloom script if need regional mask)
# combined_mask = np.zeros(lat.shape, dtype=bool)
depth_mask = dep > 0
# combined_mask &= depth_mask

# Create a new xarray Dataset with the mask
mask_ds = xr.Dataset(
    {
        # 'mask': (('row', 'col'), combined_mask)
    'mask': (('row', 'col'), depth_mask)
    },
    coords={
        'lat': (('row', 'col'), lat),
        'lon': (('row', 'col'), lon)
    }
)
# Save the mask to a new NetCDF file
#mask_ds.to_netcdf('..//..//data//evaluation//suchy_ecospace_mask.nc')
mask = mask_ds['mask'].values


v = list(v_f)[0] # hard coding to DIA

create_maps = True

# will loop over data to get true range
max_val = 0
min_val = 100000

# loop all years get the full range of phyto data
# before determining thresholds used for 'bloom' visual

print('finding max and min values for plotting vis')
skip = True
if not skip:
    for year in range(yr_strt, yr_end + 1):
        ############################################
        # pull model results for 'region' (mask corresponding to suchy)
        print(year)
        year_data = ds.sel(time=str(year))
        masked_ds = year_data.where(mask)
        if log_transform:
            masked_ds[v] = np.log(masked_ds[v]+0.01)
        masked_var = masked_ds[v]
        max_var_year = np.nanmax(masked_var.data)
        min_var_year = np.nanmin(masked_var.data)

        if max_var_year > max_val:
            max_val = max_var_year
        if min_var_year < min_val:
            min_val = min_var_year

print("done: ")
print(min_val, max_val)
# HARDCODE
min_val = 0
max_val = 3.2

# HARDCODE
yr_strt = 2009
yr_end = 2009


# Loop over the years and select data based on the 'time' coordinate
timestamps = []

for year in range(yr_strt, yr_end + 1):
    ############################################
    # pull model results for 'region' (mask corresponding to suchy)
    print(year)
    year_data = ds.sel(time=str(year))
    masked_ds = year_data.where(mask)
    if log_transform:
        masked_ds[v] = np.log(masked_ds[v]+0.01)
    masked_var = masked_ds[v]

    # select time period within each year (one per timestep)
    # yr_plt_strt = 2005
    # yr_plt_end = 2005
    mo_strt = 2
    mo_end = 2
    da_srt = 1
    da_end = 28

    # Select time period to plot figs (one per timestep)
    start_date = pd.Timestamp(year=year, month=mo_strt, day=da_srt)
    end_date = pd.Timestamp(year=year, month=mo_end, day=da_end)

    # Get the timestep start, end from dates
    time_coords = masked_ds.time.to_index()
    # ts_strt = time_coords.get_loc(start_date, method='bfill') # deprecated
    # ts_end = time_coords.get_loc(end_date, method='ffill') + 1
    ts_strt = time_coords.get_indexer([start_date], method='bfill')[0]
    ts_end = time_coords.get_indexer([end_date], method='ffill')[0] + 1

    map_num = 0
    for ts in range(ts_strt, ts_end):
    # for ts in range(ts_strt - 1, ts_end - 1):

        # print(masked_var['lat'][50,50])
        # print(masked_var['lon'][50,50])
        # print(masked_var['dep'][50,50])
        # print(dep[50,50])
        v1 = masked_var # to further develop, to loop over v's
        d1 = v1[ts]


        timestamp = pd.Timestamp(d1.time.values)
        timestamps.append(timestamp)
        year1 = timestamp.year
        month1 = timestamp.month
        day1 = timestamp.day

        # if create_maps and map_num == 0: # remove map_num if you want map per time step
        if create_maps:
            ################### FIGS PER TIME STEP #######################
            # get the ecospace output indices for start and end from desired date range to plot
            # test map of ECOSPACE domain by row /col
            # note the bottom left should be fixed (nans)
            fig, axs = plt.subplots(figsize=(5, 8))

            # bathy
            # c = plt.contourf(  # navlon, navlat,
            #     # dep.data, 25, vmin=4, vmax=430, cmap=cmocean.cm.ice_r,
            #     # NOTE MINUS VALS TO HACK COLOR
            #     dep.data, 50, vmin=-600, vmax=500, cmap=cmocean.cm.ice_r,
            #     # transform=crt.crs.PlateCarree(),
            #     zorder=1)

            # ALT WAY TO COLOR WATER ETC
            # Create a mask: depth > 0 is dark blue, depth <= 0 is light grey
            masked_deps = np.where(dep <= 0, 1, np.nan)  # 1 for dark blue (depth > 0), 0 for light grey (depth <= 0)
            masked_deps = np.where(dep <= 0, 1, 0)
            # Create a custom colormap: light grey and dark blue


            # WHY ARE THE MODEL AND DEPTH MAPS MISALIGNED
            # AND SHIFTED BY ONE PIXEL TO THE LEFT AND DOWN
            # Is it happening in contour or prep for pcolormesh?
            # or is it a bug when creating the Ecospace NC file due to indexing issue?

            print(d1['depth'][50,50])
            print(d1['lat'][50,50])
            print(d1['lon'][50,50])
            print(dep[50,50])


            # cmap = ListedColormap(['lightgrey', 'darkblue'])
            cmap = ListedColormap(['darkblue', 'lightgrey'])
            # Plot using contourf
            c = plt.contourf(
                masked_deps.data, levels=1, cmap=cmap,
                # v1['dep'], levels=1, cmap=cmap,
                zorder=1
            )


            # custom colormap to emphasise lower values
            colors = [(1, 1, 1), (0.8, 0.8, 0.8), (0.6, 0.6, 0.6), (0.4, 0.4, 0.4), (0, 0, 1), (0, 1, 0), (1, 1, 0),
                      (1, 0, 0)]
            # sc1
            # nodes = [0.0, 0.1, 0.2, 0.35, 0.6, 0.8, 0.9, 1.0] # Adjust these values to emphasize certain ranges
            # sc2
            nodes = [0.0, 0.15, 0.25, 0.3, 0.5, 0.65, 0.75, 1.0]
            custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))

            # specific month coloring experiments
            nodes = [0.0, 0.15, 0.25, 0.3, 0.69, 0.77, 0.86, 1.0]
            custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))

            # PPT colors, custom
            colors = [(1, 1, 1), (0.8, 0.8, 0.8), (0.6, 0.6, 0.6), (0.4, 0.4, 0.4), (0, 0, 1), (0, 1, 0), (1, 1, 0)]
            nodes = [0.0, 0.15, 0.25, 0.3, 0.45, 0.85, 1.0]
            custom_cmap = LinearSegmentedColormap.from_list("custom_cmap", list(zip(nodes, colors)))

            # sc 1
            vmin = 0
            vmax = 3

            # dump vals not within certain range? experiment left here
            dump_vals = False
            if dump_vals:
                threshold = 2
                d1 = d1.where(d1 > threshold, np.nan)

            # w = plt.pcolormesh(d1.data,
            #                    cmap=custom_cmap,
            #                    vmin=vmin, vmax=vmax,
            #                    # alpha=0.5,
            #                    zorder=3)

            # want effect where transparency changes depending on values
            # There seems to be a way.. https://stackoverflow.com/questions/17170229/setting-transparency-based-on-pixel-values-in-matplotlib
            # I'm now thinking though I will just use a single color like NASA anim
            vmin, vmax = min_val, max_val
            alpha_min_value = 0.5 * max_val # Value where alpha = 0 (fully transparent)
            alpha_max_value = max_val
            # Clip the data to ensure it's within the alpha bounds
            clipped_data = np.clip(d1.data, alpha_min_value, alpha_max_value)
            # Normalize the clipped data to be between 0 and 1 for alpha scaling
            normalised_data = (clipped_data - alpha_min_value) / (alpha_max_value - alpha_min_value)
            normalised_data[d1.data < alpha_min_value] = 0  # Fully transparent below alpha_min_value
            normalised_data[d1.data > alpha_max_value] = 1
            # print(np.nanmin(normalised_data), np.nanmax(normalised_data))
                # normalised_data = (d1.data - vmin) / (vmax - vmin) # normalised data
            # alpha_data = 1 - normalised_data
            alpha_data = normalised_data

            single_color = [0.0, 1.0, 0.0, 1.0]  # RGB for green with alpha = 1 (fully opaque)
            rgba_data = np.ones((d1.shape[0], d1.shape[1], 4)) * single_color  # Create an array of blue color
            rgba_data[..., -1] = alpha_data  # Set the alpha channel to vary based on data
            w = plt.pcolormesh(rgba_data, shading='auto', zorder=3)


            title = ecospace_code + "-" + str(year) + "-" + str(month1) + "-" + str(day1) + "-" + str(ts)
            title = str(year) + "-" + str(month1) + "-" + str(day1) + "-" + str(ts)
            axs.set_title('{v} - {d}'.format(v=v, d=title))

            # plt.colorbar(w, ax=axs)
            axs.invert_yaxis()

            # plt.show()

            # Save the figure in different formats
            month01 = buildSortableString(month1, 2)
            day01 = buildSortableString(day1, 2)
            # out_fig_name = "..//figs//{}_{}-{}{}{}-{}.{}"
            out_fig_name = "..//..//figs//TEMP_{}_{}-{}{}{}-{}.{}"
            # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "PDF"), format='pdf')
            out_fig_name = out_fig_name.format(ecospace_code, v, year, month01, day01, ts, "PNG")
            fig.savefig(out_fig_name, format='png')
            # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "EPS"), format='eps')
            print("Saved fig " + out_fig_name)
            map_num += 1

    exit()