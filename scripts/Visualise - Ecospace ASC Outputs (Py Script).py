# by GO 2024-05
# basing it on experiments in notebook
# 'Visuals - Analysis SoG Vars Asc'.ipynb in NEMO git

import numpy as np
import xarray as xr
import pandas as pd
import os
import netCDF4 as nc
from matplotlib.colors import LinearSegmentedColormap

import cartopy as cp
import matplotlib.pyplot as plt
from matplotlib import patches
import mpl_toolkits
from mpl_toolkits.basemap import Basemap
import cartopy
from cartopy import crs, feature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.patches import Rectangle
import cmocean as cm
import cartopy.crs as ccrs
from helpers import buildSortableString, is_leap_year
from datetime import datetime, timedelta

# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//ECOSPACE 216 2024 - 10yr1yr 2003-2018//asc//"
#ecospace_code = "Scv1-NoMultMix"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v1 - NoMultMixing - 3DAY 1yr10yr//asc//"
# ecospace_code = "Scv2-MultMix"
# path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v2 - MultMixing - 3DAY 1yr10yr//asc//"
#ecospace_code = "Scv3-MultxPAR"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v3 - MixingxPAR - 3DAY 1yr10yr//asc//"
#ecospace_code = "Scv4-MultxPARlimZ"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v4 - MixingxPARLimitZ - 3DAY 1yr10yr//asc//"
ecospace_code = "Scv5-PARMixingNut90"
#path_ecospace_out= "C://Users//Greig//Documents//EwE output//ECOSPACE_diag_LTL_v6-7-0-18060_2024_Carb_36day3day_v11//Sc216 v5 - PARMixingNut90 - 3DAY 1yr10yr//asc//"
# path_ecospace_out = "..//data//ecospace_out//"
path_ecospace_out = "C://Users//Greig//Sync//PSF//EwE//Georgia Strait 2021//LTL_model//ECOSPACE_OUT"

# mxng_p = "NEMO_prepped_as_ASC/{var}/"
# tmp_p = "NEMO_prepped_as_ASC/{var}/"
# li_p = "ECOSPACE_in_PAR3_Sal10m/{var}/"
# k_p = "ECOSPACE_in_PAR3_Sal10m/RUN203_{var}/"
# sal_p = "NEMO_prepped_as_ASC/{var}/"
# path_data2 = "../data/forcing/"
# li_p2 = "RDRS_light_monthly_ASC/{var}/"
# wi_p = "RDRS_wind_monthly_ASC/{var}/"

# template file names w/ var names

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
       "PP1-DIA": path_ecospace_out + "EcospaceMapBiomass-PP1-DIA-{}.asc"#,
       # "PP2-NAN": path_ecospace_out + "EcospaceMapBiomass-PP2-NAN-{}.asc",
       # "PP3-PIC": path_ecospace_out + "EcospaceMapBiomass-PP3-PIC-{}.asc", #time step format eg: 00620,
      }


# read the data into an xarray object which is versatile
yr_strt = 2003
yr_end = 2018
mo_strt = 1
da_strt = 2
mo_end = 12
da_end = 30 #unused?


################## GET THE NC FILE #########################
ecospace_file = ecospace_code + "_" + str(yr_strt) + "-" + str(yr_end) + ".nc"
ecospace_nc = os.path.join(path_ecospace_out, ecospace_file)
print(ecospace_nc)
ds_loaded = xr.open_dataset(ecospace_nc)

print(ds_loaded)

# load bathy and tmask for plotting
with nc.Dataset('..//data//basemap//bathy_salishsea_1500m_20210706.nc') as bath:
    navlon = bath.variables['nav_lon'][:]
    navlat = bath.variables['nav_lat'][:]
    bathy = bath.variables['Bathymetry'][:]
#     bathy=np.sum(mesh.variables['bathymetry'][0,:])#*mesh.variables['e3t_0'][0,:,:,:],0)

with nc.Dataset('..//data//basemap//mesh_mask_20210406.nc') as mesh:
    #     print(mesh.variables)
    tmask = mesh.variables['tmask'][:]

grid_p = '..//data/basemap//'
grid_f = 'bathy_salishsea_1500m_20210706.nc'
grid = xr.open_dataset(os.path.join(grid_p, grid_f), mask_and_scale=False)

################## MAP TEMPLATE ############################
def make_map(ax, grid, w_map=[-124, -123.9, 47.7, 50.6],
             rotation=39.2,
             par_inc=0.25,
             mer_inc=0.5,
             fs=14,
             bg_color='#e5e5e5'
             ):
    """
    """
    # fcolor='burlywood'
    fcolor = 'white'
    # Make projection
    m = Basemap(ax=ax,
                projection='lcc', resolution='c',
                lon_0=(w_map[1] - w_map[0]) / 2 + w_map[0] + rotation,
                lat_0=(w_map[3] - w_map[2]) / 2 + w_map[2],
                llcrnrlon=w_map[0], urcrnrlon=w_map[1],
                llcrnrlat=w_map[2], urcrnrlat=w_map[3])

    # Add features and labels
    x, y = m(grid.nav_lon.values, grid.nav_lat.values)
    ax.contourf(x, y, grid.Bathymetry, [-0.01, 0.01], colors=fcolor)
    ax.contour(x, y, grid.Bathymetry, [-0.01, 0.01], colors='black', linewidths=0.1)
    ax.contourf(x, y, grid.Bathymetry, [0.011, 500], colors=bg_color)

    #     m.drawmeridians(np.arange(-125.5, -122, mer_inc), labels=[0, 0, 0, 1], linewidth=0.2, fontsize=fs)
    m.drawmeridians(np.arange(-125.5, -122, mer_inc), labels=[0, 0, 0, 0], linewidth=0.2, fontsize=fs)
    #     m.drawparallels(np.arange(48, 51, par_inc), labels=[1, 0, 0, 0], linewidth=0.2, fontsize=fs)
    m.drawparallels(np.arange(47, 51, par_inc), labels=[0, 0, 0, 0], linewidth=0.2, fontsize=fs)

    return m


def adjust_map(ax,
               lat_bl=47.8,
               lon_bl=-123.2,
               lat_br=48.8,
               lon_br=-122.28,
               lat_tl=50.3,
               lon_tl=-124.75,
               lat_bl2=48.2,
               lon_bl2=-123.5,
               label_grid=False
               ):
    # fcolor='burlywood'
    fcolor = 'white'
    # set width using map units
    # bottom left
    x_bl, y_bl = m(lon_bl, lat_bl)
    x_br, _ = m(lon_br, lat_br)
    ax.set_xlim(x_bl, x_br)

    # top left
    x_tl, y_tl = m(lon_tl, lat_tl)
    x_bl, y_bl = m(lon_bl2, lat_bl2)
    ax.set_ylim(y_bl, y_tl)

    # fix a little path in bottom right
    lccx_TL, lccy_TL = m(-122.83, 49.4)
    lccx_BR, lccy_BR = m(-122.58, 48.7)
    lccx_BL, lccy_BL = m(-122.33, 48.7)
    lccw = lccx_BL - lccx_BR
    lcch = lccy_TL - lccy_BL

    ax.add_patch(patches.Rectangle(
        (lccx_BL, lccy_BL), lccw, lcch,
        facecolor=fcolor, edgecolor='k',
        linewidth=0,
        zorder=0))

    if label_grid:
        #         rotation = 0
        #         ax.annotate("49.5 N", xy=(m(-124.8, 49.5)), xytext=(m(-125, 49.5)),
        #                     xycoords='data', textcoords='data',
        #                     ha='right', va='center', fontsize=8, rotation=rotation)
        fsmer = 7.5
        ax.annotate("50.5N", xy=(m(-123.99, 50.43)), xytext=(m(-123.95, 50.43)),
                    xycoords='data', textcoords='data',
                    ha='left', va='center', fontsize=fsmer, rotation=-30)

        ax.annotate("123.5W", xy=(m(-123.65, 50.05)), xytext=(m(-123.65, 50.15)),
                    xycoords='data', textcoords='data',
                    ha='left', va='center', fontsize=fsmer, rotation=60)
        ax.annotate("50N", xy=(m(-123.5, 49.94)), xytext=(m(-123.42, 49.94)),
                    xycoords='data', textcoords='data',
                    ha='left', va='center', fontsize=fsmer, rotation=-30)
        ax.annotate("49.5N", xy=(m(-123.06, 49.4)), xytext=(m(-122.9, 49.4)),
                    xycoords='data', textcoords='data',
                    ha='left', va='center', fontsize=fsmer, rotation=-30)

        ax.annotate("122.5W", xy=(m(-122.8, 49.15)), xytext=(m(-122.61, 49.15)),
                    xycoords='data', textcoords='data',
                    ha='left', va='center', fontsize=fsmer, rotation=60)

        ax.annotate("49N", xy=(m(-122.65, 49)), xytext=(m(-122.4, 48.93)),
                    xycoords='data', textcoords='data',
                    ha='left', va='center', fontsize=fsmer, rotation=-30)


#         ax.annotate("", xy=(m(-123.5, 50)), xytext=(m(-123.4, 50)),
#                     xycoords='data', textcoords='data',
#                     ha='left', va='center', fontsize=8, rotation=rotation)

#         rotation = -30
#         ax.annotate("50 N", xy=(m(-124.3, 50)), xytext=(m(-124.3, 50.08)),
#                     xycoords='data', textcoords='data',
#                     ha='right', va='center', fontsize=8, rotation=rotation)
#         ax.annotate("49 N", xy=(m(-124.3, 49)), xytext=(m(-124.3, 49.08)),
#                     xycoords='data', textcoords='data',
#                     ha='right', va='center', fontsize=8, rotation=rotation)


def custom_formatter(x, pos):
    return f'{x:.2f}'


def custom_formatter2(x, pos):
    return f'{x:.1f}'

# test grey rotated map
fig, ax = plt.subplots(figsize=(5,8)) #w,h

m = make_map(ax,grid,
             w_map=[-130, -115, 41, 60],
             rotation=36,
             par_inc=0.5,
             mer_inc=1
            )

adjust_map(ax, lat_bl=46.8, lon_bl=-123.3,
           lat_br=47.8, lon_br=-121.7,
           lat_tl=50.9, lon_tl=-124.75,
           lat_bl2=47.2, lon_bl2=-123.3,
           label_grid=False
              )
# adjust_map(ax)
print("xy lats lons don't look right here but do in panel")
plt.show()
# fig.savefig("test_map.png", format='png')







################### FIGS PER TIME STEP #######################
# get the ecospace output indices for start and end from desired date range to plot
# test map of ECOSPACE domain by row /col
# note the bottom left should be fixed (nans)
# select time period to plot figs (one per timestep)
yr_plt_strt = 2005
yr_plt_end = 2005
mo_plt_strt = 1
mo_plt_end = 12
da_plt_srt = 1
da_plt_end = 30

# Select time period to plot figs (one per timestep)
start_date = pd.Timestamp(year=yr_plt_strt, month=mo_plt_strt, day=da_plt_srt)
end_date = pd.Timestamp(year=yr_plt_end, month=mo_plt_end, day=da_plt_end)

ds = ds_loaded

# Get the timestep start, end from dates
time_coords = ds.time.to_index()
# ts_strt = time_coords.get_loc(start_date, method='bfill') # deprecated
# ts_end = time_coords.get_loc(end_date, method='ffill') + 1
ts_strt = time_coords.get_indexer([start_date], method='bfill')[0]
ts_end = time_coords.get_indexer([end_date], method='ffill')[0] + 1

#ts_strt = 720 # model timestep start
#ts_end = 740  # model timestep end


for ts in range(ts_strt-1, ts_end-1):
    fig, axs = plt.subplots(figsize=(5, 8))

    v = list(v_f)[0]
    v1 = ds[v]
    d1 = v1[ts]


    timestamp = pd.Timestamp(d1.time.values)

    year1 = timestamp.year
    month1 = timestamp.month
    day1 = timestamp.day



    # custom colormap to emphasise lower values
    colors = [(1, 1, 1), (0.8, 0.8, 0.8), (0.6, 0.6, 0.6), (0.4, 0.4, 0.4), (0, 0, 1), (0, 1, 0), (1, 1, 0), (1, 0, 0)]
    # sc1
    #nodes = [0.0, 0.1, 0.2, 0.35, 0.6, 0.8, 0.9, 1.0] # Adjust these values to emphasize certain ranges
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

    w = plt.pcolormesh(d1.data, cmap=custom_cmap, vmin=vmin, vmax=vmax)
    title = ecospace_code + "-" + str(year1) + "-" + str(month1) + "-" + str(day1) + "-" + str(ts)
    title = str(year1) + "-" + str(month1) + "-" + str(day1) + "-" + str(ts)
    axs.set_title('{v} - {d}'.format(v=v, d=title))
    plt.colorbar(w, ax=axs)
    axs.invert_yaxis()
    # plt.show()

    # Save the figure in different formats
    month01 = buildSortableString(month1,2)
    day01 = buildSortableString(day1, 2)
    out_fig_name = "..//figs//{}_{}-{}{}{}-{}.{}"
    out_fig_name = "..//figs//EXPERIMENT_{}_{}-{}{}{}-{}.{}"
    # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "PDF"), format='pdf')
    out_fig_name = out_fig_name.format(ecospace_code, v, year1, month01, day01, ts, "PNG")
    fig.savefig(out_fig_name, format='png')
    # fig.savefig(out_fig_name.format(ecospace_code, v, ts, "EPS"), format='eps')
    print("Saved fig " + out_fig_name)






############### ANIMATION ################
# not done yet
do_anim = False
if not do_anim:
    exit(1)

from pathlib import Path
import matplotlib.animation as animation

# must download and point matplotlib to .exe
ffmpeg_folder = Path(
    "C://Users//Greig//Documents//github//NEMO-Salish-Sea-1500//code//ffmpeg-2023-04-03-git-6941788d24-essentials_build//bin//ffmpeg.exe")
plt.rcParams['animation.ffmpeg_path'] = ffmpeg_folder


def _update_plot(i, fig, clrmsh):
    # i = (i * 12) + 5 # if only displaying a single month each yr
    i = i + 420
    frm_data = v1[i]
    # clrmsh = plt.pcolormesh(frm_data, cmap = cm.cm.balance, vmin=min_value, vmax=max_value)
    clrmsh = plt.pcolormesh(frm_data, cmap=cm.cm.balance)
    title = u"%s — %s" % (v, str(v1.time[i].values)[:7])
    ax.set_title(title)
    return clrmsh


for v in v_f:
    # PAR-VarZ-VarK
    if (v == "PP1-DIA"):
        v1 = ds[v]
        frm_data = v1[6]
        print(v1)
        # frames = v1.time.size
        frames = 96  # override
        fps = 4

        min_value = 0.1
        max_value = 2

        fig = plt.figure(figsize=(5, 8))
        ax = fig.add_subplot()
        ax.invert_yaxis()
        title = u"%s — %s" % (v, str(v1.time[0].values)[:7])
        ax.set_title(title)

        # clrmsh = plt.pcolormesh(frm_data, cmap = cm.cm.balance, vmin=min_value, vmax=max_value)
        clrmsh = plt.pcolormesh(frm_data, cmap=cm.cm.balance)
        plt.colorbar(clrmsh, ax=ax)

        anim = animation.FuncAnimation(fig, _update_plot, fargs=(fig, clrmsh),
                                       frames=frames, interval=1)

        anim.save('{v}.mp4'.format(v=v), writer=animation.FFMpegWriter(fps=fps))

# # get the day of year in 3DAY BLOCKS
# # while dumping the remainder (blocks 121,122) b/c this is how 3D model driver data prepped
# # ie 1 year = '10 years' thus 120 time steps per real year, 5-6 days remainder
# days_of_year = range(2, 360, 3)
# date_list = []
# months = list(range(1, 12))
# pd_timestamps = []
# time_step_model = []
# i = 0
# for yr in range(yr_strt, yr_end+1):
#     for doy in days_of_year:
#         date1 = datetime(yr, 1, 1) + timedelta(days=doy - 1)
#         month1 = date1.month
#         if yr == yr_end and (month1 == mo_end + 1):
#             break
#         day_of_month = date1.day
#         date_list.append([yr, month1, day_of_month, doy])
#         pd_timestamps.append(pd.Timestamp(yr, month1, day_of_month))
#         time_step_model.append(buildSortableString(i+1, 5))
#         i += 1
#     if yr == yr_end:
#         break


# https://pandas.pydata.org/docs/user_guide/timeseries.html#timeseries-offset-aliases
# time = pd.date_range(start='{yr_strt}-{mo_strt}-{da_strt}'.format(yr_strt=yr_strt, mo_strt=mo_strt, da_strt=da_strt),
#                               end='{yr_end}-{mo_end}-{da_end}'.format(yr_end=yr_end, mo_end=mo_end, da_end=da_end),
#                               freq='3D')
# for t in range(1,len(pd_timestamps)+1):
#     time_step_model.append(buildSortableString(t,5))


# empty ds for all
# ds = xr.Dataset(
#     coords={
#         'time': pd_timestamps,
#         'row': range(1, rows + 1),
#         'col': range(1, cols + 1)
#     },
#     attrs={
#         'description': 'dataset of monthly ASC files',
#     }
# )


# create empty variable with correct shape
# for v in v_f:
#     ds[v] = xr.DataArray(
#         np.nan * np.zeros((len(pd_timestamps), rows, cols)),
#         dims=('time', 'row', 'col'),
#         attrs={'description': f'{v} data'}
#     )

# load the data
# for v in v_f:
#     attribute = v_f[v]
#     for y in range(yr_strt, yr_end + 1):
#         for m in sorted(months_d.keys()):
#             f_n = v_f[v].format(var=v, year=y, month=m)
#             with open(f_n) as f:
#                 data = np.loadtxt(f, skiprows=skiprows)
#
#                 # homogenize what nans are
#                 data[data == -9999.0] = ['nan']
#                 data[data == 0.0] = ['nan']
#
#                 # fix issue with bottom left area in map
#                 data[140:, :15] = ['nan']
#
#                 ds[f'{v}'.format(var=v)].loc[{'time': f'{y}-{m}-15'.format(year=y, month=m)}] = xr.DataArray(
#                     data,
#                     dims=('row', 'col'),
#                     attrs={'description': f'{v} data for year {y} month {m}'.format(var=v, year=y, month=m)}
#                 )

# load these ASCs into a nice xarray dataset
# for v in v_f:
#     attribute = v_f[v]
#     for t in range(0, len(time_step_model)):
#         f_n = v_f[v].format(time_step_model[t])
#         ti = pd_timestamps[t]
#         yr = ti.year
#         mo = ti.month
#         da = ti.day
#
#         with open(f_n) as f:
#             data = np.loadtxt(f, skiprows=skiprows)
#
#             # homogenize what nans are
#             data[data == -9999.0] = ['nan']
#             data[data == 0.0] = ['nan']
#
#             # fix issue with bottom left area in map
#             data[140:, :15] = ['nan']
#
#             ds[f'{v}'.format(var=v)].loc[{'time': f'{yr}-{mo}-{da}'.format(year=yr, month=mo, day=da)}] = xr.DataArray(
#                 data,
#                 dims=('row', 'col'),
#                 attrs={'description': f'{v} data for year {yr} month {mo} day {da}'.format(var=v, year=yr, month=mo, day=da)}
#             )
