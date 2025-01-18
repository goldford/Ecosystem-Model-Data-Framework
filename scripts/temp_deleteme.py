
import sys
sys.path.append("C://Users//Greig//Documents//github//HOTSSea_v1_NEMOandSupportCode//desktop//code//manuscript_figs")
# noinspection PyUnresolvedReferences
from GO_tools import buildSortableString, read_sdomains # last to from GO_helpers
import netCDF4 as nc
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.path import Path as Path
matplotlib.use('TkAgg') # https://matplotlib.org/stable/users/explain/figure/backends.html



# ////////////////////// PATHS ////////////////////////
# local test p's
nc_p = "C:/Users/Greig/Downloads/NEMO_out/" # double fwd slashes don't work with os path join
# eg SalishSea1500-RUN216_MonthlyMean_grid_T_2D_y1980m10.nc
project_p = ""
tmask_p = "C:/Users/Greig/Documents/GitHub/Ecosystem-Model-Data-Framework/data/mesh mask"
tmask_f = "mesh_mask_20210406.nc"
bathy_p = "..//..//data//bathymetry//bathy_salishsea_1500m_20210706.nc"
out_p = ""

frmask_p = "..//..//data//"
frmask_f = "FraserRiverMaskPnts.yml"
frmask_fullp = os.path.join(frmask_p, frmask_f)

# server paths
# nc_p = "DATA/SS1500-RUN216/NEMO_monthly_NC/"
# project_p = "/project/6006412/goldford/ECOSPACE/"
# tmask_p = "..//..//data//mesh mask//mesh_mask_20210406.nc"
# out_p = project_p + "DATA/SS1500-RUN{modelrun}/NEMO_monthly_NC_pacea/"


# get masks, z levs, widths
with nc.Dataset(os.path.join(tmask_p, tmask_f), 'r') as mesh:
    # shapes of (1, 40, 299, 132)
    # tmask - land mask in 3D, t pnt
    # e3t0 - depth bin widths
    # gdept_0 - centres of the depth levs
    tmask = mesh.variables['tmask'][:]
    e3t0 = mesh.variables['e3t_0'][:]
    gdept_0 = mesh.variables['gdept_0'][:]
    gdepw_0 = mesh.variables['gdepw_0'][:]
    navlat = mesh.variables['nav_lat'][:]
    navlon = mesh.variables['nav_lon'][:]

# Select the required slice
gdept_slice = gdept_0[0, 0, :, :]  # Shape: (299, 132)

# Plot as a heatmap
plt.figure(figsize=(10, 8))
plt.pcolormesh(navlon, navlat, gdept_slice, shading='auto')
plt.colorbar(label='Depth (gdept_0)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Depth Variation at gdept_0[:, 0, :, :]')
plt.show()

