# Created: July 15, 2024
# Author: G Oldford
# Purpose: Compare Ecospace and SSCast outputs
#          at monthly or seasonal scales, by phytoplankton 'provinces' (Tereza's)
# Source of Data:
#  - SalishSeaCast phyto data downloaded from url below using another script
#    https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DBiologyFields1moV19-05.html
#
import os
import numpy as np
import pandas as pd
import xarray as xr
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib.path import Path
from datetime import datetime
from helpers import read_sdomains

# File paths
ssc_f = "SalishSeaCast_biology_{}.nc"
ssc_p = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL"
ecospace_f = "Scv7-PARMixingNut90Temp_2003-2018.nc"
ecospace_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT"
ecospace_map_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"
ecospace_map_f = "Ecospace_grid_20210208_rowscols.csv"
domain_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
domain_f = "analysis_domains_jarnikova.yml"
meshm_f = 'ubcSSn3DMeshMaskV17-02.nc'
meshm2D_f = 'ubcSSn2DMeshMaskV17-02.nc'
bathy_f = 'ubcSSnBathymetryV21-08_a29d_efc9_4047.nc'

# Load subdomains
sdomains = read_sdomains(os.path.join(domain_p, domain_f))

# Time range
yr_st = 2008
yr_en = 2018

def get_sscast_data(file_path):
    if os.path.exists(file_path):
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                metadata = {
                    "Dimensions": {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions},
                    "Variables": {var: dataset.variables[var].dimensions for var in dataset.variables}
                }
                time_var = dataset.variables['time']
                time_units = time_var.units
                time_calendar = time_var.calendar if 'calendar' in time_var.ncattrs() else 'standard'
                times = nc.num2date(time_var[:], units=time_units, calendar=time_calendar)

                ciliates = dataset.variables['ciliates'][:]
                flagellates = dataset.variables['flagellates'][:]
                diatoms = dataset.variables['diatoms'][:]

            return metadata, times, ciliates, flagellates, diatoms
        except OSError as e:
            return {"Error": str(e)}, None, None, None, None
    else:
        return {"Error": "File not found."}, None, None, None, None

def process_year(year):
    metadata, times, ssc_cil, ssc_fla, ssc_dia = get_sscast_data(os.path.join(ssc_p, ssc_f.format(year)))
    if metadata.get("Error"):
        print(f"Error processing year {year}: {metadata['Error']}")
        return None, None, None

    with nc.Dataset(os.path.join(ssc_p, meshm_f)) as mesh:
        tmask = mesh.variables['tmask'][:]  # 0's where depth lev exceeds
        e3t0 = mesh.variables['e3t_0'][:]  # 'widths' or 'weights' for each depth lev
        gdept_0 = mesh.variables['gdept_0'][:] # depth under t grid points

    sum_e3t = np.nansum(e3t0[:, 0:19, :, :], axis=1)  # 19 because that's depth max of phyto in model
    e3t0_weights = np.where(sum_e3t != 0, e3t0[:, 0:19, :, :] / sum_e3t[np.newaxis, :], 0)

    weighted_dia = ssc_dia * e3t0_weights
    weighted_cil = ssc_cil * e3t0_weights
    weighted_fla = ssc_fla * e3t0_weights

    dep_avg_dia = np.sum(weighted_dia, axis=1)  # shape should now be 12, 898, 398
    dep_avg_fla = np.sum(weighted_fla, axis=1)
    dep_avg_cil = np.sum(weighted_cil, axis=1)

    return dep_avg_dia, dep_avg_fla, dep_avg_cil

def compute_anomalies(data):
    valid_data = data[np.isfinite(data)]
    mean_val = np.mean(valid_data)
    std_val = np.std(valid_data)

    anomalies = np.zeros_like(data)
    mask = np.isfinite(data)
    anomalies[mask] = (data[mask] - mean_val) / std_val

    return anomalies

def compute_subdomain_anomalies(data, subdomain_array, subdomain_name):
    mask = (subdomain_array != '') & np.isfinite(data)
    valid_data = data[mask]
    mean_val = np.nanmean(valid_data)
    std_val = np.nanstd(valid_data)
    print(subdomain_name + ":")
    print("mean over all years")
    print(mean_val)
    print("min over all years")
    print(np.nanmin(valid_data))

    mask = (subdomain_array == subdomain_name) & np.isfinite(data)
    valid_data = data[mask]

    anomalies = np.full_like(data, np.nan)
    anomalies[mask] = (valid_data - mean_val) / std_val

    return anomalies

def create_boxplots(data, title, save_path):
    months = np.arange(1, 13)
    data_by_month = [data[month - 1, :, :].flatten() for month in months]

    plt.figure(figsize=(12, 6))
    plt.boxplot(data_by_month, labels=months)
    plt.title(f'Average Anomalies Boxplot for {title}')
    plt.xlabel('Month')
    plt.ylabel('Anomaly (relative to standard deviation)')
    plt.yscale('log')
    plt.grid(True)

    plt.savefig(save_path)
    plt.close()

def compute_boxplot_stats(data):
    stats = {}
    data = data[data > -999]
    if len(data) == 0:
        return None  # Return None if no valid data

    data = np.ma.masked_invalid(data)
    q1 = np.percentile(data, 25, axis=0)
    q3 = np.percentile(data, 75, axis=0)
    median = np.percentile(data, 50, axis=0)
    iqr = q3 - q1
    lower_whisker = q1 - 1.5 * iqr
    upper_whisker = q3 + 1.5 * iqr

    lower_whisker = np.maximum(lower_whisker, np.min(data, axis=0))
    upper_whisker = np.minimum(upper_whisker, np.max(data, axis=0))

    stats['q1'] = q1
    stats['q3'] = q3
    stats['median'] = median
    stats['lower_whisker'] = lower_whisker
    stats['upper_whisker'] = upper_whisker

    return stats

def get_monthly_stats(data, tmask):
    months = np.arange(1, 13)
    monthly_stats = []

    for month in months:
        month_data = data[month - 1, :, :]
        month_data = np.where(tmask == 0, np.nan, month_data)
        month_data = month_data.flatten()
        month_stats = compute_boxplot_stats(month_data)
        monthly_stats.append(month_stats)
    return monthly_stats

def plot_custom_boxplot(stats_ecospace, stats_ssc, title, save_path):
    fig, ax = plt.subplots(figsize=(15, 8))
    months = np.arange(1, 13)
    positions_ecospace = np.arange(1, 13) - 0.2  # Adjust positions slightly to the left
    positions_ssc = np.arange(1, 13) + 0.2  # Adjust positions slightly to the right
    widths = 0.35

    for i, (month_stats_ecospace, month_stats_ssc) in enumerate(zip(stats_ecospace, stats_ssc)):
        pos_ecospace = positions_ecospace[i]
        pos_ssc = positions_ssc[i]

        if month_stats_ecospace is not None:
            ax.plot([pos_ecospace, pos_ecospace], [month_stats_ecospace['lower_whisker'], month_stats_ecospace['q1']], color='k')
            ax.plot([pos_ecospace, pos_ecospace], [month_stats_ecospace['q3'], month_stats_ecospace['upper_whisker']], color='k')
            ax.plot([pos_ecospace - widths / 2, pos_ecospace + widths / 2], [month_stats_ecospace['lower_whisker'], month_stats_ecospace['lower_whisker']], color='k')
            ax.plot([pos_ecospace - widths / 2, pos_ecospace + widths / 2], [month_stats_ecospace['upper_whisker'], month_stats_ecospace['upper_whisker']], color='k')
            ax.add_patch(plt.Rectangle((pos_ecospace - widths / 2, month_stats_ecospace['q1']), widths, month_stats_ecospace['q3'] - month_stats_ecospace['q1'], fill=True, color='pink'))
            ax.plot([pos_ecospace - widths / 2, pos_ecospace + widths / 2], [month_stats_ecospace['median'], month_stats_ecospace['median']], color='k')

        if month_stats_ssc is not None:
            ax.plot([pos_ssc, pos_ssc], [month_stats_ssc['lower_whisker'], month_stats_ssc['q1']], color='k')
            ax.plot([pos_ssc, pos_ssc], [month_stats_ssc['q3'], month_stats_ssc['upper_whisker']], color='k')
            ax.plot([pos_ssc - widths / 2, pos_ssc + widths / 2], [month_stats_ssc['lower_whisker'], month_stats_ssc['lower_whisker']], color='k')
            ax.plot([pos_ssc - widths / 2, pos_ssc + widths / 2], [month_stats_ssc['upper_whisker'], month_stats_ssc['upper_whisker']], color='k')
            ax.add_patch(plt.Rectangle((pos_ssc - widths / 2, month_stats_ssc['q1']), widths, month_stats_ssc['q3'] - month_stats_ssc['q1'], fill=True, color='lightblue'))
            ax.plot([pos_ssc - widths / 2, pos_ssc + widths / 2], [month_stats_ssc['median'], month_stats_ssc['median']], color='k')

    ax.set_xticks(months)
    ax.set_xticklabels(months)
    ax.set_title(title)
    ax.set_xlabel('Month')
    ax.set_ylabel('Anomaly (relative to standard deviation)')
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
# def plot_custom_boxplot(stats, title, save_path):
#     fig, ax = plt.subplots(figsize=(12, 6))
#     positions = np.arange(1, 13)
#     widths = 0.6
#
#     for i, month_stats in enumerate(stats):
#         if month_stats is None:
#             continue
#
#         ax.plot([positions[i], positions[i]], [month_stats['lower_whisker'], month_stats['q1']], color='k')
#         ax.plot([positions[i], positions[i]], [month_stats['q3'], month_stats['upper_whisker']], color='k')
#         ax.plot([positions[i] - widths/2, positions[i] + widths/2], [month_stats['lower_whisker'], month_stats['lower_whisker']], color='k')
#         ax.plot([positions[i] - widths/2, positions[i] + widths/2], [month_stats['upper_whisker'], month_stats['upper_whisker']], color='k')
#         ax.add_patch(plt.Rectangle((positions[i] - widths/2, month_stats['q1']), widths, month_stats['q3'] - month_stats['q1'], fill=True, color='lightblue'))
#         ax.plot([positions[i] - widths/2, positions[i] + widths/2], [month_stats['median'], month_stats['median']], color='k')
#
#     ax.set_xticks(positions)
#     ax.set_xticklabels(positions)
#     ax.set_title(title)
#     ax.set_xlabel('Month')
#     ax.set_ylabel('Anomaly (relative to standard deviation)')
#     ax.grid(True)
#
#     plt.savefig(save_path)
#     plt.close()

def get_sscast_reg_idx(file_path):
    if os.path.exists(file_path):
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                metadata = {
                    "Dimensions": {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions},
                    "Variables": {var: dataset.variables[var].dimensions for var in dataset.variables}
                }
                lats = dataset.variables['latitude'][:]
                lons = dataset.variables['longitude'][:]
                gridY = dataset.variables['gridY'][:]
                gridX = dataset.variables['gridX'][:]
                bathy = dataset.variables['bathymetry'][:]

            return metadata, lats, lons, gridY, gridX, bathy
        except OSError as e:
            return {"Error": str(e)}, None
    else:
        return {"Error": "File not found."}, None

def get_ecospace_map(file_path):
    if os.path.exists(file_path):
        df = pd.read_csv(file_path)
        metadata = {
            "Column Names": df.columns.tolist(),
            "Number of Rows": len(df)
        }
        return metadata, df
    else:
        return {"Error": "File not found."}, None

def get_ecospace_data(file_path):
    if os.path.exists(file_path):
        try:
            with nc.Dataset(file_path, 'r') as dataset:
                metadata = {
                    "Dimensions": {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions},
                    "Variables": {var: dataset.variables[var].dimensions for var in dataset.variables}
                }
                time_var = dataset.variables['time']
                time_units = time_var.units
                time_calendar = time_var.calendar if 'calendar' in time_var.ncattrs() else 'standard'
                times = nc.num2date(time_var[:], units=time_units, calendar=time_calendar)

                PZ1_CIL = dataset.variables['PZ1-CIL'][:]
                PP1_DIA = dataset.variables['PP1-DIA'][:]
                PP2_NAN = dataset.variables['PP2-NAN'][:]
                PP3_PIC = dataset.variables['PP3-PIC'][:]

                EWE_row = dataset.variables['EWE_row'][:]
                EWE_col = dataset.variables['EWE_col'][:]

            return metadata, times, PZ1_CIL, PP1_DIA, PP2_NAN, PP3_PIC, EWE_row, EWE_col
        except OSError as e:
            return {"Error": str(e)}, None, None, None, None
    else:
        return {"Error": "File not found."}, None, None, None, None

def convert_times(times, time_units, time_calendar):
    times = nc.num2date(times, units=time_units, calendar=time_calendar)
    times = np.array([np.datetime64(t) for t in times])
    return times

def reshape_to_year_month(data, years, months):
    unique_years = np.unique(years)
    unique_months = np.unique(months)

    year_dim = len(unique_years)
    month_dim = len(unique_months)

    reshaped_data = np.empty((year_dim, month_dim, data.shape[1], data.shape[2]))

    for i, year in enumerate(unique_years):
        for j, month in enumerate(unique_months):
            mask = (years == year) & (months == month)
            reshaped_data[i, j, :, :] = data.sel(time=mask).mean(dim='time')

    return reshaped_data

def get_ecospace_data_2(file_path, start_year=None, end_year=None, monthly_averages=False):
    if os.path.exists(file_path):
        try:
            with (nc.Dataset(file_path, 'r') as dataset):
                metadata = {
                    "Dimensions": {dim: len(dataset.dimensions[dim]) for dim in dataset.dimensions},
                    "Variables": {var: dataset.variables[var].dimensions for var in dataset.variables}
                }
                time_var = dataset.variables['time']
                time_units = time_var.units
                time_calendar = time_var.calendar if 'calendar' in time_var.ncattrs() else 'standard'
                times = convert_times(time_var[:], time_units, time_calendar)

                if start_year is not None and end_year is not None:
                    start_date = np.datetime64(f'{start_year}-01-01')
                    end_date = np.datetime64(f'{end_year}-12-31')
                    time_mask = (times >= start_date) & (times <= end_date)
                    times = times[time_mask]
                else:
                    time_mask = slice(None)

                PZ1_CIL = dataset.variables['PZ1-CIL'][time_mask]
                PP1_DIA = dataset.variables['PP1-DIA'][time_mask]
                PP2_NAN = dataset.variables['PP2-NAN'][time_mask]
                PP3_PIC = dataset.variables['PP3-PIC'][time_mask]

                EWE_row = dataset.variables['row'][:]
                EWE_col = dataset.variables['col'][:]
                EWE_depth = dataset.variables['depth'][:]
                row_indices = EWE_row
                col_indices = EWE_col

                ds = xr.Dataset({
                    'PZ1_CIL': (['time', 'EWE_row', 'EWE_col'], PZ1_CIL),
                    'PP1_DIA': (['time', 'EWE_row', 'EWE_col'], PP1_DIA),
                    'PP2_NAN': (['time', 'EWE_row', 'EWE_col'], PP2_NAN),
                    'PP3_PIC': (['time', 'EWE_row', 'EWE_col'], PP3_PIC)
                }, coords={'time': times, 'EWE_row': row_indices, 'EWE_col': col_indices})

                if monthly_averages:
                    ds = ds.resample(time='1M').mean()

                reshaped_data = {}
                for var in ds.data_vars:
                    data = ds[var]
                    years = data['time.year']
                    months = data['time.month']
                    reshaped_data[var] = reshape_to_year_month(data, years, months)

                return metadata, reshaped_data, row_indices, col_indices, EWE_depth
        except OSError as e:
            return {"Error": str(e)}, None, None, None, None
    else:
        return {"Error": "File not found."}, None, None, None, None

def prepare_data_for_plotting(subdomain_array_ssc, subdomain_array_ecospace, ssc_data, ecospace_data, subdomains, tmask_ssc, tmask_ecospace):
    results = {}

    for subdomain in subdomains:
        results[subdomain] = {
            'ecospace': {},
            'ssc': {}
        }

        for data_name, data_values in ssc_data.items():
            anomalies = compute_subdomain_anomalies(data_values, subdomain_array_ssc, subdomain)
            avg_anomalies = np.nanmean(anomalies, axis=0)
            stats = get_monthly_stats(avg_anomalies, tmask_ssc)
            results[subdomain]['ssc'][data_name] = stats

        for data_name, data_values in ecospace_data.items():
            anomalies = compute_subdomain_anomalies(data_values, subdomain_array_ecospace, subdomain)
            avg_anomalies = np.nanmean(anomalies, axis=0)
            stats = get_monthly_stats(avg_anomalies, tmask_ecospace)
            results[subdomain]['ecospace'][data_name] = stats

    return results

# def prepare_data_for_plotting(subdomain_array_ssc, subdomain_array_ecospace, ssc_data, ecospace_data, subdomains, tmask_ssc, tmask_ecospace):
#     results = {}
#
#     for subdomain in subdomains:
#         results[subdomain] = {
#             'ecospace': {},
#             'ssc': {}
#         }
#
#         for data_name, data_values in ssc_data.items():
#             anomalies = compute_subdomain_anomalies(data_values, subdomain_array_ssc, subdomain)
#             avg_anomalies = np.nanmean(anomalies, axis=0)
#             stats = get_monthly_stats(avg_anomalies, tmask_ssc)
#             results[subdomain]['ssc'][data_name] = stats
#
#         for data_name, data_values in ecospace_data.items():
#             anomalies = compute_subdomain_anomalies(data_values, subdomain_array_ecospace, subdomain)
#             avg_anomalies = np.nanmean(anomalies, axis=0)
#             stats = get_monthly_stats(avg_anomalies, tmask_ecospace)
#             results[subdomain]['ecospace'][data_name] = stats
#
#     return results

# def prepare_data_for_plotting(subdomain_array_ssc, subdomain_array_ecospace, ssc_data, ecospace_data, subdomains, tmask_ssc, tmask_ecospace):
#     results = {}
#
#     for subdomain in subdomains:
#         results[subdomain] = {}
#         for data_name, data_values in ssc_data.items():
#             anomalies = compute_subdomain_anomalies(data_values, subdomain_array_ssc, subdomain)
#             avg_anomalies = np.nanmean(anomalies, axis=0)
#             stats = get_monthly_stats(avg_anomalies, tmask_ssc)
#             results[subdomain][data_name] = stats
#
#         for data_name, data_values in ecospace_data.items():
#             anomalies = compute_subdomain_anomalies(data_values, subdomain_array_ecospace, subdomain)
#             avg_anomalies = np.nanmean(anomalies, axis=0)
#             stats = get_monthly_stats(avg_anomalies, tmask_ecospace)
#             results[subdomain][data_name] = stats
#
#     return results

def plot_results(results, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for subdomain, data in results.items():
        for data_name, stats in data.items():
            plot_custom_boxplot(stats, f'{data_name} in {subdomain}', os.path.join(output_dir, f'{data_name}_anomalies_{subdomain}_boxplot.png'))

def main():
    # Load SSC data
    ssc_data = {
        "diatoms": [],
        "flagellates": [],
        "ciliates": []
    }

    for year in range(yr_st, yr_en + 1):
        dep_avg_dia, dep_avg_fla, dep_avg_cil = process_year(year)
        if dep_avg_dia is not None:
            ssc_data["diatoms"].append(dep_avg_dia)
            ssc_data["flagellates"].append(dep_avg_fla)
            ssc_data["ciliates"].append(dep_avg_cil)

    for key in ssc_data:
        ssc_data[key] = np.stack(ssc_data[key], axis=0)

    # Load Ecospace data
    ecospace_md, ecospace_data, row_indices, col_indices, EWE_depth = get_ecospace_data_2(
        os.path.join(ecospace_p, ecospace_f),
        start_year=2008,
        end_year=2018,
        monthly_averages=True
    )

    # Prepare data for plotting
    with nc.Dataset(os.path.join(ssc_p, meshm2D_f)) as mesh:
        tmask_ssc = mesh.variables['tmaskutil'][:]  # dry land mask (1, 898, 398), 0's 1's

    metadata_ssc, lats_ssc, lons_ssc, gridY_ssc, gridX_ssc, bathy_ssc = get_sscast_reg_idx(os.path.join(ssc_p, bathy_f))

    subdomain_paths = {name: Path(coords) for name, coords in sdomains.items()}
    subdomain_array_ssc = np.full((lats_ssc.shape[0], lats_ssc.shape[1]), '', dtype=object)

    for i in range(lats_ssc.shape[0]):
        for j in range(lats_ssc.shape[1]):
            if tmask_ssc[0, i, j] == 1:
                point = (lats_ssc[i, j], lons_ssc[i, j], )
                for name, path in subdomain_paths.items():
                    if path.contains_point(point):
                        subdomain_array_ssc[i, j] = name
                        break

    # Prepare Ecospace subdomain array and mask
    ecospace_coords_md, ecospace_coords_df = get_ecospace_map(os.path.join(ecospace_map_p, ecospace_map_f))
    ecospace_coords_df['sdomain'] = ''
    for index, row in ecospace_coords_df.iterrows():
        pos = (row['lat'], row['lon'])
        sd_es = ""
        for sd in sdomains.keys():
            if Path(sdomains[sd]).contains_point(pos):
                sd_es = sd
                break
        ecospace_coords_df.at[index, 'sdomain'] = sd_es

    subdomain_array_ecospace = np.full((151, 93), '', dtype=object)
    for _, row in ecospace_coords_df.iterrows():
        ewe_row = int(row['EWE_row'])
        ewe_col = int(row['EWE_col'])
        sdomain = row['sdomain']
        subdomain_array_ecospace[ewe_row - 1, ewe_col - 1] = sdomain

    EWE_landmask = np.where(EWE_depth > 0, 1, 0)

    results = prepare_data_for_plotting(subdomain_array_ssc, subdomain_array_ecospace, ssc_data, ecospace_data,
                                        ['SGN', 'SGS', 'SGI'], tmask_ssc, EWE_landmask)

    for subdomain, data in results.items():
        plot_custom_boxplot(data['ecospace']['PP1_DIA'], data['ssc']['diatoms'], f'Diatoms in {subdomain}',
                            f'../figs/diatoms_ecospace_vs_ssc_anomalies_{subdomain}_boxplot.png')


if __name__ == "__main__":
    main()