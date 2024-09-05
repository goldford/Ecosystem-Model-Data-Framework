# Created: Sep 1, 2024
# Author: G Oldford
# Purpose: Plot Ecospace Community Composition for North and South SoG and
#          compare against community composition from literature
#
# Input files:
#  - Ecospace output NC file (e.g., Scv51-RSPI_AllTemp_Wind_2000-2018.nc)
#  -
# Output files / figs:
#
# Sources of Data:
#  -

# File name examples:
#  - SScast: SalishSeaCast_biology_2008.nc
#  -
#
# Notes:
#   -

# compute Ecospace seasonal average biomasses grouped by north and south SoG
# separately for phyto and zoop
# compare against the estimates of McEwan

from helpers import read_sdomains, get_ecospace_data_3
import os
import netCDF4 as nc
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt
import pandas as pd

# Paths
# ecospace
ecospace_p = "C:/Users/Greig/Sync/PSF/EwE/Georgia Strait 2021/LTL_model/ECOSPACE_OUT"
# file_Ecospace = "Scv7-PARMixingNut90Temp_2003-2018.nc"
# ecospace_f = "Scv51-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC51'
# ecospace_f = "Scv53-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC53' # same as 52 but with temp PROD response applied to DIA (proxy for nutrients)
# ecospace_f = "Scv54-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC54' # same as 53 but with temp MORT response applied to DIA (proxy for nutrients)
# ecospace_f = "Scv56-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC56' # same as 54 and 55 but with temp MORT response applied to DIA (proxy for nutrients)
#                        # AND hab cap driver (PAR) - this seems to make variability higher
# ecospace_f = "Scv58-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC58' # same as 54 and 55 but with temp MORT response applied to DIA (proxy for nutrients)
#                        # AND hab cap driver (PAR) - this seems to make variability higher
#                        # with temp penalty response driver for DIA removed to compared variability to 56
# ecospace_f = "Scv59-RSPI_AllTemp_Wind_2001-2018.nc"
# ecospace_code = 'SC59' # same as 58 with temp MORT response applied to DIA (proxy for nutrients)
#                        # NO hab capacity driver data connected and with SPIN UP data properly repeated
# ecospace_f = "Scv64-RSPI_AllTemp_Wind_2001-2018.nc"
# ecospace_code = 'SC64' #
# ecospace_f = "Scv56_2-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC56_2' #
# ecospace_f = "Scv56_3-RSPI_AllTemp_Wind_2000-2018.nc"
# ecospace_code = 'SC56_3' #
ecospace_f = "Scv56_5-RSPI_AllTemp_Wind_2000-2018.nc"
ecospace_code = 'SC56_5' #


# path_Nemcek = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/Phytoplankton Salish Sea Nemcek2023 2015-2019/MODIFIED"

# SSC
path_SSC = "C:/Users/Greig/Sync/6. SSMSP Model/Model Greig/Data/28. Phytoplankton/SalishSeaCast_BioFields_2008-2018/ORIGINAL"
file_SSC_mo = "SalishSeaCast_biology_2008.nc"
# https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSnBathymetryV21-08.nc?bathymetry[(0.0):1:(897.0)][(0.0):1:(397.0)],latitude[(0.0):1:(897.0)][(0.0):1:(397.0)],longitude[(0.0):1:(897.0)][(0.0):1:(397.0)]
file_SSC_grd = "ubcSSnBathymetryV21-08_a29d_efc9_4047.nc"

# subregions for aggregated analysis
subregions_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
mcewan_subregions_f = "analysis_domains_mcewan.yml"

# Ecospace map
ecospace_map_p = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/basemap"
ecospace_map_f = "Ecospace_grid_20210208_rowscols.csv"

# output
path_evalout = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"

yr_st = 2008
yr_en = 2018


def main():
    # Load Ecospace data as monthly averages
    ecospace_md, ecospace_data, row_indices, col_indices, EWE_depth = get_ecospace_data_3(
        os.path.join(ecospace_p, ecospace_f),
        start_year=yr_st,
        end_year=yr_en,
        monthly_averages=True
    )

    sdomains = read_sdomains(os.path.join(subregions_p, mcewan_subregions_f))
    subdomain_paths = {name: Path(coords) for name, coords in sdomains.items()}

    # Get Ecospace map
    if os.path.exists(os.path.join(ecospace_map_p, ecospace_map_f)):
        ecospace_coords_df = pd.read_csv(os.path.join(ecospace_map_p, ecospace_map_f))
        ecospace_coords_md = {
            "Column Names": ecospace_coords_df.columns.tolist(),
            "Number of Rows": len(ecospace_coords_df)
        }
    else:
        print("Error: Ecospace map not found.")

    # assign domains
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

    results = {}

    for subdomain in sdomains:
        if subdomain not in results:
            results[subdomain] = {}

        for data_name, data_values in ecospace_data.items():

            if data_name not in results[subdomain]:
                results[subdomain][data_name] = {}

            # mask = (subdomain_array_ecospace != '') & np.isfinite(data_values)
            # valid_data = data_values[mask]

            mask = (subdomain_array_ecospace == subdomain) & np.isfinite(data_values)
            valid_sd_data = np.where(mask, data_values, np.nan)

            seasons = ['winter', 'spring', 'summer']
            for season in seasons:

                if season not in results[subdomain][data_name]:
                    results[subdomain][data_name][season] = {}

                  # the minus ones make it easier to cross check
                # spring(Mar - May)
                # summer(Jun - Oct)
                # winter(Nov - Feb)
                if season == 'winter':
                    season_data = valid_sd_data[:, [11 - 1, 12 - 1, 1 - 1, 2 - 1], :, :]
                elif season == 'spring':
                    season_data = valid_sd_data[:, [3 - 1, 4 - 1, 5 - 1], :, :]
                elif season == 'summer':
                    season_data = valid_sd_data[:, [6 - 1, 7 - 1, 8 - 1, 9 - 1, 10 - 1], :, :]

                season_data = np.where(EWE_landmask == 0, np.nan, season_data)
                season_data = season_data.flatten()

                # compute stats on the seasonal data for this subregion
                stats = {}
                season_data = season_data[season_data > -999]
                if len(season_data) == 0:
                    return None  # Return None if no valid data

                season_data = np.ma.masked_invalid(season_data)

                # enough stats for boxplots
                q1 = np.percentile(season_data, 25, axis=0)
                q3 = np.percentile(season_data, 75, axis=0)
                median = np.percentile(season_data, 50, axis=0)
                mean = np.nanmean(season_data, axis=0)
                std_dev = np.nanstd(season_data, axis=0)
                iqr = q3 - q1
                lower_whisker = q1 - 1.5 * iqr
                upper_whisker = q3 + 1.5 * iqr
                lower_whisker = np.maximum(lower_whisker, np.min(season_data, axis=0))
                upper_whisker = np.minimum(upper_whisker, np.max(season_data, axis=0))
                stats['q1'] = q1
                stats['q3'] = q3
                stats['median'] = median
                stats['mean'] = mean
                stats['std'] = std_dev
                stats['lower_whisker'] = lower_whisker
                stats['upper_whisker'] = upper_whisker

                results[subdomain][data_name][season] = stats

    print(results)

    # Bar plot
    bar_width = 0.35
    index = np.arange(len(seasons))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Function to plot a stacked bar chart on a given axis
    def plot_stacked_bars(ax, region, title):
        print('region')
        print(region)
        PP_PZ_keys = [key for key in results[region] if key.startswith('PP') or key.startswith('PZ')]

        # Stacking the bars
        bottom = np.zeros(len(seasons))
        for key in PP_PZ_keys:
            print(key)

            # Adjust for proportion mixotrophs
            if key == 'PZ1_CIL':
                PP_scale = 0.125
            elif key == 'PZ2_DIN':
                PP_scale = 0.125
            elif key == 'PZ3_HNF':
                continue
            else:
                PP_scale = 1

            values = []
            stds = []
            for season in seasons:
                mean = results[region][key][season]['mean']
                std = results[region][key][season]['std']  # Assuming your data structure has 'std'
                mean_adjusted = PP_scale * mean
                std_adjusted = PP_scale * std # wrong, to do - fix
                values.append(mean_adjusted)
                stds.append(std_adjusted)

            print(values)
            ax.bar(index, values, bar_width, bottom=bottom, yerr=stds, label=key)
            bottom += np.array(values)

        ax.set_xlabel('Season')
        ax.set_ylabel('Mean Value')
        ax.set_title(title)
        ax.set_xticks(index)
        ax.set_xticklabels(seasons)
        ax.legend()

    # Plot for 'NSoG'
    plot_stacked_bars(ax1, 'NSoG', 'Panel (a): NSoG')

    # Plot for 'SSoG'
    plot_stacked_bars(ax2, 'SSoG', 'Panel (b): SSoG')

    plt.tight_layout()

    # Save the plot to a file
    plt.savefig('..//..//figs//seasonal_B_plots_NS_S_' + ecospace_code + '.png', bbox_inches='tight')

    plt.show()

if __name__ == "__main__":
    main()