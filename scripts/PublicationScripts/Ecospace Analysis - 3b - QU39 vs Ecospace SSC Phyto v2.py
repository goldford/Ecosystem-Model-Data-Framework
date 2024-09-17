# CHAT GPT ATTEMPT to revise to v2 - not working


import pandas as pd
import calendar
import numpy as np
import matplotlib.pyplot as plt

# Paths
# paths
model_run = 'SC39' # older run
model_run = 'SC51_4_2' # this is short run but same as full key run
include_SSC = True
pathfile_QU39_SSC_Ecospace = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED\QU39_joined_matchtoEcospace_SSC_' + model_run + '.csv'
out_p = "..//..//figs//"


QU39_df = pd.read_csv(pathfile_QU39_SSC_Ecospace, sep=',', engine='python')
QU39_df['DateTime'] = pd.to_datetime(QU39_df['DateTime'], format='%Y-%m-%d %H:%M:%S')

# Define the columns to group by and the columns to sum
groupby_columns = ['eventID', 'year', 'month','class']
sum_columns = ['PP1-DIA', 'PP2-NAN', 'PP3-PIC', 'PZ1-CIL', 'PZ2-DIN',
               'ssc-DIA', 'ssc-FLA', 'ssc-CIL',
               'measurementValue'
               ]

# Group by the specified columns and compute the sum for the specified fields
QU39_df = QU39_df.groupby(groupby_columns)[sum_columns].sum().reset_index()

QU39_df['QU39'] = QU39_df['measurementValue']
QU39_df['logQU39'] = np.log(QU39_df['measurementValue'])

QU39_df['log_PP1-DIA'] = np.log(QU39_df['PP1-DIA'])
QU39_df['log_PP2-NAN'] = np.log(QU39_df['PP2-NAN'])
QU39_df['log_PP3-PIC'] = np.log(QU39_df['PP3-PIC'])
QU39_df['log_PZ2-DIN'] = np.log(QU39_df['PZ2-DIN']) # added 2024-09

QU39_df['log_ssc-DIA'] = np.log(QU39_df['ssc-DIA'])
QU39_df['log_ssc-FLA'] = np.log(QU39_df['ssc-FLA'])
QU39_df['log_ssc-CIL'] = np.log(QU39_df['ssc-CIL'])

# Function to calculate anomalies
def calculate_anomalies(df, valuefield, groupon):
    df[valuefield + '_anomaly_fr_mean'] = -999
    df[valuefield + '_anomaly_fr_median'] = -999
    df[valuefield + '_anomaly_fr_mean_std_norm'] = -999
    df[valuefield + '_anomaly_fr_median_std_norm'] = -999
    df[valuefield + '_log_anomaly_fr_mean_std_norm'] = -999

    df_filt = df[df[valuefield] > -999]
    monthly_year_means = df_filt.groupby(['year', 'month', groupon])[valuefield].mean().reset_index()
    monthly_year_medians = df_filt.groupby(['year', 'month', groupon])[valuefield].median().reset_index()

    # Compute annual mean and std
    monthly_means = monthly_year_means.groupby([groupon, 'month'])[valuefield].mean()
    annual_means = monthly_means.groupby([groupon]).mean()
    annual_std = monthly_means.groupby([groupon]).std()

    df = df.merge(annual_means, on=[groupon], suffixes=('', '_annual_mean'))
    df = df.merge(annual_std, on=[groupon], suffixes=('', '_annual_std'))

    # Calculate the anomaly
    df[valuefield + '_anomaly_fr_mean'] = -999
    df[valuefield + '_anomaly_fr_mean_std_norm'] = -999
    condition = df[valuefield] > -999

    df.loc[condition, valuefield + '_anomaly_fr_mean'] = (df[valuefield] - df[valuefield + '_annual_mean'])
    df.loc[condition, valuefield + '_anomaly_fr_mean_std_norm'] = (
                df[valuefield + '_anomaly_fr_mean'] / df[valuefield + '_annual_std'])

    return df


# Function to plot boxplots
def plot_boxplots(data, field_codes, colors, source_names, positions, months, widths, output_path, plot_label,
                  fig_label, include_SSC=False):
    plt.figure(figsize=(15, 10) if include_SSC else (8, 4))
    fig, ax = plt.subplots()

    for month in months:
        # Plotting data for each functional group
        for idx, field_code in enumerate(field_codes):
            model_data = data[data['month'] == month][field_code + '_anomaly_fr_mean_std_norm'].dropna().values
            plt.boxplot(model_data, positions=[positions[idx][month - 1]], widths=widths, patch_artist=True,
                        boxprops=dict(facecolor=colors[field_code], color='black'),
                        medianprops=dict(color='black'), whiskerprops=dict(color='black'),
                        capprops=dict(color='black'), showfliers=False)

    plt.xlabel('Month')
    plt.ylabel('Anomaly')
    plt.xticks(months, [calendar.month_abbr[m] for m in months])

    handles = [plt.Line2D([0], [0], color=colors[code], lw=4) for code in field_codes]
    plt.legend(handles, source_names, loc='upper right')
    plt.text(0.02, 0.95, fig_label, transform=ax.transAxes, fontsize=14, verticalalignment='top')

    plt.tight_layout()
    plt.savefig(output_path + plot_label + ".png")
    plt.show()


# Configuration of functional groups
functional_groups = {
    "PP1-DIA": {"field_codes": ['log_PP1-DIA', 'logQU39', 'log_ssc-DIA'], "taxon": ['Bacillariophyceae'],
                "plot_label": "Diatoms", "fig_label": "(a)"},
    "PP2-NAN": {"field_codes": ['log_PP2-NAN', 'logQU39', 'log_ssc-FLA'], "taxon": ['Choanoflagellatea', 'Dictyochophyceae',
                                                                                    'Cryptophyceae', 'Metromonadea',
                                                                                    'Chrysophyceae', 'Telonemea',
                                                                                    'Chlorodendrophyceae', 'Bicosoecophyceae',
                                                                                    'Xanthophyceae', 'Coccolithophyceae',
                                                                                    'Euglenophyceae', 'Raphidophyceae'],
                "plot_label": "Nanophytoplankton", "fig_label": "(b)"},
    "PP3-PIC": {"field_codes": ['log_PP3-PIC', 'logQU39'], "taxon": ['Pyramimonadophyceae'], "plot_label": "Picoplankton",
                "fig_label": "(c)"}
}

# Colors configuration
colors = {'log_PP1-DIA': 'blue', 'log_ssc-DIA': 'orange', 'logQU39': 'pink',
          'log_PP2-NAN': 'blue', 'log_ssc-FLA': 'pink',
          'log_PP3-PIC': 'blue'}

# Boxplot positions
months = np.arange(1, 13)
positions_ecospace = months - 0.3
positions_ssc = months
positions_qu39 = months + 0.3
positions = [positions_ecospace, positions_ssc, positions_qu39]
widths = 0.25

# Process and plot each functional group
for group, config in functional_groups.items():
    # Calculate anomalies
    for valuefield in config['field_codes']:
        QU39_df_anom = calculate_anomalies(QU39_df, valuefield, 'class')

    # Filter data
    class_data = QU39_df_anom[QU39_df_anom['class'].isin(config['taxon'])]

    print(class_data.head)
    print(class_data[config['field_codes'][0] + '_anomaly_fr_mean'])

    # why is this empty?
    condition = class_data[config['field_codes'][0] + '_anomaly_fr_mean'] > -999
    class_data = class_data.loc[condition]
    print("three")
    print(class_data)

    # class_data = QU39_df_anom[QU39_df_anom['class'] == config['taxon']]

    # print(group)

    # Plot boxplots
    plot_boxplots(class_data, config['field_codes'], colors, ['Ecospace', 'SalishSeaCast', 'QU39'], positions, months,
                  widths, out_p, config['plot_label'], config['fig_label'], include_SSC)
