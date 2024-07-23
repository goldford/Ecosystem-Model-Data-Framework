# Created by G Oldford
# July 19 2024
# Purpose: use previously prepped data from del bel belluz at stn QU39 (2016 - 2021)
#           matched in previous script to Ecospace and SSCast outputs (SMELT)
# Source:
#    Del Bel Belluz, J. (2024). Protistan plankton time series from the northern Salish Sea and central coast,
#    British Columbia, Canada—Ocean Biodiversity Information System (a62c37c4-6fbe-4a93-8c87-c67dda36436c)
#    [Dataset]. https://obis.org/dataset/071a38b3-0d30-47ba-85aa-d98313460075
# Input:
#    file (csv) with bottle sampled ID'd phytoplankton (microscopy) matched to Ecospace and SalishSeaCast (SMELT) out
#
# Output:
#    plots of data desriptors and distributions, monthly boxplots of model vs obs phytoplankton
#
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# paths
pathfile_QU39_SSC_Ecospace = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED\QU39_joined_matchtoEcospaceAndSSC.csv'

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
QU39_df['log_ssc-DIA'] = np.log(QU39_df['ssc-DIA'])
QU39_df['log_ssc-FLA'] = np.log(QU39_df['ssc-FLA'])
QU39_df['log_ssc-CIL'] = np.log(QU39_df['ssc-CIL'])

QU39_df['ecospace_PP1-DIA_anomaly_fr_mean'] = 0.0
QU39_df['ecospace_PP1-DIA_anomaly_fr_median'] = 0.0
QU39_df['ecospace_PP1-DIA_annual_mean'] = 0.0
QU39_df['ecospace_PP1-DIA_annual_median'] = 0.0
QU39_df['ecospace_PP1-DIA_annual_std'] = 0.0

# QU39_df['SSC_ssc-DIA_anomaly_fr_mean'] = 0.0
# QU39_df['SSC_ssc-DIA_anomaly_fr_median'] = 0.0
# QU39_df['SSC_ssc-DIA_annual_mean'] = 0.0
# QU39_df['SSC_ssc-DIA_annual_median'] = 0.0
# QU39_df['SSC_ssc-DIA_annual_std'] = 0.0


# anomalies were already calculated for the measurements in a previous step
# do it here for ecospace, ssc
# prefixes: ecospace, ssc; valuefield: e.g., PP1-DIA, ssc-dia (original was measurementValue
def calculate_anomalies(df, valuefield, groupon):

    # Compute monthly means for each year and class
    df_filt = df[df[valuefield] > -999]

    monthly_year_means = df_filt.groupby(['year', 'month', groupon])[valuefield].mean().reset_index()
    monthly_year_medians = df_filt.groupby(['year', 'month', groupon])[valuefield].median().reset_index()

    # Compute annual mean and standard deviation from the monthly means for each class
    monthly_means = monthly_year_means.groupby([groupon, 'month'])[valuefield].mean()
    monthly_medians = monthly_year_medians.groupby([groupon, 'month'])[valuefield].median()
    annual_means = monthly_means.groupby([groupon]).mean()
    annual_std = monthly_means.groupby([groupon]).std()
    annual_medians = monthly_medians.groupby([groupon]).median()

    # Merge the annual means and standard deviations back into the original dataframe
    df = df.merge(annual_means, on=[groupon], suffixes=('', '_annual_mean'))
    df = df.merge(annual_std, on=[groupon], suffixes=('', '_annual_std'))
    df = df.merge(annual_medians, on=[groupon], suffixes=('', '_annual_median'))

    # Calculate the anomaly
    df[valuefield + '_anomaly_fr_mean'] = -999
    df[valuefield + '_anomaly_fr_median'] = -999
    df[valuefield + '_anomaly_fr_mean_std_norm'] = -999
    df[valuefield + '_anomaly_fr_median_std_norm'] = -999
    df[valuefield + '_log_anomaly_fr_mean_std_norm'] = -999

    condition = df[valuefield] > -999

    df.loc[condition, valuefield + '_anomaly_fr_mean'] = (df[valuefield] - df[valuefield + '_annual_mean'])
    df.loc[condition, valuefield + '_anomaly_fr_median'] = (df[valuefield] - df[valuefield + '_annual_median'])
    df.loc[condition, valuefield + '_anomaly_fr_mean_std_norm'] = (df[valuefield] - df[valuefield + '_annual_mean']) / df[valuefield + '_annual_std']
    df.loc[condition, valuefield + '_anomaly_fr_median_std_norm'] = (df[valuefield] - df[valuefield + '_annual_median']) / df[
        valuefield + '_annual_std']
    df.loc[condition, valuefield + '_log_anomaly_fr_mean_std_norm'] = np.log(df[valuefield + '_anomaly_fr_mean_std_norm'] + 1 - df[valuefield + '_anomaly_fr_mean_std_norm'].min())

    # df[df[valuefield] == -999] = -999

    return df

# anomalies for diatoms for ssc and ecospace
groupon = 'class'
# groupon = 'family'

# TO DO - SUM FIRST  - grouped BY SAMPLING EVENT and 'groupon'
# print(QU39_df.columns)
# df = QU39_df.groupby(['eventID', 'year','month', groupon])[valuefield].sum().reset_index()


valuefield = 'ssc-DIA'
valuefield = 'log_ssc-DIA'
print(QU39_df.columns)
QU39_df = calculate_anomalies(QU39_df, valuefield, groupon)

valuefield = 'log_PP1-DIA'
QU39_df = calculate_anomalies(QU39_df, valuefield, groupon)

# I calculated anomalies previously but am doing it again here to get it right
valuefield = 'QU39'
valuefield = 'logQU39'
QU39_df = calculate_anomalies(QU39_df, valuefield, groupon)


# Configuration
taxon = 'Bacillariophyceae'
# taxon = 'Chaetocerotaceae'
# taxon = 'Thalassiosiraceae'

# field_codes = ['PP1-DIA', 'ssc-DIA', 'QU39']
field_codes = ['log_PP1-DIA', 'log_ssc-DIA', 'logQU39']
# source_names = ['Ecospace', 'SalishSeaCast', 'QU39']
source_names = ['logEcospace', 'logSalishSeaCast', 'logQU39']

# Colors for each model_code
# colors = {'PP1-DIA': 'blue', 'ssc-DIA': 'orange', 'logQU39': 'pink'}
colors = {'log_PP1-DIA': 'blue', 'log_ssc-DIA': 'orange', 'logQU39': 'pink'}

months = np.arange(1, 13)
positions_ecospace = months - 0.3 # adjust slightly right
positions_ssc = months
positions_qu39 = months + 0.3  # Adjust positions slightly to the left
widths = 0.25

plt.figure(figsize=(15, 10))

# Data filtering
class_data = QU39_df[QU39_df[groupon] == taxon]
condition = class_data[field_codes[0] + '_anomaly_fr_mean'] > -999
class_data = class_data.loc[condition]

print(len(class_data['eventID'].unique()))

# Create a dictionary to hold data by month for each model
data_by_month = {month: {} for month in months}

# Extract data for each model and each month
for field_codes in field_codes:
    for month in months:
        month_data = class_data[class_data['month'] == month]
        data_by_month[month][field_codes] = month_data[field_codes + '_anomaly_fr_mean_std_norm'].dropna().values
        # data_by_month[month][field_codes] = month_data[field_codes + '_log_anomaly_fr_mean_std_norm'].dropna().values

# Plot boxplots for both models side-by-side for each month
for month in months:
    print("month: " + str(month))

    # Data for both models for the current month
    ecospace_data = data_by_month[month].get('log_PP1-DIA', [])
    qu39_data = data_by_month[month].get('logQU39', [])
    ssc_data = data_by_month[month].get('log_ssc-DIA', [])

    # Plot the first model\
    plt.boxplot(ecospace_data, positions=[positions_ecospace[month - 1]], widths=widths, patch_artist=True,
                boxprops=dict(facecolor=colors['log_PP1-DIA'], color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                showfliers=False)

    # # Plot the second model
    plt.boxplot(ssc_data, positions=[positions_ssc[month - 1]], widths=widths, patch_artist=True,
                boxprops=dict(facecolor=colors['log_ssc-DIA'], color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                showfliers=False)

    plt.boxplot(qu39_data, positions=[positions_qu39[month - 1]], widths=widths, patch_artist=True,
                boxprops=dict(facecolor=colors['logQU39'], color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                showfliers=False)

plt.title(f'Anomaly from Means Rel to Std Dev by Month for {taxon}')
plt.suptitle('')  # Suppress the default title to keep it clean
plt.xlabel('Month')
plt.ylabel('Mean')

# Adjust x-axis ticks and labels
plt.xticks(months, [f'Month {m}' for m in months])

# Create custom legend
# leg_codes = ['PP1-DIA', 'ssc-DIA', 'logQU39']
leg_codes = ['log_PP1-DIA', 'log_ssc-DIA', 'logQU39']
print(leg_codes)

for leg_code in leg_codes:
    print(leg_code)
    print(colors[leg_code])
handles = [plt.Line2D([0], [0], color=colors[leg_code], lw=4) for leg_code in leg_codes]
plt.legend(handles, source_names, loc='upper right')

plt.tight_layout()
plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_vs_Ecospace_vs_SSC_2016_2018.png')
plt.show()

##########################

# Configuration
taxon = 'Bacillariophyceae'
# taxon = 'Chaetocerotaceae'
# taxon = 'Thalassiosiraceae'

# field_codes = ['PP1-DIA', 'ssc-DIA', 'QU39']
field_codes = ['log_PP1-DIA', 'logQU39']
# source_names = ['Ecospace', 'SalishSeaCast', 'QU39']
source_names = ['logEcospace', 'logQU39']

# Colors for each model_code
# colors = {'PP1-DIA': 'blue', 'ssc-DIA': 'orange', 'logQU39': 'pink'}
colors = {'log_PP1-DIA': 'blue', 'logQU39': 'orange'}

months = np.arange(1, 13)
positions_ecospace = months - 0.2 # adjust slightly right
# positions_ssc = months
positions_qu39 = months + 0.2  # Adjust positions slightly to the left
widths = 0.35

plt.figure(figsize=(15, 10))

# Data filtering
class_data = QU39_df[QU39_df[groupon] == taxon]
condition = class_data[field_codes[0] + '_anomaly_fr_mean'] > -999
class_data = class_data.loc[condition]

print(len(class_data['eventID'].unique()))

# Create a dictionary to hold data by month for each model
data_by_month = {month: {} for month in months}

# Extract data for each model and each month
for field_codes in field_codes:
    for month in months:
        month_data = class_data[class_data['month'] == month]
        data_by_month[month][field_codes] = month_data[field_codes + '_anomaly_fr_mean_std_norm'].dropna().values
        # data_by_month[month][field_codes] = month_data[field_codes + '_log_anomaly_fr_mean_std_norm'].dropna().values

# Plot boxplots for both models side-by-side for each month
for month in months:
    print("month: " + str(month))

    # Data for both models for the current month
    ecospace_data = data_by_month[month].get('log_PP1-DIA', [])
    qu39_data = data_by_month[month].get('logQU39', [])
    # ssc_data = data_by_month[month].get('log_ssc-DIA', [])

    # Plot the first model\
    plt.boxplot(ecospace_data, positions=[positions_ecospace[month - 1]], widths=widths, patch_artist=True,
                boxprops=dict(facecolor=colors['log_PP1-DIA'], color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                showfliers=False)

    # # Plot the second model
    # plt.boxplot(ssc_data, positions=[positions_ssc[month - 1]], widths=widths, patch_artist=True,
    #             boxprops=dict(facecolor=colors['log_ssc-DIA'], color='black'),
    #             medianprops=dict(color='black'),
    #             whiskerprops=dict(color='black'),
    #             capprops=dict(color='black'),
    #             showfliers=False)

    plt.boxplot(qu39_data, positions=[positions_qu39[month - 1]], widths=widths, patch_artist=True,
                boxprops=dict(facecolor=colors['logQU39'], color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                showfliers=False)

plt.title(f'Anomaly from Means Rel to Std Dev by Month for {taxon}')
plt.suptitle('')  # Suppress the default title to keep it clean
plt.xlabel('Month')
plt.ylabel('Mean')

# Adjust x-axis ticks and labels
plt.xticks(months, [f'Month {m}' for m in months])

# Create custom legend
# leg_codes = ['PP1-DIA', 'ssc-DIA', 'logQU39']
leg_codes = ['log_PP1-DIA', 'logQU39']
print(leg_codes)

for leg_code in leg_codes:
    print(leg_code)
    print(colors[leg_code])
handles = [plt.Line2D([0], [0], color=colors[leg_code], lw=4) for leg_code in leg_codes]
plt.legend(handles, source_names, loc='upper right')

plt.tight_layout()
plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_vs_Ecospace_2016_2018.png')
plt.show()

##########################

# Plot the histogram of measurement values by class
plt.figure(figsize=(15, 10))
classes = class_data['class'].unique()
for taxonomic_class in classes:
    if taxonomic_class == 'Bacillariophyceae':
        class_data = class_data[class_data['class'] == taxonomic_class]
        bins = np.logspace(np.log10(class_data['QU39'].min() + 1),
                           np.log10(class_data['QU39'].max() + 1), 50)
        plt.hist(class_data['QU39'] + 1, bins=bins, alpha=0.5,
                 label=taxonomic_class)  # Shift by 1 to avoid log(0)
        plt.xscale('log')

        # plt.hist(class_data['measurementValue'], bins=50, alpha=0.5, label=taxonomic_class, log=True)
        # plt.xscale('log')

plt.xlabel('Measurement Value')
plt.ylabel('Frequency')
plt.title('Distribution of Measurement Values by Class - QU39')
plt.legend(title='Taxonomic Class')
plt.grid(True)
plt.savefig('..//figs//' + 'Histo_Diatoms_PP_QU39_2016_2019.png')
plt.show()

plt.figure(figsize=(15, 10))
for taxonomic_class in classes:
    if taxonomic_class == 'Bacillariophyceae':
        class_data = class_data[class_data['class'] == taxonomic_class]
        bins = np.logspace(np.log10(class_data['ssc-DIA'].min() + 1),
                           np.log10(class_data['ssc-DIA'].max() + 1), 50)
        plt.hist(class_data['ssc-DIA'] + 1, bins=bins, alpha=0.5,
                 label=taxonomic_class)  # Shift by 1 to avoid log(0)
        plt.xscale('log')

        # plt.hist(class_data['measurementValue'], bins=50, alpha=0.5, label=taxonomic_class, log=True)
        # plt.xscale('log')

plt.xlabel('Measurement Value')
plt.ylabel('Frequency')
plt.title('Distribution of Measurement Values by Class - SalishSeaCast (SMELT)')
plt.legend(title='Taxonomic Class')
plt.grid(True)
plt.savefig('..//figs//' + 'Histo_Diatoms_PP_SSCDIA_2016_2019.png')
plt.show()

plt.figure(figsize=(15, 10))
for taxonomic_class in classes:
    if taxonomic_class == 'Bacillariophyceae':
        class_data = class_data[class_data['class'] == taxonomic_class]
        bins = np.logspace(np.log10(class_data['PP1-DIA'].min() + 1),
                           np.log10(class_data['PP1-DIA'].max() + 1), 50)
        plt.hist(class_data['PP1-DIA'] + 1, bins=bins, alpha=0.5,
                 label=taxonomic_class)  # Shift by 1 to avoid log(0)
        plt.xscale('log')

        # plt.hist(class_data['measurementValue'], bins=50, alpha=0.5, label=taxonomic_class, log=True)
        # plt.xscale('log')

plt.xlabel('Measurement Value')
plt.ylabel('Frequency')
plt.title('Distribution of Measurement Values by Class - Ecospace')
plt.legend(title='Taxonomic Class')
plt.grid(True)
plt.savefig('..//figs//' + 'Histo_Diatoms_PP_PP1DIA_2016_2019.png')
plt.show()

##################### histogram by month ##############
# Create a dictionary to hold the count of records by month
count_by_month = {month: 0 for month in months}

# Count the number of records for each month
for month in months:
    count_by_month[month] = len(class_data[class_data['month'] == month])

# Convert the count dictionary to a list for plotting
counts = [count_by_month[month] for month in months]

# Plot the histogram
plt.figure(figsize=(10, 6))
plt.bar(months, counts, color='skyblue', edgecolor='black')
plt.title('Count of Records by Month')
plt.xlabel('Month')
plt.ylabel('Count')
plt.xticks(months, [f'Month {m}' for m in months])
plt.grid(axis='y')

plt.tight_layout()
plt.show()

# # visuals, monthly pairwise boxplots
# taxonomic_class = 'Bacillariophyceae'
# model_code = 'PP1-DIA'
# plt.figure(figsize=(15, 10))
# # plt.subplot(1, 1, i)
# class_data = QU39_df[QU39_df['class'] == taxonomic_class]
# condition = QU39_df[model_code + '_anomaly_fr_mean'] > -999
# class_data = class_data.loc[condition]
# exclude_huge = False
# # if exclude_huge:
# #     median = class_data['measurementValue'].median()
# #     class_data = class_data[class_data['measurementValue'] <= median * 5]
# #     print(median)
# class_data.boxplot(column= model_code + '_anomaly_fr_mean_std_norm', by='month', grid=False, showfliers=False)
# plt.yscale('log')
# plt.title(f'Anomalies from mean by Month for {model_code} ({taxonomic_class})')
# plt.suptitle('')  # Suppress the default title to keep it clean
# plt.xlabel('Month')
# plt.ylabel('Anomaly (from Mean normalised w/ std)')
#
# # plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_2016_2019.png')
# plt.tight_layout()
# plt.show()
#
#
# # visuals, monthly pairwise boxplots
# taxonomic_class = 'Bacillariophyceae'
# model_code = 'PP1-DIA'
# plt.figure(figsize=(15, 10))
# # plt.subplot(1, 1, i)
# class_data = QU39_df[QU39_df['class'] == taxonomic_class]
# condition = QU39_df[model_code + '_anomaly_fr_mean'] > -999
# class_data = class_data.loc[condition]
# exclude_huge = False
# # if exclude_huge:
# #     median = class_data['measurementValue'].median()
# #     class_data = class_data[class_data['measurementValue'] <= median * 5]
# #     print(median)
# class_data.boxplot(column= model_code, by='month', grid=False, showfliers=False)
# # plt.yscale('log')
# plt.title(f'Means by Month for {model_code} ({taxonomic_class})')
# plt.suptitle('')  # Suppress the default title to keep it clean
# plt.xlabel('Month')
# plt.ylabel('Mean')
#
# # plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_2016_2019.png')
# plt.tight_layout()
# plt.show()
#
#
# # visuals, monthly pairwise boxplots
# taxonomic_class = 'Bacillariophyceae'
# model_code = 'PP1-DIA'
# plt.figure(figsize=(15, 10))
# # plt.subplot(1, 1, i)
# class_data = QU39_df[QU39_df['class'] == taxonomic_class]
# condition = QU39_df[model_code + '_anomaly_fr_mean'] > -999
# class_data = class_data.loc[condition]
# exclude_huge = False
# # if exclude_huge:
# #     median = class_data['measurementValue'].median()
# #     class_data = class_data[class_data['measurementValue'] <= median * 5]
# #     print(median)
# class_data.boxplot(column= model_code + '_anomaly_fr_mean', by='month', grid=False, showfliers=False)
# # plt.yscale('log')
# plt.title(f'Anomaly from Means by Month for {model_code} ({taxonomic_class})')
# plt.suptitle('')  # Suppress the default title to keep it clean
# plt.xlabel('Month')
# plt.ylabel('Mean')
#
# # plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_2016_2019.png')
# plt.tight_layout()
# plt.show()
#
#
#
# # visuals, monthly pairwise boxplots
# taxonomic_class = 'Bacillariophyceae'
# model_code = 'PP1-DIA'
# plt.figure(figsize=(15, 10))
# # plt.subplot(1, 1, i)
# class_data = QU39_df[QU39_df['class'] == taxonomic_class]
# condition = QU39_df[model_code + '_anomaly_fr_mean'] > -999
# class_data = class_data.loc[condition]
# exclude_huge = False
# # if exclude_huge:
# #     median = class_data['measurementValue'].median()
# #     class_data = class_data[class_data['measurementValue'] <= median * 5]
# #     print(median)
# class_data.boxplot(column= model_code + '_anomaly_fr_mean_std_norm', by='month', grid=False, showfliers=False)
# # plt.yscale('log')
# plt.title(f'Anomaly from Means Rel to Std Dev by Month for {model_code} ({taxonomic_class})')
# plt.suptitle('')  # Suppress the default title to keep it clean
# plt.xlabel('Month')
# plt.ylabel('Mean')
#
# # plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_2016_2019.png')
# plt.tight_layout()
# plt.show()


# months = np.arange(1, 13)
# positions_ecospace = np.arange(1, 13) - 0.2  # Adjust positions slightly to the left
# positions_qu39 = np.arange(1, 13) + 0.2  # Adjust positions slightly to the right
# widths = 0.35
#
# for i, (month_stats_ecospace, month_stats_ssc) in enumerate(zip(stats_ecospace, stats_ssc)):
#     pos_ecospace = positions_ecospace[i]
#     pos_qu39 = positions_qu39[i]
###################

###################
# Sample data preparation (replace this with your actual data)
# QU39_df = your_dataframe

########################################
# Sample data preparation (replace this with your actual data)
# QU39_df = your_dataframe
