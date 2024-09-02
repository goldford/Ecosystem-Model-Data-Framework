# Created by G Oldford
# July 18 2024
# Purpose: join multiple tables from dataset of phytoplankton concentrations by group
# Source:
#    Del Bel Belluz, J. (2024). Protistan plankton time series from the northern Salish Sea and central coast,
#    British Columbia, Canada—Ocean Biodiversity Information System (a62c37c4-6fbe-4a93-8c87-c67dda36436c)
#    [Dataset]. https://obis.org/dataset/071a38b3-0d30-47ba-85aa-d98313460075
# Input:
#    Tables from dataset, microscopy ID'd phytoplankton with concentration (cells per litre)
# Output:
#    wide and long joined table, anomalies calculated, only for station QU39
#

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the data from the files
occurrence_path = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\ORIGINAL\occurrence.txt'
extendedmeasurementorfact_path = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\ORIGINAL\extendedmeasurementorfact.txt'
event_path = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\ORIGINAL\event.txt'
out_path = r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED\QU39_joined.csv'

occurrence_df = pd.read_csv(occurrence_path, sep='\t')
extendedmeasurementorfact_df = pd.read_csv(extendedmeasurementorfact_path, sep='\t')
event_df = pd.read_csv(event_path, sep='\t')

# Filter the data
filtered_event_df = event_df[event_df['locationID'] == 'QU39']
filtered_measurement_df = extendedmeasurementorfact_df[extendedmeasurementorfact_df['measurementType'] == 'Abundance of phytoplankton']
filtered_measurement_df['measurementValue'] = pd.to_numeric(filtered_measurement_df['measurementValue'], errors='coerce')

# Merge the data
merged_df = pd.merge(filtered_measurement_df, occurrence_df, on='occurrenceID')
merged_df = pd.merge(merged_df, filtered_event_df, on='eventID')

# Convert eventDate to datetime and extract month and year
merged_df['eventDate'] = pd.to_datetime(merged_df['eventDate'])
merged_df['year'] = merged_df['eventDate'].dt.year
merged_df['month'] = merged_df['eventDate'].dt.month

# Define seasons
def get_season(month):
    if month in [12, 1, 2]:
        return 'Winter'
    elif month in [3, 4, 5]:
        return 'Spring'
    elif month in [6, 7, 8]:
        return 'Summer'
    elif month in [9, 10, 11]:
        return 'Fall'

merged_df['season'] = merged_df['month'].apply(get_season)

# Function to calculate anomalies
def calculate_anomalies(df):
    # Compute monthly means for each year and class
    monthly_year_means = df.groupby(['year', 'month', 'class'])['measurementValue'].mean().reset_index()
    monthly_year_medians = df.groupby(['year', 'month', 'class'])['measurementValue'].median().reset_index()

    # Compute annual mean and standard deviation from the monthly means for each class
    monthly_means = monthly_year_means.groupby(['class', 'month'])['measurementValue'].mean()
    monthly_medians = monthly_year_medians.groupby(['class', 'month'])['measurementValue'].median()
    annual_means = monthly_means.groupby(['class']).mean()
    annual_std = monthly_means.groupby(['class']).std()
    annual_medians = monthly_medians.groupby(['class']).median()

    # Merge the annual means and standard deviations back into the original dataframe
    df = df.merge(annual_means, on=['class'], suffixes=('', '_annual_mean'))
    df = df.merge(annual_std, on=['class'], suffixes=('', '_annual_std'))
    df = df.merge(annual_medians, on=['class'], suffixes=('', '_annual_median'))


    # Calculate the anomaly
    df['anomaly_fr_mean'] = (df['measurementValue'] - df['measurementValue_annual_mean']) / df['measurementValue_annual_std']
    df['anomaly_fr_median'] = (df['measurementValue'] - df['measurementValue_annual_median']) / df[
        'measurementValue_annual_std']
    df['log_anomaly_fr_mean'] = np.log(df['anomaly_fr_mean'] + 1 - df['anomaly_fr_mean'].min())

    return df

# Calculate anomalies for the merged data
merged_df = calculate_anomalies(merged_df)
merged_df.to_csv(out_path)

# Summarize the count of sampling events by each month
sampling_events_count = merged_df.groupby('month')['eventID'].nunique().reset_index()
sampling_events_count.columns = ['month', 'sampling_events_count']

# Plot the histogram of sampling events by month
plt.figure(figsize=(10, 6))
plt.bar(sampling_events_count['month'], sampling_events_count['sampling_events_count'], color='blue', edgecolor='black')
plt.xlabel('Month')
plt.ylabel('Number of Sampling Events')
plt.title('Number of Sampling Events by Month')
plt.xticks(sampling_events_count['month'])
plt.grid(axis='y')
plt.show()

# Summarize the count of sampling events by each season
sampling_events_count = merged_df.groupby('season')['eventID'].nunique().reset_index()
sampling_events_count.columns = ['season', 'sampling_events_count']
print(sampling_events_count)

# Plot the histogram of sampling events by month
plt.figure(figsize=(10, 6))
plt.bar(sampling_events_count['season'], sampling_events_count['sampling_events_count'], color='blue', edgecolor='black')
plt.xlabel('Season')
plt.ylabel('Number of Sampling Events')
plt.title('Number of Sampling Events by Season')
plt.xticks(sampling_events_count['season'])
plt.grid(axis='y')
plt.show()

# Plot the histogram of measurement values by class
plt.figure(figsize=(15, 10))
classes = merged_df['class'].unique()
for taxonomic_class in classes:
    if taxonomic_class == 'Bacillariophyceae':
        class_data = merged_df[merged_df['class'] == taxonomic_class]

        bins = np.logspace(np.log10(class_data['measurementValue'].min() + 1),
                           np.log10(class_data['measurementValue'].max() + 1), 50)
        plt.hist(class_data['measurementValue'] + 1, bins=bins, alpha=0.5,
                 label=taxonomic_class)  # Shift by 1 to avoid log(0)
        plt.xscale('log')

        # plt.hist(class_data['measurementValue'], bins=50, alpha=0.5, label=taxonomic_class, log=True)
        # plt.xscale('log')

plt.xlabel('Measurement Value')
plt.ylabel('Frequency')
plt.title('Distribution of Measurement Values by Class')
plt.legend(title='Taxonomic Class')
plt.grid(True)
plt.savefig('..//figs//' + 'Histo_Diatoms_PP_QU39_2016_2019.png')
plt.show()


# Plot anomalies by season for each class
classes = merged_df['class'].unique()
plt.figure(figsize=(15, 10))
for i, taxonomic_class in enumerate(classes, 1):
    if taxonomic_class == 'Bacillariophyceae':
        plt.subplot(len(classes), 1, i)
        class_data = merged_df[merged_df['class'] == taxonomic_class]
        exclude_huge = False
        if exclude_huge:
            median = class_data['measurementValue'].median()
            class_data = class_data[class_data['measurementValue'] <= median * 5]
            print(median)
        class_data.boxplot(column='log_anomaly_fr_mean', by='season', grid=False, showfliers=False)
        # plt.yscale('log')
        plt.title(f'Log-Anomalies by Season for {taxonomic_class}')
        plt.suptitle('')  # Suppress the default title to keep it clean
        plt.xlabel('Season')
        plt.ylabel('Anomaly')

# Plot anomalies by month for each class
classes = merged_df['class'].unique()
plt.figure(figsize=(15, 10))
for i, taxonomic_class in enumerate(classes, 1):
    if taxonomic_class == 'Bacillariophyceae':
        plt.subplot(len(classes), 1, i)
        class_data = merged_df[merged_df['class'] == taxonomic_class]
        exclude_huge = False
        if exclude_huge:
            median = class_data['measurementValue'].median()
            class_data = class_data[class_data['measurementValue'] <= median * 5]
            print(median)
        class_data.boxplot(column='anomaly_fr_median', by='month', grid=False, showfliers=False)
        # plt.yscale('log')
        plt.title(f'Log-Anomalies by Month for {taxonomic_class}')
        plt.suptitle('')  # Suppress the default title to keep it clean
        plt.xlabel('Month')
        plt.ylabel('Anomaly (from Median)')

plt.savefig('..//figs//' + 'Diatoms_Monthly_QU39_2016_2019.png')
plt.tight_layout()
plt.show()

