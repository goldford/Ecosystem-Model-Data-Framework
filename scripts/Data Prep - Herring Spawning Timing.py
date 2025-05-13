

import pandas as pd

# Load the CSV file
file_path = "C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//3. Herring//Spawn Records Pearsall 2017//MODIFIED//spawning_dates.csv"
df = pd.read_csv(file_path)

# Display basic information and first few rows
df_info = df.info()
df_head = df.head()

# Convert 'Start' to datetime format
df['Start'] = pd.to_datetime(df['Start'], errors='coerce')

# Drop rows with missing or invalid 'Start' dates
df_valid = df.dropna(subset=['Start'])

# Add a day-of-year column
df_valid['DayOfYear'] = df_valid['Start'].dt.dayofyear

# Group by year and calculate the average day of year
average_doy_per_year = df_valid.groupby('Year')['DayOfYear'].mean().reset_index()
average_doy_per_year.columns = ['Year', 'Average_Start_DOY']

# Add a month column based on the 'Start' date
df_valid['Month'] = df_valid['Start'].dt.month

# Group by year and calculate both average DOY and average month
aggregated = df_valid.groupby('Year').agg({
    'DayOfYear': 'mean',
    'Month': 'mean'
}).reset_index()

aggregated.columns = ['Year', 'Average_Start_DOY', 'Average_Start_Month']


import matplotlib.pyplot as plt

# Plotting average DOY and average month over time
fig, ax1 = plt.subplots(figsize=(12, 6))

# Plot average DOY
ax1.plot(aggregated['Year'], aggregated['Average_Start_DOY'], label='Average DOY', marker='o')
ax1.set_xlabel('Year')
ax1.set_ylabel('Average Day of Year', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

# Create a second y-axis for the average month
ax2 = ax1.twinx()
ax2.plot(aggregated['Year'], aggregated['Average_Start_Month'], label='Average Month', color='tab:orange', linestyle='--', marker='s')
ax2.set_ylabel('Average Month (numeric)', color='tab:orange')
ax2.tick_params(axis='y', labelcolor='tab:orange')

# Title and legends
plt.title('Trend of Average Spawning Start DOY and Month Over Time')
fig.tight_layout()
plt.grid(True)
plt.show()