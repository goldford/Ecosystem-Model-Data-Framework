import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data
file_path = 'C://Users//Greig//Sync//6. SSMSP Model//Model Greig//Data//28. Phytoplankton//Phytoplankton Salish Sea Nemcek2023 2015-2019//MODIFIED//Nemcek_Supp_Data.csv'
df = pd.read_csv(file_path)

# Convert 'Date.Time' to datetime
df['Date.Time'] = pd.to_datetime(df['Date.Time'])

# Extract month and year from the date
df['Month'] = df['Date.Time'].dt.month
df['Year'] = df['Date.Time'].dt.year

# Combine 'Prasinophytes', 'Cryptophytes', 'Haptophytes', 'Dictyochophytes' into 'PP2-NAN'
df['PP2-NAN'] = df[['Prasinophytes', 'Cryptophytes', 'Haptophytes', 'Dictyochophytes', 'Raphidophytes']].sum(axis=1)

# List of functional groups to include, including the new 'PP2-NAN' group
functional_groups = [
    'PP2-NAN', 'Dinoflagellates',
    'Cyanobacteria', 'total diatoms'
]

# Filter out rows where 'sdomain' is empty
df_filtered = df[df['sdomain'].isin(['SGS', 'SGN', 'SGI'])]

# Extract month from Date.Time
df_filtered['Month'] = df_filtered['Date.Time'].dt.month

# Group by 'Month' and 'sdomain' and count the number of samples
samples_by_month_sdomain = df_filtered.groupby(['Month', 'sdomain']).size().reset_index(name='Sample Count')



# Bar Plot
fig, ax = plt.subplots(figsize=(12, 6))

for sdomain in samples_by_month_sdomain['sdomain'].unique():
    subset = samples_by_month_sdomain[samples_by_month_sdomain['sdomain'] == sdomain]
    ax.bar(subset['Month'], subset['Sample Count'], label=sdomain)

ax.set_title('Number of Samples by Month and Sdomain')
ax.set_xlabel('Month')
ax.set_ylabel('Sample Count')
ax.legend(title='Sdomain')
ax.set_xticks(range(1, 13))
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])

plt.show()

# Heatmap
samples_pivot = samples_by_month_sdomain.pivot(index='sdomain', columns='Month', values='Sample Count')

fig, ax = plt.subplots(figsize=(12, 6))
cax = ax.matshow(samples_pivot, cmap='YlGnBu')

fig.colorbar(cax)

ax.set_title('Heatmap of Samples by Month and Sdomain')
ax.set_xlabel('Month')
ax.set_ylabel('Sdomain')
ax.set_xticks(range(len(samples_pivot.columns)))
ax.set_xticklabels(['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'])
ax.set_yticks(range(len(samples_pivot.index)))
ax.set_yticklabels(samples_pivot.index)

plt.show()