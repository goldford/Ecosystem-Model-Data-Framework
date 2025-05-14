# -------------------------------------------------------------------
# Script Name: Ecospace Analysis - QU39 vs Ecospace vs SSC
# Author: G. Oldford
# Last Updated: May 2025
#
# Description:
#   This script compares observational phytoplankton data from QU39
#   (2016–2019) with matched model outputs from:
#     - Ecospace (EwE) – functional group biomass
#     - SalishSeaCast (SSC/SMELT) – phytoplankton concentration
#
#   The comparison is visualized using boxplots of normalized anomalies
#   (relative to climatological monthly means and standard deviations).
#
# Data Source:
#   Del Bel Belluz, J. (2024). Protistan plankton time series from the
#   northern Salish Sea and central coast, BC, Canada.
#   OBIS Dataset: https://obis.org/dataset/071a38b3-0d30-47ba-85aa-d98313460075
#
# Inputs:
#   - CSV file (QU39_joined_matchtoEcospace_SSC_<model_run>.csv)
#     containing bottle-sampled phytoplankton counts matched to
#     Ecospace and SSC model outputs.
#
# Outputs:
#   - Monthly boxplots comparing QU39 vs Ecospace vs SSC (per group)
#   - Histograms of raw and log-transformed concentration distributions
#   - Count of records by month
#
# Notes:
#   - Anomalies are computed as:
#       (value - monthly mean) / monthly std deviation
#   - Log-transformed values are used for interpretability
#   - Figures are saved to the `figs/` directory
# -------------------------------------------------------------------


# ----------------------------- Configuration ---------------------------- #
import os
import pandas as pd
import calendar
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_squared_error, mean_absolute_error
from scipy.stats import pearsonr


# Toggle whether SSC outputs are included in the comparison
include_SSC = True

# Model run code (used in file paths and plot names)
model_run = 'SC88_1'

# Input data path
pathfile_QU39_SSC_Ecospace = (
    r'C:\Users\Greig\Sync\6. SSMSP Model\Model Greig\Data\28. Phytoplankton'
    r'\Phyto Concent del Bel Belluz 2024 2016 - 2019\MODIFIED'
    r'\QU39_joined_matchtoEcospace_SSC_' + model_run + '.csv'
)

# Output figure path
out_p = os.path.normpath("..//..//figs//") + os.sep

# Standard placeholder for bad data
FILL_VALUE = -999

# Month ordering
months = np.arange(1, 13)
# ------------------------------------------------------------------------ #


def calculate_anomalies(df, valuefield, groupon):
    """
    Calculate anomalies of a value field normalized by climatological mean and std dev.

    Parameters:
        df : pd.DataFrame
            The dataframe with columns including valuefield, year, month, and groupon.
        valuefield : str
            The column to compute anomalies on.
        groupon : str
            The taxonomic grouping field (e.g., 'class' or 'family').

    Returns:
        pd.DataFrame with new columns for anomalies and climatology.
    """
    df = df.copy()

    if 'level_0' in df.columns:
        df = df.drop(columns=['level_0'])

    if 'eventID' in df.columns:
        df = df.drop(columns=['eventID'])

    # Filter valid values
    df = df[df[valuefield] > FILL_VALUE]

    # Compute monthly climatology grouped by year/month and groupon
    numeric_cols = df.select_dtypes(include=np.number).columns.difference(
        ['year', 'month', groupon, 'closest_ecospace_time'])
    monthly_stats = df.groupby(['year', 'month', groupon, 'closest_ecospace_time'])[numeric_cols].mean().reset_index()

    # Compute monthly and annual climatologies
    monthly_means = monthly_stats.groupby([groupon, 'month'])[valuefield].mean()
    monthly_medians = monthly_stats.groupby([groupon, 'month'])[valuefield].median()
    annual_mean = monthly_means.groupby(groupon).mean()
    annual_std = monthly_means.groupby(groupon).std()
    annual_median = monthly_medians.groupby(groupon).median()

    # Merge into main dataframe
    df = df.merge(annual_mean.rename(valuefield + '_annual_mean'), on=groupon, how='left')
    df = df.merge(annual_std.rename(valuefield + '_annual_std'), on=groupon, how='left')
    df = df.merge(annual_median.rename(valuefield + '_annual_median'), on=groupon, how='left')

    # Anomaly calculations
    df[valuefield + '_anomaly_fr_mean'] = df[valuefield] - df[valuefield + '_annual_mean']
    df[valuefield + '_anomaly_fr_median'] = df[valuefield] - df[valuefield + '_annual_median']
    df[valuefield + '_anomaly_fr_mean_std_norm'] = (
            df[valuefield + '_anomaly_fr_mean'] / df[valuefield + '_annual_std']
    )
    df[valuefield + '_anomaly_fr_median_std_norm'] = (
            df[valuefield + '_anomaly_fr_median'] / df[valuefield + '_annual_std']
    )
    df[valuefield + '_log_anomaly_fr_mean_std_norm'] = np.log1p(
        df[valuefield + '_anomaly_fr_mean_std_norm'] - df[valuefield + '_anomaly_fr_mean_std_norm'].min()
    )

    return df

def aggregate_monthly_means(df, value_fields, groupby_time='closest_ecospace_time'):
    """
    Aggregate the log-transformed fields by year and month.

    Parameters:
        df : pd.DataFrame
        value_fields : list of str
            Columns to aggregate (e.g., ['logQU39', 'log_PP2-NAN', 'log_ssc-FLA'])
        groupby_time : str
            Time field to group by (typically 'closest_ecospace_time')

    Returns:
        pd.DataFrame with columns: year, month, mean_<fieldname> for each field
    """
    df = df.copy()
    df[groupby_time] = pd.to_datetime(df[groupby_time])
    df['year'] = df[groupby_time].dt.year
    df['month'] = df[groupby_time].dt.month

    grouped = df.groupby(['year', 'month'])[value_fields].mean().reset_index()
    grouped.columns = ['year', 'month'] + [f'mean_{col}' for col in value_fields]
    return grouped

def compute_model_stats(obs, model):
    """
    Compute common statistical comparisons for model skill assessment.

    Parameters:
        obs : np.ndarray
        model : np.ndarray

    Returns:
        dict with Pearson R, std devs, centered RMSE, RMSE, bias, MAE, WSS, N
    """
    mask = ~np.isnan(obs) & ~np.isnan(model)
    obs = obs[mask]
    model = model[mask]

    obs_std = np.std(obs)
    model_std = np.std(model)
    r, _ = pearsonr(obs, model)

    bias = np.mean(model - obs)
    mae = mean_absolute_error(obs, model)
    rmse = np.sqrt(np.mean((model - obs) ** 2))

    # Centered RMSE (remove means)
    centered_rmse = np.sqrt(np.mean(((model - np.mean(model)) - (obs - np.mean(obs))) ** 2))

    # Willmott’s Skill Score (WSS)
    obs_mean = np.mean(obs)
    denom = np.sum((np.abs(model - obs_mean) + np.abs(obs - obs_mean)) ** 2)
    wss = 1 - (np.sum((model - obs) ** 2) / denom) if denom != 0 else np.nan

    return {
        'R': r,
        'obs_std': obs_std,
        'model_std': model_std,
        'centered_RMSE': centered_rmse,
        'RMSE': rmse,
        'Bias': bias,
        'MAE': mae,
        'WSS': wss,
        'N': len(obs)
    }


def run_model_statistics(df, model_obs_pairs):
    results = []

    for pair in model_obs_pairs:
        taxon = pair['taxon_filter']
        group = pair['group']
        obs_field = pair['obs']
        ecospace_field = pair['ecospace']
        ssc_field = pair['ssc']

        # Filter for this group
        df_sub = df[df['class'].isin(taxon if isinstance(taxon, list) else [taxon])]
        df_sub = df_sub[pd.to_datetime(df_sub['closest_ecospace_time']).dt.year >= 1980]

        # Aggregate to year-month
        monthly = aggregate_monthly_means(df_sub, [obs_field, ecospace_field] + ([ssc_field] if ssc_field else []))

        # Calculate stats
        ecospace_stats = compute_model_stats(monthly[f'mean_{obs_field}'].values, monthly[f'mean_{ecospace_field}'].values)
        ecospace_stats.update({'Group': group, 'Model': 'Ecospace'})

        results.append(ecospace_stats)

        if ssc_field:
            ssc_stats = compute_model_stats(monthly[f'mean_{obs_field}'].values, monthly[f'mean_{ssc_field}'].values)
            ssc_stats.update({'Group': group, 'Model': 'SSC'})
            results.append(ssc_stats)

    return pd.DataFrame(results)


def plot_monthly_boxplot(
    df,
    field_codes,
    taxon_name,
    groupon,
    plot_label,
    source_names,
    color_map,
    output_file=None,
    include_ssc=False,
    label_position='(a)',
    use_log_anomaly=False
):
    """
    Generate a monthly boxplot comparing anomalies from different sources.

    Parameters:
        df : pd.DataFrame
            Data with anomaly fields already calculated.
        field_codes : list of str
            Fields to plot (e.g., ['log_PP1-DIA', 'logQU39', 'log_ssc-DIA']).
        taxon_name : str
            Name of the taxon to filter for (used with groupon column).
        groupon : str
            Column to group/filter by (e.g., 'class').
        plot_label : str
            Text label for the plot title and output filename.
        source_names : list of str
            Names to show in the legend (should match field_codes order).
        color_map : dict
            Dictionary mapping field_codes to plot colors.
        output_file : str or None
            Full file path to save the plot. If None, does not save.
        include_ssc : bool
            If True, includes the SSC model in plot alignment.
        label_position : str
            Text label to show in plot corner (e.g., '(a)', '(b)').
        use_log_anomaly : bool
            If True, use the log-transformed normalized anomaly field.
    """

    # Setup
    months = np.arange(1, 13)
    if include_ssc and len(field_codes) == 3:
        positions = [months - 0.25, months + 0.25, months]
    elif len(field_codes) == 2:
        positions = [months - 0.2, months + 0.2]
    else:
        positions = [months]

    widths = 0.35 if len(field_codes) == 2 else 0.2

    fig, ax = plt.subplots(figsize=(10, 6) if len(field_codes) == 2 else (15, 10))

    # Filter for taxon of interest
    subset = df[df[groupon].isin([taxon_name] if isinstance(taxon_name, str) else taxon_name)]
    condition = subset[field_codes[0] + '_anomaly_fr_mean'] > FILL_VALUE
    subset = subset.loc[condition]

    # Organize data for plotting
    data_by_month = {month: {} for month in months}
    for field in field_codes:
        for month in months:
            mdata = subset[subset['month'] == month]
            field_to_plot = field + ('_log_anomaly_fr_mean_std_norm' if use_log_anomaly else '_anomaly_fr_mean_std_norm')
            data_by_month[month][field] = mdata[field_to_plot].dropna().values

    # Plot each group
    for i, field in enumerate(field_codes):
        for month in months:
            plt.boxplot(
                data_by_month[month].get(field, []),
                positions=[positions[i][month - 1]],
                widths=widths,
                patch_artist=True,
                boxprops=dict(facecolor=color_map[field], color='black'),
                medianprops=dict(color='black'),
                whiskerprops=dict(color='black'),
                capprops=dict(color='black'),
                showfliers=False
            )

    # Labels and legend
    plt.title(f'Anomaly from Means Rel to Std Dev by Month for {plot_label}')
    plt.xlabel('Month')
    plt.ylabel('Normalized Anomaly')
    plt.xticks(months, [calendar.month_abbr[m] for m in months])
    handles = [plt.Line2D([0], [0], color=color_map[f], lw=4) for f in field_codes]
    plt.legend(handles, source_names, loc='upper right')
    plt.text(0.02, 0.95, label_position, transform=ax.transAxes, fontsize=14, verticalalignment='top')
    plt.tight_layout()

    if output_file:
        plt.savefig(output_file)
        print(f"[Saved] {output_file}")

    plt.show()


def plot_log_histogram(df, field, group_filter, group_column, title, output_file):
    """
    Plot log-scaled histogram for a given field and taxonomic group.

    Parameters:
        df : pd.DataFrame
        field : str - column name for value (e.g., 'QU39', 'PP1-DIA')
        group_filter : str - taxonomic group (e.g., 'Bacillariophyceae')
        group_column : str - column name for group ('class', etc.)
        title : str - title of the plot
        output_file : str - path to save the plot
    """
    subset = df[df[group_column] == group_filter]
    values = subset[field] + 1  # avoid log(0)
    bins = np.logspace(np.log10(values.min()), np.log10(values.max()), 50)

    plt.figure(figsize=(10, 6))
    plt.hist(values, bins=bins, alpha=0.7, label=group_filter)
    plt.xscale('log')
    plt.xlabel('Measurement Value')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"[Saved] {output_file}")
    plt.show()


def plot_monthly_sample_counts(df, group_filter, group_column, output_file):
    """
    Plot a bar chart of sample count by month for a given group.

    Parameters:
        df : pd.DataFrame
        group_filter : str - e.g. 'Bacillariophyceae'
        group_column : str - e.g. 'class'
        output_file : str - path to save the plot
    """
    subset = df[df[group_column] == group_filter]
    count_by_month = subset['month'].value_counts().reindex(np.arange(1, 13), fill_value=0)

    plt.figure(figsize=(10, 5))
    plt.bar(count_by_month.index, count_by_month.values, color='skyblue', edgecolor='black')
    plt.title(f'Count of Records by Month for {group_filter}')
    plt.xlabel('Month')
    plt.ylabel('Sample Count')
    plt.xticks(count_by_month.index, [calendar.month_abbr[m] for m in count_by_month.index])
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(output_file)
    print(f"[Saved] {output_file}")
    plt.show()


def main():
    print("Loading dataset...")
    df = pd.read_csv(pathfile_QU39_SSC_Ecospace)
    df['DateTime'] = pd.to_datetime(df['DateTime'], format='%Y-%m-%d %H:%M:%S')
    df = df[pd.to_datetime(df['closest_ecospace_time']).dt.year >= 1980]

    # Derive values
    df['QU39'] = df['measurementValue']
    df['logQU39'] = np.log(df['measurementValue'].clip(lower=1))  # +1 for safety
    df['log_PP1-DIA'] = np.log(df['PP1-DIA'].clip(lower=1))
    df['log_PP2-NAN'] = np.log(df['PP2-NAN'].clip(lower=1))
    df['log_PP3-PIC'] = np.log(df['PP3-PIC'].clip(lower=1))
    df['log_PZ2-DIN'] = np.log(df['PZ2-DIN'].clip(lower=1))
    df['log_ssc-DIA'] = np.log(df['ssc-DIA'].clip(lower=1))
    df['log_ssc-FLA'] = np.log(df['ssc-FLA'].clip(lower=1))
    df['log_ssc-CIL'] = np.log(df['ssc-CIL'].clip(lower=1))

    # Compute anomalies
    df = calculate_anomalies(df, 'logQU39', 'class')
    df = calculate_anomalies(df, 'log_PP1-DIA', 'class')
    df = calculate_anomalies(df, 'log_ssc-DIA', 'class')
    df = calculate_anomalies(df, 'log_PP2-NAN', 'class')
    df = calculate_anomalies(df, 'log_ssc-FLA', 'class')
    df = calculate_anomalies(df, 'log_PP3-PIC', 'class')
    df = calculate_anomalies(df, 'log_PZ2-DIN', 'class')

    # Diatoms
    plot_monthly_boxplot(
        df, ['log_PP1-DIA', 'log_ssc-DIA', 'logQU39'], 'Bacillariophyceae', 'class',
        'Diatoms', ['Ecospace', 'SMELT', 'QU39'],
        {'log_PP1-DIA': 'blue', 'log_ssc-DIA': 'pink', 'logQU39': 'orange'},
        output_file=os.path.join(out_p, model_run + '_Diatoms_Monthly_QU39_vs_Ecospace_vs_SSC_2016_2018.png'),
        include_ssc=True, label_position='(a)'
    )

    # Nanoflagellates
    nanogroup = [
        'Choanoflagellatea', 'Dictyochophyceae', 'Cryptophyceae', 'Metromonadea',
        'Chrysophyceae', 'Telonemea', 'Chlorodendrophyceae', 'Bicosoecophyceae',
        'Xanthophyceae', 'Coccolithophyceae', 'Euglenophyceae', 'Raphidophyceae'
    ]
    df_nan = df[df['class'].isin(nanogroup)]
    plot_monthly_boxplot(
        df_nan, ['log_PP2-NAN', 'log_ssc-FLA', 'logQU39'], nanogroup, 'class',
        'Nanophytoplankton', ['Ecospace', 'SMELT', 'QU39'],
        {'log_PP2-NAN': 'blue', 'log_ssc-FLA': 'pink', 'logQU39': 'orange'},
        output_file=os.path.join(out_p, 'Nanophytoplankton_' + model_run + '_Monthly_QU39_vs_Ecospace_vs_SSC_2016_2018.png'),
        include_ssc=True, label_position='(b)'
    )

    # Picoplankton
    plot_monthly_boxplot(
        df, ['log_PP3-PIC', 'logQU39'], 'Pyramimonadophyceae', 'class',
        'Picoplankton', ['Ecospace', 'QU39'],
        {'log_PP3-PIC': 'blue', 'logQU39': 'orange'},
        output_file=os.path.join(out_p, 'Picoplankton_' + model_run + '_Monthly_QU39_vs_Ecospace_2016_2018.png'),
        include_ssc=False, label_position='(c)'
    )

    # Dinoflagellates
    plot_monthly_boxplot(
        df, ['log_PZ2-DIN', 'logQU39'], 'Dinophyceae', 'class',
        'Dinoflagellates', ['Ecospace', 'QU39'],
        {'log_PZ2-DIN': 'blue', 'logQU39': 'orange'},
        output_file=os.path.join(out_p, 'Dinoflagellates_' + model_run + '_Monthly_QU39_vs_Ecospace_2016_2018.png'),
        include_ssc=False, label_position='(d)'
    )

    # Histograms
    plot_log_histogram(df, 'QU39', 'Bacillariophyceae', 'class', 'Distribution - QU39', out_p + 'Histo_QU39_Diatoms.png')
    plot_log_histogram(df, 'PP1-DIA', 'Bacillariophyceae', 'class', 'Distribution - Ecospace', out_p + 'Histo_PP1_DIA.png')
    plot_log_histogram(df, 'ssc-DIA', 'Bacillariophyceae', 'class', 'Distribution - SSC', out_p + 'Histo_SSC_DIA.png')

    # Monthly sample count
    plot_monthly_sample_counts(df, 'Bacillariophyceae', 'class', out_p + 'Monthly_Sample_Count_Diatoms.png')

    print("✅ All plots generated.")
    # Observation and model comparison definitions
    model_obs_pairs = [
        {
            'group': 'Diatoms',
            'taxon_filter': 'Bacillariophyceae',
            'obs': 'logQU39',
            'ecospace': 'log_PP1-DIA',
            'ssc': 'log_ssc-DIA'
        },
        {
            'group': 'Nanophytoplankton',
            'taxon_filter': [
                'Choanoflagellatea', 'Dictyochophyceae', 'Cryptophyceae', 'Metromonadea',
                'Chrysophyceae', 'Telonemea', 'Chlorodendrophyceae', 'Bicosoecophyceae',
                'Xanthophyceae', 'Coccolithophyceae', 'Euglenophyceae', 'Raphidophyceae'
            ],
            'obs': 'logQU39',
            'ecospace': 'log_PP2-NAN',
            'ssc': 'log_ssc-FLA'
        },
        {
            'group': 'Picoplankton',
            'taxon_filter': 'Pyramimonadophyceae',
            'obs': 'logQU39',
            'ecospace': 'log_PP3-PIC',
            'ssc': None
        },
        {
            'group': 'Dinoflagellates',
            'taxon_filter': 'Dinophyceae',
            'obs': 'logQU39',
            'ecospace': 'log_PZ2-DIN',
            'ssc': None
        }
    ]
    stats_df = run_model_statistics(df, model_obs_pairs)


    print('stats generated!')
    print(stats_df)
    print('DONE')


if __name__ == "__main__":
    main()