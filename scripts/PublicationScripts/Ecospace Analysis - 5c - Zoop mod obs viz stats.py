import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# === Configuration ===
path_evalout = "C:/Users/Greig/Documents/github/Ecosystem-Model-Data-Framework/data/evaluation"
ecospace_code = 'Scv88_1'
input_file = f"Zooplankton_matched_to_model_out_{ecospace_code}.csv"

# === Load Data ===
df = pd.read_csv(f"{path_evalout}/{input_file}")

# === Define groups ===
# omitting Jellyfish!
obs_cols = [col for col in df.columns if col in [
    'ZC1-EUP', 'ZC2-AMP', 'ZC3-DEC', 'ZC4-CLG', 'ZC5-CSM',
    'ZS2-CTH', 'ZS3-CHA', 'ZS4-LAR', 'ZF1-ICT']]
model_cols = [f"EWE-{col}" for col in obs_cols]

# === Replace zeros using Perry et al. 2021 technique ===
for col in obs_cols:
    nonzero_vals = df[col][df[col] > 0]
    if not nonzero_vals.empty:
        min_nonzero = nonzero_vals.min()
        df[col] = df[col].apply(lambda x: np.random.uniform(0, 0.5 * min_nonzero) if x == 0 else x)

# === Step 1: Aggregate by closest_ecospace_time ===
df_agg = df.groupby(['closest_ecospace_time', 'Index'])[obs_cols + model_cols].mean().reset_index()

# === Step 2: Compute Total Biomass Fields ===
df_agg['TOT'] = df_agg[obs_cols].sum(axis=1)
df_agg['EWE-TOT'] = df_agg[model_cols].sum(axis=1)

# === Step 3: Scatter Plots ===
def plot_log_scatter(x, y, label, ax):
    ax.scatter(np.log10(x + 1e-6), np.log10(y + 1e-6), alpha=0.6)
    ax.plot([-3, 3], [-3, 3], 'k--', lw=1)
    ax.set_xlabel("log10(Observed Biomass + 1e-6)")
    ax.set_ylabel("log10(Modelled Biomass + 1e-6)")
    ax.set_title(label)
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)

# Point-wise total scatter
fig, ax = plt.subplots(figsize=(6,6))
plot_log_scatter(df_agg['TOT'], df_agg['EWE-TOT'], 'Total Biomass', ax)
plt.tight_layout()
plt.show()

# Group-wise scatter plots
n = len(obs_cols)
fig, axs = plt.subplots(nrows=int(np.ceil(n/2)), ncols=2, figsize=(12, 2.5 * int(np.ceil(n/2))))
axs = axs.flatten()

for i, col in enumerate(obs_cols):
    plot_log_scatter(df_agg[col], df_agg[f"EWE-{col}"], col, axs[i])

for j in range(i+1, len(axs)):
    fig.delaxes(axs[j])

plt.tight_layout()
plt.show()
print('done the one-to-one comparison')


# === Additional Aggregation by Station, Season, Year, Region ===
df_station = df.groupby(['Station', 'Year', 'Season', 'Region'])[obs_cols + model_cols].mean().reset_index()
df_station['TOT'] = df_station[obs_cols].sum(axis=1)
df_station['EWE-TOT'] = df_station[model_cols].sum(axis=1)

# Apply log10 transformation before computing skill statistics
for col in obs_cols + ['TOT']:
    df_station[f'log10_{col}'] = np.log10(df_station[col] + 1e-6)
    model_col = f"EWE-{col}" if col != 'TOT' else 'EWE-TOT'
    df_station[f'log10_{model_col}'] = np.log10(df_station[model_col] + 1e-6)

# Station-level total scatter by Region
fig, ax = plt.subplots(figsize=(6,6))
sns.scatterplot(
    x=df_station['log10_TOT'],
    y=df_station['log10_EWE-TOT'],
    hue=df_station['Region'],
    alpha=0.7,
    ax=ax
)
ax.plot([-3, 3], [-3, 3], 'k--', lw=1)
ax.set_xlabel("log10(Observed Biomass + 1e-6)")
ax.set_ylabel("log10(Modelled Biomass + 1e-6)")
ax.set_title("Total Biomass (Station-Averaged by Region)")
ax.set_xlim(-3, 3)
ax.set_ylim(-3, 3)
plt.tight_layout()
plt.show()

# Station-level group-wise scatter plots by Region
fig, axs = plt.subplots(nrows=int(np.ceil(n/2)), ncols=2, figsize=(12, 2.5 * int(np.ceil(n/2))))
axs = axs.flatten()

for i, col in enumerate(obs_cols):
    sns.scatterplot(
        x=df_station[f'log10_{col}'],
        y=df_station[f'log10_EWE-{col}'],
        hue=df_station['Region'],
        alpha=0.7,
        ax=axs[i]
    )
    axs[i].plot([-3, 3], [-3, 3], 'k--', lw=1)
    axs[i].set_xlabel("log10(Obs + 1e-6)")
    axs[i].set_ylabel("log10(Mod + 1e-6)")
    axs[i].set_title(col + " (Station-Averaged)")
    axs[i].set_xlim(-3, 3)
    axs[i].set_ylim(-3, 3)

for j in range(i+1, len(axs)):
    fig.delaxes(axs[j])

plt.tight_layout()
plt.show()

# === Compute model skill statistics ===
def compute_skill_statistics(obs, mod):
    obs = np.array(obs)
    mod = np.array(mod)
    mask = ~np.isnan(obs) & ~np.isnan(mod)
    obs = obs[mask]
    mod = mod[mask]

    obs_std = np.std(obs)
    mod_std = np.std(mod)
    correlation = np.corrcoef(obs, mod)[0, 1]
    rmse = np.sqrt(np.mean((obs - mod) ** 2))
    bias = np.mean(mod - obs)
    wss = 1 - (np.sum((obs - mod)**2) / np.sum((np.abs(obs - np.mean(obs)) + np.abs(mod - np.mean(obs)))**2))

    return {
        'obs_std': obs_std,
        'mod_std': mod_std,
        'R': correlation,
        'RMSE': rmse,
        'Bias': bias,
        'WSS': wss
    }

# Compute and save skill statistics to CSV
skill_rows = []

for col in obs_cols + ['TOT']:
    obs_log = df_station[f'log10_{col}']
    model_col = f"EWE-{col}" if col != 'TOT' else 'EWE-TOT'
    mod_log = df_station[f'log10_{model_col}']
    stats = compute_skill_statistics(obs_log, mod_log)
    stats['Group'] = col
    skill_rows.append(stats)

skill_df = pd.DataFrame(skill_rows)
skill_df = skill_df.sort_values(by='WSS', ascending=False)
skill_df.to_csv(f"{path_evalout}/model_skill_stats_log10.csv", index=False)

print("\n=== Skill Statistics Saved to model_skill_stats_log10.csv ===")
print(skill_df)

print('done')