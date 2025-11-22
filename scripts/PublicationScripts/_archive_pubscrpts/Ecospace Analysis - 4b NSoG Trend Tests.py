# Created by G Oldford, Sep 23, 2024
# Purpose: Evaluate Trends in NSoG from Ecospace Sims
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import mannkendall as mk
import xarray as xr
import datetime as dt

bloom_p = "..//..//data//evaluation//"
ecospace_bloom_timing_df = pd.read_csv(os.path.join(bloom_p, 'ecospace_bloom_timing_NSoG.csv'))

print(ecospace_bloom_timing_df)

# prep - need to convert to units etc that stats functions expect
print(ecospace_bloom_timing_df)
# ecospace_bloom_timing_df['Bloom Date'] = pd.to_datetime(print(ecospace_bloom_timing_df['Bloom Date']))


df = ecospace_bloom_timing_df
bd_col = 'Bloom Date'
df[bd_col] = pd.to_datetime(df[bd_col])

# Create the numpy array of datetime objects
dat_pytime = np.array([dt.datetime(year, month, day)
                       for year, month, day in zip(df[bd_col].dt.year.values,
                                                   df[bd_col].dt.month.values,
                                                   df[bd_col].dt.day.values)])

xarray_da = df.set_index(bd_col).to_xarray()
dat = xarray_da['Day of Year']

print(dat)
print(dat_pytime)

# evaluate trend. will need to check residuals for autocorr

def do_ss(dat, dat_pytime, resolution):
    # created by G Oldford Jan 2024
    # calling mks.sens_slope is not usually called directly in mk lib,
    # so this code was borrowed from mannkendall lib
    # Intercept calculated using Conover, W.J. (1980) method (from pymannkendall)
    #
    # inputs
    #   dat              - depth-averaged xarray with one variable and time dimension
    #   dat_pytime       - xarray or numpy array, the time counter var from dat_davg as py datetime, required by mk
    #   resolution       - the precision of the measurements, significant digits
    # returns
    #   list of dictionary elements, slope and upper conf. limit, lower conf, intercept - annual
    #
    #  note all slopes are converted to units of years but sen_slope returns units of seconds

    SS_slopes = []
    alpha_cl = 90

    dat_davg_np = np.asarray(dat)

    print(dat_davg_np)


    # this is not correct! should be slopes not data used as input
    # it is as implemented in mk lib but seems incorrect when compared to glibert p. 227
    t = mk.mkt.nb_tie(dat_davg_np, resolution)
    (s, n) = mk.mks.s_test(dat_davg_np, dat_pytime)
    k_var = mk.mkt.kendall_var(dat_davg_np, t, n)
    z = mk.mks.std_normal_var(s, k_var)
    mk_out = mk.mks.sen_slope(dat_pytime, dat_davg_np,k_var,alpha_cl=alpha_cl) # units of second!
    mult = 3600 * 24 * 365.25
    slope_yr = mk_out[0] * mult; slope_yr_min = mk_out[1] * mult; slope_yr_max = mk_out[2] * mult
    inter_yr = np.nanmedian(dat_davg_np) - np.median(np.arange(len(dat_davg_np))[~np.isnan(dat_davg_np.flatten())]) * slope_yr
    # or median(x) - (n-1)/2 *slope
    SS_slopes.append({"annual": {"slope": slope_yr,
                                      "slope_LCL": slope_yr_min,
                                      "slope_UCL": slope_yr_max,
                                      "inter":inter_yr}})
    return SS_slopes


# compute Thiel-Sen aka Sens Slope
SS_slopes = do_ss(dat, dat_pytime, 0.01)

print(SS_slopes)

# alt way following modded code used in do_ss_seasonal in utilities of HOTSSea
SS_slopes_2 = []
dat_np = np.asarray(dat)
result, s, vari, z = mk.compute_mk_stat(dat_pytime, dat_np, 0.01,
                                                    alpha_mk=0.95, alpha_cl=0.9)
# this returns annualised slopes using the pyttime
inter_yr = np.nanmedian(dat_np) - np.median(
                np.arange(len(dat_np))[~np.isnan(dat_np.flatten())]) * result['slope']
SS_slopes_2.append({"annual": {"slope": result['slope'],
                               "lcl": result['lcl'],
                               "ucl": result['ucl'],
                               "p": result['p'],
                               "ss": result['ss'],
                               "inter": inter_yr}})

# the ucl and lcl differ... between the methods but slopes etc same
print(SS_slopes_2)



from scipy import stats
from scipy.stats import linregress
# can either use TS (sens slope) or use LR (start w/ LR)
# if resids normal, then may as well use LR (MK w/ SS not necessary)
obs_years = np.array([item.year for item in dat_pytime])
result_LR = linregress(obs_years, dat_np) # using 'result' returns st err of intercept
# LR_slope, LR_inter, LR_rval, LR_pval, LR_stderr = linregress(obs_years, dat_np)
residuals = dat_np - (result_LR.slope * obs_years + result_LR.intercept)
stat_sw, p_value_sw = stats.shapiro(residuals)  # test for normality (small p val = small prob of normal)
# _, p_value_whites, _, _ = het_white(residuals, exog=add_constant(dat_1seas_timenum))  # test for heterosk
# _, p_value_goldfield, _ = het_goldfeldquandt(residuals, add_constant(dat_1seas_timenum))  # test for heterosk
# note that low p-vals in Het-White and Goldfeld-Quandt tests suggests evidence against homoscedasticity





# Durbin-Watson test
import statsmodels.api as sm
dw_statistic = sm.stats.durbin_watson(residuals)
print(f"Durbin-Watson statistic: {dw_statistic}") # if 2 no autocorr, if 4 perfect negative, if 0 positive autocorr
# alternative, more robust approach:
# autocorrelation
# from statsmodels.tsa.stattools import acf
# alpha_acf = 0.05
# nlags = len(residuals) - 1
# acf_values = acf(residuals, adjusted=True, qstat=True, fft=False, alpha=alpha_acf, nlags=nlags,
#                              missing="conservative")
# AR1_coef = acf_values[0][1]
# AR1_confint= acf_values[1]
# AR1_qstat = acf_values[2]
# AR1_pval = acf_values[3] # low means unlikely random (thus autocorr)
# print(f"AR1_coeff: {AR1_coef}")
# print(f"AR1_confint: {AR1_confint}")
# print(f"AR1_qstat: {AR1_qstat}")
# print(f"AR1_pval: {AR1_pval}")

# Plotting the autocorrelation function (ACF)
sm.graphics.tsa.plot_acf(residuals, lags=5)
# plt.show()


# Two-sided inverse Students t-distribution (from linregress source examp)
# p - probability, df - degrees of freedom
from scipy.stats import t
res = result_LR
print(f"R-squared: {res.rvalue**2:.6f}")

# Calculate 95% confidence interval on slope and intercept:
# Two-sided inverse Students t-distribution
# p - probability, df - degrees of freedom
from scipy.stats import t
tinv = lambda p, df: abs(t.ppf(p/2, df))
ts = tinv(0.05, len(obs_years)-2)
print(f"slope (95%): {res.slope:.6f} +/- {ts*res.stderr:.6f}")
print(f"intercept (95%): {res.intercept:.6f}" f" +/- {ts*res.intercept_stderr:.6f}")
print(f"p-val non-zero slope: {res.pvalue:.6f}") # lower more confident
res = linregress(obs_years, dat_np, alternative='less')
print(f"p-val negative slope: {res.pvalue:.6f}") # lower more confident

# print(LR_slope, LR_inter, LR_rval, LR_pval)
print(p_value_sw, stat_sw)
print("done")

slope_ci = ts * res.stderr
intercept_ci = ts * res.intercept_stderr




# to do - different time and time mult if lr
#######################################
################# PLOT ################
yr_strt = 1980
yr_end = 2018
scenario = 'FULLKEY_SC51_5'
xlim_min = 1979.5
xlim_max = 2018.5
log_transform = True
min_y_tick = 30
fig_width = 10
fig_height = 6
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
df = ecospace_bloom_timing_df

ax.errorbar(df['Year'], df['Day of Year'], yerr=1.5,
            fmt='s', color='black',
            label='Ecospace Model',
            capsize=0, marker='o', elinewidth=1,
            markersize=4)
# ax.axhline(y=68, color='black', linestyle='--')
# ax.axhline(y=108, color='black', linestyle='--')
x_min, x_max = ax.get_xlim()

# suchy et al 2022 definition of bloom timing early, avg, late
# ax.fill_between([xlim_min, x_max], 68, 108, color='lightgrey', alpha=0.5, zorder=0)
ax.set_xlim([xlim_min, xlim_max])
y_min, y_max = ax.get_ylim()
ax.set_ylim([min_y_tick, y_max])
ax.set_xlabel('Year')
ax.set_ylabel('Day of Year')
# ax.set_title('Diatom Bloom Timing Comparison')


# Plot the regression line
x_vals = np.array(ax.get_xlim())
y_vals = res.intercept + res.slope * x_vals
ax.plot(x_vals, y_vals, linestyle='--',
        color='black', label='Trend line')


#########
### alt way
# Calculate the 95% confidence interval
x_vals = np.linspace(df['Year'].min()-0.5, df['Year'].max()+0.5, 100)
y_vals = res.intercept + res.slope * x_vals
n = len(obs_years)
mean_x = np.mean(obs_years)
t_stat = t.ppf(0.975, df=n - 2)  # 95% confidence level
SE = np.sqrt(np.sum((dat_np - (res.intercept + res.slope * obs_years))**2) / (n - 2))
SE_fit = SE * np.sqrt(1/n + (x_vals - mean_x)**2 / np.sum((obs_years - mean_x)**2))
y_upper_alt = y_vals + t_stat * SE_fit
y_lower_alt = y_vals - t_stat * SE_fit


y_lower = (res.intercept - intercept_ci) + (res.slope - slope_ci) * x_vals
y_upper = (res.intercept + intercept_ci) + (res.slope + slope_ci) * x_vals
# ax.fill_between(x_vals, y_lower, y_upper, color='red', alpha=0.2, label='95% CI')
ax.fill_between(x_vals, y_lower_alt, y_upper_alt, color='grey', alpha=0.2, label='95% CI')

ax.legend()

#
# # Create a secondary y-axis
# ax2 = ax.twinx()
#
# # Converting the trend line to secondary y-axis
# y_min, y_max = ax.get_ylim()
# trend_min = res.slope * x_min + res.intercept
# trend_max = res.slope * x_max + res.intercept
# ax2.set_ylim((trend_min, trend_max) if log_transform else (trend_min, trend_max))
# ax2.set_ylabel('Trend value', color='red')
# ax2.tick_params(axis='y', labelcolor='red')
#
# # Adding a dashed grey horizontal line at the zero level for the secondary y-axis
# ax2.axhline(y=0, color='grey', linestyle='--', linewidth=1)



if log_transform:
    trnsfrm = "logDIA"
else:
    trnsfrm = ""

plt.title(scenario + " " + trnsfrm + " NSOG")
# plt.savefig('..//..//figs//' + 'BloomTiming_Ecospace_vs_Suchy_2003_2016_altmethod_' + scenario + '_' + trnsfrm + '.png')
plt.savefig('..//..//figs//' + 'BloomTiming_Ecospace_NSOG_wTrend_' + str(yr_strt) + '_' + str(yr_end) + '_suchymethod_' + scenario + '_' + trnsfrm + '.png')
plt.show()


