# sockeye salmon script
# by: G Oldford
# Nov 2024
# simple monte carlo to estimate CI of B, G, etc
# for sources etc see params and write up in tech report
#
# Eqns
# main equation:
# B_juve_mean = N_0 * w_0 * (e^(G-M) - 1) / (G-M)
# M = M_tot_smolt_to_offshore
# N_0 = n_adults / M_tot_smolt_to_adult
# w_0 = a_lw_juve * l_smlt_entry ^ b_lw_juve
# G = LN(w_juve_offshoremig / w_0)
# w_juve_offshoremig = a_lw_juve * l_juve_offshoremig ^ b_lw_juve
# assume a constant: a_lw_juve, b_lw_juve
# assume lognormal dist: l_juve_offshoremig, l_smlt_entry, n_adults
# assume normal dist:  M, M_tot_smolt_to_adult


import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg') # https://matplotlib.org/stable/users/explain/figure/backends.html


modelled_area = 11274
t_year = 40 / 365 # time spent in area

# adult fish params
n_adults_mean = 4380000 # based on total run size
n_adults_min = math.exp(math.log(n_adults_mean) * 0.9) # arbitrary guess - 13 mil
n_adults_max = math.exp(math.log(n_adults_mean) * 1.1)  #


# Recalculate the log-space standard deviation
n_adults_min_log = np.log(n_adults_min)
n_adults_max_log = np.log(n_adults_max)
n_adults_std_log = (n_adults_max_log - n_adults_min_log) / (2 * 1.96)
w_adults_mean = 2.84
w_adults_min = 2.78
w_adults_max = 2.9
w_adults_std_log = (np.log(n_adults_max) - np.log(n_adults_min)) / (2 * 1.96)

M_tot_smolt_to_adult_mean = -2.7 # not annualised!
M_tot_smolt_to_adult_min = -2.2
M_tot_smolt_to_adult_max = -3.2
M_tot_smolt_to_adult_std = (M_tot_smolt_to_adult_max - M_tot_smolt_to_adult_min) / (2 * 1.96)

# juve fish params
a_lw_juve = 2.4E-06  # healey, g mm-1
b_lw_juve = 3.32

l_smlt_entry_mean = 100 # mm; length at time of estuary entry
# l_smlt_entry_min = 28 # # did this indirectly, below
# l_smlt_entry_max = 60
l_smlt_entry_mean_log = math.log(l_smlt_entry_mean)
l_smlt_entry_min_log = l_smlt_entry_mean_log * 0.95
l_smlt_entry_max_log = l_smlt_entry_mean_log * 1.05
l_smlt_entry_min = math.exp(l_smlt_entry_min_log)
l_smlt_entry_max = math.exp(l_smlt_entry_max_log)
l_smlt_entry_std_log = (l_smlt_entry_max_log - l_smlt_entry_min_log) / (2 * 1.96) # assume normal dist
#
l_juve_offshoremig_mean = 132 # mm; length at time of offshore mig
l_juve_offshoremig_mean_log = math.log(l_juve_offshoremig_mean)
l_juve_offshoremig_min_log = l_juve_offshoremig_mean_log * 0.98
l_juve_offshoremig_max_log = l_juve_offshoremig_mean_log * 1.02
l_juve_offshoremig_std_log = (l_juve_offshoremig_max_log - l_juve_offshoremig_min_log) / (2 * 1.96) # assumed lognormal
l_juve_offshoremig_min = math.exp(l_juve_offshoremig_min_log)
l_juve_offshoremig_max = math.exp(l_juve_offshoremig_max_log)
# w_smlt_entry_mean = a_lw_juve * l_smlt_entry_mean ** b_lw_juve
# w_smlt_entry_min = a_lw_juve * l_smlt_entry_min ** b_lw_juve
# w_smlt_entry_max = a_lw_juve * l_smlt_entry_max ** b_lw_juve
w_juve_offshoremig_mean = a_lw_juve * l_juve_offshoremig_mean ** b_lw_juve
w_juve_offshoremig_min = a_lw_juve * l_juve_offshoremig_min ** b_lw_juve
w_juve_offshoremig_max = a_lw_juve * l_juve_offshoremig_max ** b_lw_juve
w_juve_offshoremig_std = (w_juve_offshoremig_max - w_juve_offshoremig_min) / (2 * 1.96)

m_daily_smolt_to_offshore_mean = 0.023 # daily mort %
t_entry_to_offshore = 40 # days; time from estuary entry to offshore migration
M_tot_smolt_to_offshore_mean = math.log(1-m_daily_smolt_to_offshore_mean) * t_entry_to_offshore
M_tot_smolt_to_offshore_min = M_tot_smolt_to_offshore_mean * 0.6
M_tot_smolt_to_offshore_max = M_tot_smolt_to_offshore_mean * 1.4
# ensure normally distributed in exponential space
m_daily_smolt_to_offshore_min = 1 - math.exp(M_tot_smolt_to_offshore_min / 90) # just to check
m_daily_smolt_to_offshore_max = 1 - math.exp(M_tot_smolt_to_offshore_max / 90) # just to check
M_tot_smolt_to_offshore_std = (M_tot_smolt_to_offshore_max - M_tot_smolt_to_offshore_min) / (2 * 1.96)
check = (M_tot_smolt_to_offshore_max + M_tot_smolt_to_adult_min) / 2


# Monte Carlo simulation
n_simulations = 50000

B_means = []
n_adults_all = []
w_adults_all = []
B_adults_all = [] # B of adults not used for calcs but saved for the model parmaeterisation
N_juve_all = []
l_smlts_all = []
w_smlts_all = []
l_juve_offshore_all = []
w_juve_offshore_all = []
G_all = []
M_all = []
G_delta_M_all = []


for _ in range(n_simulations):

    # Sample parameters
    n_adults = np.random.lognormal(np.log(n_adults_mean), n_adults_std_log)
    w_adults = np.random.lognormal(np.log(w_adults_mean), w_adults_std_log)
    B_adults = n_adults * w_adults
    n_adults_all.append(n_adults)
    w_adults_all.append(w_adults)
    B_adults_all.append(B_adults)

    l_smlt_entry = np.random.lognormal(l_smlt_entry_mean_log, l_smlt_entry_std_log)
    l_juve_offshoremig = np.random.lognormal(l_juve_offshoremig_mean_log, l_juve_offshoremig_std_log)
    M_tot_smolt_to_adult = np.random.normal(M_tot_smolt_to_adult_mean, abs(M_tot_smolt_to_adult_std))
    M_tot_smolt_to_offshore = np.random.normal(M_tot_smolt_to_offshore_mean, abs(M_tot_smolt_to_offshore_std))

    # Calculate juvenile weights
    w_smlt_entry = a_lw_juve * l_smlt_entry ** b_lw_juve
    w_juve_offshoremig = a_lw_juve * l_juve_offshoremig ** b_lw_juve

    # growth
    G = np.log(w_juve_offshoremig / w_smlt_entry)
    # Init abundance
    N_0 = n_adults / np.exp(M_tot_smolt_to_adult)

    # B_mean calculation
    # changing sign of M because sign is already in M calc above
    B_juve_mean = N_0 * w_smlt_entry * (np.exp(G + M_tot_smolt_to_offshore) - 1) / (G + M_tot_smolt_to_offshore)

    l_smlts_all.append(l_smlt_entry)
    w_smlts_all.append(w_smlt_entry)
    l_juve_offshore_all.append(l_juve_offshoremig)
    w_juve_offshore_all.append(w_juve_offshoremig)
    N_juve_all.append(N_0)
    G_all.append(G)
    M_all.append(M_tot_smolt_to_offshore)
    G_delta_M_all.append(G + M_tot_smolt_to_offshore)

    B_means.append(B_juve_mean)

# Calculate mean and confidence intervals
# assumes we explored data and determined log-normal distribution
N_mean_estimate = np.exp(np.mean(np.log(np.array(N_juve_all))))
N_mean_ci_lower = np.exp(np.percentile(np.log(np.array(N_juve_all)), 2.5))
N_mean_ci_upper = np.exp(np.percentile(np.log(np.array(N_juve_all)), 97.5))

print(f"Estimated number of juveniles: {N_mean_estimate:.2f}")
print(f"95% Confidence Interval: [{N_mean_ci_lower:.2f}, {N_mean_ci_upper:.2f}]")
print("")

# Calculate mean and confidence intervals
B_mean_estimate = np.mean(B_means) * 1/1000 * 1/1000 # convert to mt
B_mean_ci_lower = np.percentile(B_means, 2.5) * 1/1000 * 1/1000
B_mean_ci_upper = np.percentile(B_means, 97.5) * 1/1000 * 1/1000

print(f"Estimated B_mean juveniles: {B_mean_estimate:.2f} mt")
print(f"95% Confidence Interval: [{B_mean_ci_lower:.2f} mt, {B_mean_ci_upper:.2f}] mt")
print("")

B_dens_estimate = B_mean_estimate / modelled_area
B_dens_ci_lower = B_mean_ci_lower / modelled_area
B_dens_ci_upper = B_mean_ci_upper / modelled_area

print(f"Estimated juve B density for model: {B_dens_estimate:.4f} mt km-2")
print(f"95% Confidence Interval: [{B_dens_ci_lower:.4f} mt km-2, {B_dens_ci_upper:.4f}] mt km-2")
print("")

# adjusting for seasonal
B_dens_estimate_seas = B_dens_estimate * t_year
B_dens_ci_lower_seas = B_dens_ci_lower * t_year
B_dens_ci_upper_seas = B_dens_ci_upper * t_year
print(f"Estimated juve B density for model, seasonal adjusted: {B_dens_estimate_seas:.5f} mt km-2")
print(f"95% Confidence Interval: [{B_dens_ci_lower_seas:.5f} mt km-2, {B_dens_ci_upper_seas:.5f}] mt km-2")
print("")

#####################
##### ADULT FISH ####
n_adults_est = np.exp(np.mean(np.log(n_adults_all)))
n_adults_ci_lower = np.exp(np.percentile(np.log(n_adults_all), 2.5))
n_adults_ci_upper = np.exp(np.percentile(np.log(n_adults_all), 97.5))
print(f"Estimated number of adult fish: {n_adults_est:.5f} ")
print(f"95% Confidence Interval: [{n_adults_ci_lower:.5f}, {n_adults_ci_upper:.5f}]")
print("")

w_adults_est = np.mean(w_adults_all)
w_adults_ci_lower = np.percentile(w_adults_all, 2.5)
w_adults_ci_upper = np.percentile(w_adults_all, 97.5)
print(f"Estimated weight of adult fish: {w_adults_est:.5f} kg")
print(f"95% Confidence Interval: [{w_adults_ci_lower:.5f} kg, {w_adults_ci_upper:.5f}] kg")

B_adults_all = np.array(B_adults_all)
B_adults_est = np.exp(np.mean(np.log(B_adults_all)))
B_adults_ci_lower = np.exp(np.percentile(np.log(B_adults_all), 2.5))
B_adults_ci_upper = np.exp(np.percentile(np.log(B_adults_all), 97.5))
print(f"Estimated total B of adult fish: {B_adults_est:.5f} kg")
print(f"95% Confidence Interval: [{B_adults_ci_lower:.5f} kg, {B_adults_ci_upper:.5f}] kg")
print("")


# B adults total mt
plt.figure(figsize=(10, 6))
plt.hist(B_adults_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of adult B (kg) from Monte Carlo Simulations")
plt.xlabel("kg")
plt.ylabel("Frequency")
plt.show()

# adult B den
Bden_adults_all = np.array(B_adults_all) * 1/1000 * 2/12 * 1/modelled_area
Bden_adults_all_logmean = np.mean(np.log(Bden_adults_all))
log_std = np.std(np.log(Bden_adults_all))
lower_CI_log = Bden_adults_all_logmean - 1.96 * log_std
upper_CI_log = Bden_adults_all_logmean + 1.96 * log_std

plt.figure(figsize=(10, 6))
plt.hist(Bden_adults_all, bins=50, color='skyblue', edgecolor='black')
plt.axvline(math.exp(Bden_adults_all_logmean), color='red')
plt.axvline(math.exp(lower_CI_log), color='red', linestyle='--')
plt.axvline(math.exp(upper_CI_log), color='red', linestyle='--')
plt.title("Histogram of adult B den from Monte Carlo Simulations")
plt.xlabel(r"mt km$^-$" + r"$^2$")
plt.ylabel("Frequency")
plt.show()

Bden_adults_est = B_adults_est * 1/1000 * 2/12 * 1/modelled_area#  convert to mt and seasonality
print(f"Estimated adult biomass density, seasonal adjusted: {Bden_adults_est:.5f} mt km^-2")
print(f"95% Confidence Interval: [{math.exp(lower_CI_log):.5f} mt km^-2, {math.exp(upper_CI_log):.5f}] mt km^-2")
print("")


###### JUVENILES ######
# # l at entry
# plt.figure(figsize=(10, 6))
# plt.hist(l_smlts_all, bins=50, color='skyblue', edgecolor='black')
# plt.title("Histogram of smolt lengths at estuary entry from Monte Carlo Simulations")
# plt.xlabel("length (mm)")
# plt.ylabel("Frequency")
# plt.show()
#
# # l at exit
# plt.figure(figsize=(10, 6))
# plt.hist(l_juve_offshore_all, bins=50, color='skyblue', edgecolor='black')
# plt.title("Histogram of smolt lengths at open ocean migration from Monte Carlo Simulations")
# plt.xlabel("length (mm)")
# plt.ylabel("Frequency")
# plt.show()
#
# # w at entry
plt.figure(figsize=(10, 6))
plt.hist(w_smlts_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of smolt weights at estuary entry from Monte Carlo Simulations")
plt.xlabel("weight (g)")
plt.ylabel("Frequency")
plt.show()

w_smlts_all = np.array(w_smlts_all) # assume log-normal
log_w_smlts_all = np.log(w_smlts_all)
log_mean = np.mean(log_w_smlts_all)
log_std = np.std(log_w_smlts_all)
lower_limit_log = log_mean - 2 * log_std
upper_limit_log = log_mean + 2 * log_std
lower_limit = np.exp(lower_limit_log)
upper_limit = np.exp(upper_limit_log)
# Filter
# w_smlts_means_truncated = [value for value in B_mean_estimate if lower_limit <= value <= upper_limit]

# since distribution is log-normal. Take the log-normal mean and CI's
lower_CI_log = log_mean - 1.96 * log_std
upper_CI_log = log_mean + 1.96 * log_std
print(f"Estimated w of smolts at o-e: {math.exp(log_mean):.2f} g")
print(f"95% Confidence Interval: [{math.exp(lower_CI_log):.2f} g, {math.exp(upper_CI_log):.2f}] g")
print("")



# w at exit
w_juve_offshore_all = np.array(w_juve_offshore_all) # assume log-normal
log_w_juve_offshore_all = np.log(w_juve_offshore_all)
log_mean = np.mean(log_w_juve_offshore_all)
log_std = np.std(log_w_juve_offshore_all)
lower_limit_log = log_mean - 2 * log_std
upper_limit_log = log_mean + 2 * log_std
lower_limit = np.exp(lower_limit_log)
upper_limit = np.exp(upper_limit_log)
# Filter
# w_smlts_means_truncated = [value for value in B_mean_estimate if lower_limit <= value <= upper_limit]

# since distribution is log-normal. Take the log-normal mean and CI's
lower_CI_log = log_mean - 1.96 * log_std
upper_CI_log = log_mean + 1.96 * log_std
print(f"Estimated w of smolts at offshore migration: {math.exp(log_mean):.2f} g")
print(f"95% Confidence Interval: [{math.exp(lower_CI_log):.2f} g, {math.exp(upper_CI_log):.2f}] g")
print("")


plt.figure(figsize=(10, 6))
plt.hist(w_juve_offshore_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of smolt weights at open ocean migration from Monte Carlo Simulations")
plt.xlabel("weight (g)")
plt.ylabel("Frequency")
plt.show()

#
# n juves
plt.figure(figsize=(10, 6))
plt.hist(N_juve_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of number of fish entering estuary entry from Monte Carlo Simulations")
plt.xlabel("n")
plt.ylabel("Frequency")
plt.show()

n_juve_est = np.mean(N_juve_all)
n_juve_est_med = np.median(N_juve_all)
n_juve_ci_lower = np.percentile(N_juve_all, 2.5)
n_juve_ci_upper = np.percentile(N_juve_all, 97.5)
print(f"Estimated number of juve fish (mean): {n_juve_est:.5f} ")
print(f"Estimated number of juve fish (median): {n_juve_est_med:.5f} ")
print(f"95% Confidence Interval: [{n_juve_ci_lower:.5f}, {n_juve_ci_upper:.5f}]")

# print(f"Estimated N juves: {math.exp(log_mean):.2f} mt")
# print(f"95% Confidence Interval: [{math.exp(lower_CI_log):.2f} mt, {math.exp(upper_CI_log):.2f}] mt")


# M
plt.figure(figsize=(10, 6))
plt.hist(G_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of instantaneous mortality, M, over early marine period from Monte Carlo Simulations")
plt.xlabel("/ mo.")
plt.ylabel("Frequency")
M_mean = np.mean(M_all)
plt.axvline(M_mean, color='red')
M_std = np.std(M_all)
M_upper_CI = M_mean + 1.96 * M_std
M_lower_CI = M_mean - 1.96 * M_std
plt.axvline(M_lower_CI, color='red', linestyle='--')
plt.axvline(M_upper_CI, color='red', linestyle='--')
plt.show()


# G
plt.figure(figsize=(10, 6))
plt.hist(G_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of instantaneous growth, G, over early marine period from Monte Carlo Simulations")
plt.xlabel("/ mo.")
plt.ylabel("Frequency")
G_mean = np.mean(G_all)
plt.axvline(G_mean, color='red')
G_std = np.std(G_all)
G_upper_CI = G_mean + 1.96 * G_std
G_lower_CI = G_mean - 1.96 * G_std
plt.axvline(G_lower_CI, color='red', linestyle='--')
plt.axvline(G_upper_CI, color='red', linestyle='--')
plt.show()

plt.figure(figsize=(10, 6))
plt.hist(G_delta_M_all, bins=50, color='skyblue', edgecolor='black')
plt.title("Histogram of instantaneous growth, G, minus instantaneous M over early marine period from Monte Carlo Simulations")
plt.xlabel("/ mo.")
plt.ylabel("Frequency")
plt.show()

g_day = (math.exp(G_mean / 30) - 1) * 100
g_day_upper_CI = (math.exp(G_upper_CI / 30) - 1) * 100
g_day_lower_CI = (math.exp(G_lower_CI / 30) - 1) * 100

print(f"Estimated g daily % of JUVENILES: {g_day:.2f}")
print(f"95% Confidence Interval: [{math.exp(g_day_upper_CI/3):.2f}, {math.exp(g_day_lower_CI/3):.2f}]")
print("")

print(f"Estimated G yr-1 of JUVENILES: {G_mean*12:.2f}")
print(f"95% Confidence Interval: [{G_lower_CI*12:.2f}, {G_upper_CI*12:.2f}]")
print("")


# final B
B_mean_estimate = np.array(B_means) * 1/1000 * 1/1000 # convert to mt
log_B_means = np.log(B_mean_estimate)
log_mean = np.mean(log_B_means)
log_std = np.std(log_B_means)
lower_limit_log = log_mean - 2 * log_std
upper_limit_log = log_mean + 2 * log_std
lower_limit = np.exp(lower_limit_log)
upper_limit = np.exp(upper_limit_log)
# Filter
B_means_truncated = [value for value in B_mean_estimate if lower_limit <= value <= upper_limit]

# since distribution is obviously log-normal. Take the log-normal mean and CI's
lower_CI_log = log_mean - 1.96 * log_std
upper_CI_log = log_mean + 1.96 * log_std

plt.figure(figsize=(10, 6))
plt.hist(B_means_truncated, bins=50, color='skyblue', edgecolor='black')
# plt.xlim(0,1300)
plt.axvline(math.exp(log_mean), color='red')
plt.axvline(math.exp(lower_CI_log), color='red', linestyle='--')
plt.axvline(math.exp(upper_CI_log), color='red', linestyle='--')
plt.title("Histogram of Biomass of juvenile pink salmon entering estuary entry from Monte Carlo Simulations")
plt.xlabel("mt")
plt.ylabel("Frequency")
plt.show()



print(f"Estimated B_mean of JUVENILES: {math.exp(log_mean):.2f} mt")
print(f"95% Confidence Interval: [{math.exp(lower_CI_log):.2f} mt, {math.exp(upper_CI_log):.2f}] mt")
print("")

B_dens_estimate = math.exp(log_mean) / modelled_area
B_dens_ci_lower = math.exp(lower_CI_log) / modelled_area
B_dens_ci_upper = math.exp(upper_CI_log) / modelled_area
B_dens_estimate_seas = B_dens_estimate * t_year
B_dens_ci_lower_seas = B_dens_ci_lower * t_year
B_dens_ci_upper_seas = B_dens_ci_upper * t_year

print(f"Estimated B density for model (Juveniles): {B_dens_estimate_seas:.5f} mt km-2")
print(f"95% Confidence Interval: [{B_dens_ci_lower_seas:.5f} mt km-2, {B_dens_ci_upper_seas:.5f}] mt km-2")
print("")