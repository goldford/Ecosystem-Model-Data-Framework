#resident KW

import numpy as np

modelled_area = 11274
mass_mt_min = 2.3
mass_mt_max = 3.1 # from life tables in 2008 from Ford et al 2010, also couture 2024

######################### SRKW ###########################
n_srkw_min = 60
n_srkw_max = 80
p_yr_srkw_min = 0.2
p_yr_srkw_max = 0.6
p_area_srkw_min = 0.25
p_area_srkw_max = 0.5

x1_min, x1_max = mass_mt_min, mass_mt_max
x2_min, x2_max = n_srkw_min, n_srkw_max
x3_min, x3_max = p_yr_srkw_min, p_yr_srkw_max
x4_min, x4_max = p_area_srkw_min, p_area_srkw_max

# Number of simulations
n_simulations = 10000

# Sample from uniform distribution within the given ranges
x1_samples = np.random.uniform(x1_min, x1_max, n_simulations)
x2_samples = np.random.uniform(x2_min, x2_max, n_simulations)
x3_samples = np.random.uniform(x3_min, x3_max, n_simulations)
x4_samples = np.random.uniform(x4_min, x4_max, n_simulations)

y_samples = x1_samples * x2_samples * x3_samples * x4_samples

# Calculate the 95% confidence interval
ci_lower_B = np.percentile(y_samples, 2.5)
ci_upper_B = np.percentile(y_samples, 97.5)

ci_lower_B_SRKW = ci_lower_B / modelled_area
ci_upper_B_SRKW = ci_upper_B / modelled_area

# Output the results
print(f"95% Confidence Interval for SRKW: [{ci_lower_B_SRKW}, {ci_upper_B_SRKW}]")

######################### NRKW ###########################
mass_mt_min = 1.700
mass_mt_max = 3.0
n_srkw_min = 90
n_srkw_max = 140
p_present_min = 0.01
p_present_max = 0.04

x1_min, x1_max = mass_mt_min, mass_mt_max
x2_min, x2_max = n_srkw_min, n_srkw_max
x3_min, x3_max = p_present_min, p_present_max

# Number of simulations
n_simulations = 10000

# Sample from uniform distribution within the given ranges
x1_samples = np.random.uniform(x1_min, x1_max, n_simulations)
x2_samples = np.random.uniform(x2_min, x2_max, n_simulations)
x3_samples = np.random.uniform(x3_min, x3_max, n_simulations)

y_samples = x1_samples * x2_samples * x3_samples

# Calculate the 95% confidence interval
ci_lower_B = np.percentile(y_samples, 2.5)
ci_upper_B = np.percentile(y_samples, 97.5)

ci_lower_B_NRKW = ci_lower_B / modelled_area
ci_upper_B_NRKW = ci_upper_B / modelled_area

# Output the results
print(f"95% Confidence Interval for NRKW: [{ci_lower_B_NRKW}, {ci_upper_B_NRKW}]")

tot_kw_min = ci_lower_B_SRKW + ci_lower_B_NRKW
tot_kw_max = ci_upper_B_SRKW + ci_upper_B_NRKW
print(f"95% Confidence Interval for Both: [{tot_kw_min}, {tot_kw_max}]")

mean = (tot_kw_min + tot_kw_max)/ 2
print(f"Expected B, best est, total: [{mean}]")
print(f"Min and max as % of best est: -{(mean-tot_kw_min) / mean} to +{(tot_kw_max-mean) / mean}")
