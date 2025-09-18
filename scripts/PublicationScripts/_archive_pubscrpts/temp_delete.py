nutrient_decline_fraction = 0.60  # 60% decline in nutrients under productive conditions
biomass_increase_fraction = 0.80  # 80% increase in phytoplankton biomass under productive conditions

# Let N_eq be equilibrium nutrient concentration
# Let P_eq be equilibrium phytoplankton biomass

# Under productive conditions:
# Nutrients = 0.4 * N_eq
# Phytoplankton = 1.8 * P_eq

# Total nitrogen is conserved:
# N_eq + P_eq = 0.4 * N_eq + 1.8 * P_eq

# Rearrange to solve for N_eq / P_eq:
# N_eq - 0.4 * N_eq = 1.8 * P_eq - P_eq
# 0.6 * N_eq = 0.8 * P_eq
# N_eq / P_eq = 0.8 / 0.6 = 4/3
# FIX SO IT'S READING FROM FIRST TIME STEP
INITIAL_BIOMASS = 8.5  # g C m^-2, for initial phytoplankton biomass only
MAX_BIOMASS = 15.3 # under bloom conditions
# New input assumptions
INITIAL_NUTRIENT_CONC = 18  # umol N L^-1 annual avg
MIN_N_CONC = 7.2       # μmol N L^-1, observed minimum concentration


# Calculation of 'Free' Nutrients (THIS SHOULD MATCH ECOSIM SCENARIO)
B_MAX_INCREASE_PROP = (MAX_BIOMASS - INITIAL_BIOMASS) / INITIAL_BIOMASS
N_MAX_DECREASE_PROP = 1 - (MIN_N_CONC / INITIAL_NUTRIENT_CONC)
N_eq_to_P_eq = B_MAX_INCREASE_PROP / N_MAX_DECREASE_PROP

# Calculate total nitrogen in the system
# T = N_eq + P_eq
# Let P_eq = 1 unit (arbitrary), then N_eq = N_eq_to_P_eq * P_eq
P_eq = 1  # define as 1 unit for proportional calculation
N_eq = N_eq_to_P_eq * P_eq
T = N_eq + P_eq

# Calculate the proportion that is 'free' (i.e. as dissolved nutrients)
PROP_N_FREE = N_eq / T
print(f"Proportion of nitrogen that is free (dissolved) under equilibrium conditions: {PROP_N_FREE:.3f}")