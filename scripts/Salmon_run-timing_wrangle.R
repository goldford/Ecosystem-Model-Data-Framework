# process lifecycle data 
require(readr)
require(tidyverse)
require(ggplot2)

timing <- read_csv("data/3Life_cycle_timing_by_CU.csv")

# one column per species per timing-type 
# e.g. juvenile marine entry & return 
# date column formatted as 2017-08
# within a given year the values should sum to 1 - kind of like a partial probabilty per month of the timing 
# chum, sokeye, coho, chinook, pink 
# chinook differentiated to stream type and ocean type 

table(timing$species)

#CK = chinook
#CM = chum
#CO = coho
#PKE = Pink even
#PKO = Pink odd 
#SEL 
#SER 
#SH = SH

timing %>% filter(species == "SEL")

# want to look at fm_peak & rt_peak

# species summaries 
table(timing$region)

# check species variable
str(timing$species)
timing$species <- factor(timing$species) 
levels(timing$species)

mean_timings <- timing %>%
  filter(region %in% c("fraser","vancouver island and mainland inlets")) %>%
  group_by(species) %>%
  summarize(mean_fm = mean(fm_peak,na.rm=T),
         sd_fm = sd(fm_peak,na.rm=T),
         mean_rt = mean(rt_peak,na.rm=T),
         sd_rt = sd(rt_peak,na.rm=T))

timing %>%
  filter(region %in% c("fraser","vancouver island and mainland inlets")) %>%
  filter(species == "CM")

mean_timings$mean_fm_month <- format(as.Date(mean_timings$mean_fm, origin = "2023-01-01"), "%m")
mean_timings$mean_rt_month <- format(as.Date(mean_timings$mean_rt, origin = "2023-01-01"), "%m")

# Convert all days to months 

monthly_timing <- timing %>%
  filter(region %in% c("fraser","vancouver island and mainland inlets")) %>%
  mutate(oe_start = as.numeric(format(as.Date(oe_start, origin = "2023-01-01"), "%m")),
         oe_peak = as.numeric(format(as.Date(oe_peak, origin = "2023-01-01"), "%m")),
         oe_end = as.numeric(format(as.Date(oe_end, origin = "2023-01-01"), "%m")),
         rt_start = as.numeric(format(as.Date(rt_start, origin = "2023-01-01"), "%m")),
         rt_peak = as.numeric(format(as.Date(rt_peak, origin = "2023-01-01"), "%m")),
         rt_end = as.numeric(format(as.Date(rt_end, origin = "2023-01-01"), "%m"))) %>%
  select(region, species, culabel, oe_start, oe_peak, oe_end, rt_start, rt_peak, rt_end) %>%
  filter(species != "SH") %>%
  mutate(species = case_when(
    species %in% c("SEL", "SER") ~ "Sockeye",
    TRUE ~ species))

month_summary <- monthly_timing %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))


# Partial probabilities --------------------------------------

# Combine timing data across all stocks for each species and run type
monthly_probs <- monthly_timing %>%
  pivot_longer(cols = c(oe_start, oe_peak, oe_end, rt_start, rt_peak, rt_end),
               names_to = "timing_type", values_to = "month") %>%
  mutate(type = ifelse(grepl("^oe", timing_type), "oe", "rt")) %>%
  group_by(species, type) %>%
  summarise(
    mu = mean(month, na.rm = TRUE),
    sigma = sd(month, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    month_probs = list({
      months <- 1:12
      probs <- dnorm(months, mean = mu, sd = sigma)
      probs / sum(probs)  # Normalize to sum to 1
    })
  ) %>%
  ungroup() %>%
  unnest_wider(month_probs, names_sep = "_") %>%
  pivot_longer(cols = starts_with("month_probs_"),
               names_prefix = "month_probs_",
               names_to = "month",
               values_to = "probability") %>%
  mutate(month = as.integer(month))

# plot distributions 
ggplot(monthly_probs, aes(x = month, y = probability, color = species, linetype = type)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    #title = "Monthly Probability Distributions by Species and Run Type",
    x = "Month",
    y = "Probability",
    color = "Species",
    linetype = "Run Type"
  ) +
  #scale_color_manual(labels = c("CK" = "Chinook", "CM" = "Chum", "CO" = "Coho", "PKE" = "Pink (even)", "PKO" = "Pink (odd)")) +
  #scale_linetype_manual(labels = c("oe" = "Ocean Entry", "rt" = "Return")) + 
  theme_minimal()

#write_csv(monthly_probs, "data/migration_probabilities.csv")



























