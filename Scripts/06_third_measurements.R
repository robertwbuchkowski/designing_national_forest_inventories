# A script to produce figure 6 showing within and between site variation
# A: Robert Buchkowski & Alex Polussa

# Load packages:
library(lubridate)
library(lme4)
library(tidyverse)

# Load in the necessary data:
source("Scripts/00_load_data.R")
source("Scripts/99_functions.R")

# List of third measurement plots:

third_plots = gpremeas %>%
  filter(meas_num == 3) %>%
  select(nfi_plot, Horizon) %>%
  distinct()


# 1191211 1349216

# Create a data frame of change:

gpremeas_calc_all = gpremeas %>%
  right_join(
    third_plots,by = join_by(nfi_plot, Horizon)
  ) %>%
  filter(meas_num < 3) %>%
  group_by(nfi_plot, pit_num, Horizon) %>%
  arrange(meas_date, .by_group = TRUE) %>%
  summarise(
    first_date = first(meas_date),
    last_date = last(meas_date),
    first_TC = first(TC),
    last_TC = last(TC),
    first_perC = first(perC),
    last_perC = last(perC),
    first_BD = first(BD),
    last_BD = last(BD),
    first_CF = first(CF),
    last_CF = last(CF),
    years_diff = as.numeric(difftime(last_date, first_date, units = "days")) / 365.25,
    TC_change_per_year = (last_TC - first_TC) / years_diff,
    perC_change_per_year = (last_perC - first_perC) / years_diff,
    BD_change_per_year = (last_BD - first_BD) / years_diff,
    CF_change_per_year = (last_CF - first_CF) / years_diff,
    .groups = "drop"
  ) %>%
  filter(years_diff > 0)


gpremeas_calc_all_2 = gpremeas %>%
  right_join(
    third_plots,by = join_by(nfi_plot, Horizon)
  ) %>%
  filter(meas_num > 1) %>%
  group_by(nfi_plot, pit_num, Horizon) %>%
  arrange(meas_date, .by_group = TRUE) %>%
  summarise(
    first_date = first(meas_date),
    last_date = last(meas_date),
    first_TC = first(TC),
    last_TC = last(TC),
    first_perC = first(perC),
    last_perC = last(perC),
    first_BD = first(BD),
    last_BD = last(BD),
    first_CF = first(CF),
    last_CF = last(CF),
    years_diff = as.numeric(difftime(last_date, first_date, units = "days")) / 365.25,
    TC_change_per_year = (last_TC - first_TC) / years_diff,
    perC_change_per_year = (last_perC - first_perC) / years_diff,
    BD_change_per_year = (last_BD - first_BD) / years_diff,
    CF_change_per_year = (last_CF - first_CF) / years_diff,
    .groups = "drop"
  ) %>%
  filter(years_diff > 0)

gpremeas_calc_all_all = gpremeas_calc_all %>%
  mutate(Change = "First") %>%
  bind_rows(
    gpremeas_calc_all_2 %>%
      mutate(Change = "Second")
  )


p1 = gpremeas_calc_all_all %>%
  select(nfi_plot, Horizon, Change, contains("change_per")) %>%
  pivot_longer(contains("change_per")) %>%
  group_by(nfi_plot, Horizon, Change, name) %>%
  summarize(value = mean(value)) %>%
  pivot_wider(names_from = Change, values_from = value) %>%
  mutate(name = case_match(name,
                           "BD_change_per_year"~"Change in Fine Fraction Bulk Density",
                           "CF_change_per_year"~"Change in Coarse Fraction Bulk Density",
                           "perC_change_per_year"~"Change in Percent Carbon",
                           "TC_change_per_year"~"Change in Carbon stock")) %>%
  ggplot(aes(x = First, y = Second, color = Horizon)) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point() + facet_wrap(Horizon~name, scales = "free", nrow = 4, ncol = 2,dir = "v") + theme_classic() + xlab("First Interval Change (T0 to T1)") + ylab("Second Interval Change (T1 to T2)") + scale_color_manual(name = "Horizon", values = c("grey", "brown"), guide = "none")

# Two NFI plots are dropped from this figure because measurements at time 1 and time 2 are from different microplots, so the paired comparison being conducted is not possible.

png("Plots/Figure6.png", width = 6, height = 6, units = "in", res = 600)
p1
dev.off()
