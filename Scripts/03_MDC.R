# A script to produce figure 3
# A: Robert Buchkowski & Alex Polussa

# Load packages:
library(lubridate)
library(lme4)
library(tidyverse)

library(patchwork)

# parallel
library(furrr)
library(progressr)

# Load in the necessary data:
source("Scripts/00_load_data.R")
source("Scripts/99_functions.R")


# Spatial MDD -------------------------------------------------------------

# Split up dataset by horizon and time point and calculate MDD

bwvr_min_1 <- bwvr_alex(Mgpremeas %>% filter(meas_num == 1), vars = "TC", unit_col = "nfi_plot")
bwvr_min_2 <- bwvr_alex(Mgpremeas %>% filter(meas_num == 2), vars = "TC", unit_col = "nfi_plot")
bwvr_org_1 <- bwvr_alex(Ogpremeas %>% filter(meas_num == 1), vars = "TC", unit_col = "nfi_plot")
bwvr_org_2 <- bwvr_alex(Ogpremeas %>% filter(meas_num == 2), vars = "TC", unit_col = "nfi_plot")

# Combine
MDD_spatial <-
  rbind(
    mdc_from_bwvr(bwvr_min_1) %>% mutate(soil = "Mineral", Time = 1),
    mdc_from_bwvr(bwvr_min_2) %>% mutate(soil = "Mineral", Time = 2),
    mdc_from_bwvr(bwvr_org_1) %>% mutate(soil = "Organic", Time = 1),
    mdc_from_bwvr(bwvr_org_2) %>% mutate(soil = "Organic", Time = 2))

# Plot (Spatial)
plot_spatial <- MDD_spatial %>%
  mutate(label = "MDD with Carbon stock spatial variation") %>%
  ggplot(aes(
    x = paste0(soil, "\n(Time - ", as.factor(Time-1), ")"),
    y = MDD_pop_alpha_power,
    fill = soil
  )) +
  geom_col(position = "dodge", color = "black") +
  facet_grid(~ label) +
  scale_fill_manual(name = "Horizon", values = c("grey", "brown") ) +
  labs(    x = NULL,    y = expression(MDD~(Mg~C~ha^-1))  ) +
  ylim(0, 5) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "white", color = "black", size = 0.9),
    strip.text = element_text(size = 12),
    plot.margin = margin(5, 5, 5, 5)
  )

plot_spatial


# Temporal MDD ------------------------------------------------------------


MDD_temporal <-   rbind( mdc_pop_full_yr(Mgpremeas %>% filter(meas_num != 3),
                                         unit_col   = "nfi_plot",
                                         time_col   = "meas_num",
                                         date_col   = "meas_date",
                                         value_col  = "TC",
                                         repl_col   = "pit_num"    ,
                                         alpha = 0.05,
                                         power = 0.80,
                                         horizon_years = 10)  %>% mutate(soil = "Mineral", Time = "Temporal"),
                         mdc_pop_full_yr(Ogpremeas %>% filter(meas_num != 3),
                                         unit_col   = "nfi_plot",
                                         time_col   = "meas_num",
                                         date_col   = "meas_date",
                                         value_col  = "TC",
                                         repl_col   = "pit_num",
                                         alpha = 0.05,
                                         power = 0.80,
                                         horizon_years = 10)  %>% mutate(soil = "Organic", Time = "Temporal")
)

# Plot (temporal)
plot_temp <- MDD_temporal %>%
  mutate(label = "MDD with Carbon Stock Change") %>%
  ggplot(aes(x = soil, y = MDD_horizon_a_p, fill = soil)) +
  geom_col(color = "black") +
  facet_grid(~ label) +
  scale_fill_manual(name = "Horizon",
                    values = c("grey", "brown")) +
  labs(x = NULL, y = expression(MDD~(Mg~C~ha^-1~10~yr^-1))) +
  ylim(0, 5) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "white", color = "black", size = 0.8),
    strip.text = element_text(size = 12),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

plot_temp

# MDD for spatial over 10 year horizon:
MDD_temporal %>% dplyr::select(MDD_horizon_a_p, soil)

# 1            4.08 Mineral
# 2            4.61 Organic

MDD_spatial %>% dplyr::select(MDD_pop_alpha_power, soil)

# 1                2.91 Mineral
# 2                3.16 Mineral
# 3                4.49 Organic
# 4                4.31 Organic


# Plot Both
plot_spatial + plot_temp

ggsave("Plots/03_MDC_temporal_spatial.png", width = 9, height = 4)

# Calculate the total change that could be detected:

# Area from: https://nfi.nfis.org/resources/general/summaries/t1/en/NFI/html/nfi_t4_for_area_en.html

forested_hectares = 369019.35*1000
4.08*forested_hectares/(10^9) #mineral_detect Pg
4.61*forested_hectares/(10^9) #organic_detect Pg


