# A script to produce figure 4 sampling effort
# A: Robert Buchkowski & Alex Polussa

# Load packages:
library(lubridate)
library(lme4)
library(tidyverse)

library(patchwork)

# parallel
library(purrr)
library(furrr)
library(lubridate)

# Load in the necessary data:
source("Scripts/00_load_data.R")
source("Scripts/99_functions.R")

# only mineral horizon
mineral_df <- gpremeas %>% filter(Horizon == "Mineral" & meas_num !=3)


# Simulation -------------------------------------------------

# We need within and between spatial variation of the sites for the simulation
# taking this from 02_within_between.R

sumdata <- read_csv("Data/sumdata.csv")
sumdata %>% filter(Type == "Total carbon" & Horizon == "Mineral")

# 2 First  Total carbon Mineral Between  19.0
# 3 First  Total carbon Mineral Within   17.3

sumdata %>% filter(Type == "Change in Total carbon" & Horizon == "Mineral")

# First Change in Total carbon Mineral Mean    0.230

# multiply by 10 to get a change over 10 years


## ---- knobs ---------------------------------------------------------------
R <- 200
Ns <- c(seq(20, 1000, by = 20))   # sites
Reps <- seq(2,10, by = 2)                                               # cores per site

# simulation parameters (tune as desired)
sigma_b    <- 19.0
sigma_w    <- 17.3
delta_mean <- 2.30 # change in carbon after 10 years
sigma_delta <- 18.4 # between SD pooled for total carbon change
rho_core   <- 0.135
paired     <- TRUE

# parallel plan
workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "8"))
plan(multisession, workers = workers)

## ---- helper: one run -----------------------------------------------------
run_one_sim <- function(N_sites, n_rep, iter) {
  # simulate core-level T1/T2 data
  sim <- simulate_two_time(
    n_sites    = N_sites,
    n_rep      = n_rep,
    sigma_b    = sigma_b,
    sigma_w    = sigma_w,
    delta_mean = delta_mean,
    sigma_delta= sigma_delta,
    paired     = paired,
    rho_core   = rho_core
    # seed intentionally not set here; reproducibility handled by furrr seed
  )

  # reshape to the columns your mdc_pop_full_yr() expects
  sim_long <- sim %>%
    pivot_longer(c(T1, T2), names_to = "meas_label", values_to = "TC") %>%
    mutate(
      meas_num  = ifelse(meas_label == "T1", 1L, 2L),
      meas_date = as.Date("2020-01-01") + (meas_num - 1) * 365,
      nfi_plot  = as.integer(site),
      pit_num   = as.integer(core)
    ) %>%
    select(nfi_plot, pit_num, meas_num, meas_date, TC)

  # run estimator; guard if pairing fails
  out <- try(
    mdc_pop_full_yr(
      df        = sim_long,
      unit_col  = "nfi_plot",
      time_col  = "meas_num",
      date_col  = "meas_date",
      value_col = "TC",
      repl_col  = "pit_num",
      alpha     = 0.05,
      power     = 0.80,
      horizon_years = 10   # set if you want 10-year MDC/MDD
    ),
    silent = TRUE
  )

  if (inherits(out, "try-error")) {
    # return a single NA row with identifiers so the grid stays rectangular
    tibble(
      N_sites = N_sites, n_rep = n_rep, iter = iter,
      m_units = NA_integer_, mean_change = NA_real_,
      lmm_slope_est_rate = NA_real_, lmm_slope_se = NA_real_, lmm_slope_p = NA_real_,
      within_sd_pooled = NA_real_, between_sd_pooled = NA_real_,
      sd_rate_obs = NA_real_, meas_var_rate_bar = NA_real_,
      within_sd_rate = NA_real_, between_sd_rate = NA_real_,
      SE_pop_rate = NA_real_, MDC_rate_pop_alpha = NA_real_,
      MDD_rate_pop_a_p = NA_real_, rho_core = NA_real_, rho_plot = NA_real_,
      rho_used_for_meas = NA_real_, n_core_pairs = NA_integer_, n_units_with_pairs = NA_integer_,
      MDC_horizon_alpha = NA_real_, MDD_horizon_a_p = NA_real_, horizon_years = NA_real_
    )
  } else {
    out %>%
      mutate(N_sites = N_sites, n_rep = n_rep, iter = iter) %>%
      # keep a tidy subset of key outputs; add/adjust as needed
      select(
        N_sites, n_rep, iter,
        m_units, mean_change,
        lmm_slope_est_rate, lmm_slope_se, lmm_slope_p,
        within_sd_pooled, between_sd_pooled,
        sd_rate_obs, meas_var_rate_bar, within_sd_rate, between_sd_rate,
        SE_pop_rate, MDC_rate_pop_alpha, MDD_rate_pop_a_p,
        rho_core, rho_plot, rho_used_for_meas, n_core_pairs, n_units_with_pairs,
        starts_with("MDC_horizon"), starts_with("MDD_horizon"), horizon_years
      )
  }
}

## ---- parameter grid + parallel map ---------------------------------------
param_grid <- tidyr::crossing(N_sites = Ns, n_rep = Reps, iter = seq_len(R))

results <- future_pmap_dfr(
  param_grid,
  function(N_sites, n_rep, iter) run_one_sim(N_sites, n_rep, iter),
  .options  = furrr_options(seed = TRUE),  # reproducible RNG per worker
  .progress = TRUE
)

## ---- save / quick looks ---------------------------------------------------
saveRDS(results, "Data/Simulation/Fig4_Decreasing_sites_mdd.RDS")
results <-  readRDS("Data/Simulation/Fig4_Decreasing_sites_mdd.RDS")


# Summaries for both variables
summ_change <- results %>%
  group_by(N_sites, n_rep) %>%
  summarise(
    lo = quantile(lmm_slope_est_rate, 0.05, na.rm = TRUE),
    hi = quantile(lmm_slope_est_rate, 0.95, na.rm = TRUE),
    med = median(lmm_slope_est_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "Estimated Change")

summ_mdd <- results %>%
  group_by(N_sites, n_rep) %>%
  summarise(
    lo = quantile(MDD_rate_pop_a_p, 0.05, na.rm = TRUE),
    hi = quantile(MDD_rate_pop_a_p, 0.95, na.rm = TRUE),
    med = median(MDD_rate_pop_a_p, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(metric = "MDD (α=0.05, power=0.8)")


# where does MDD cross true effect
# beyond this crossing threshold, we have >= 80% power to detect the true effect.
# cannot rely on MDD to interpret a single dataset


# crossing where median MDD <= true effect
crossings <- summ_mdd %>%
  # filter(n_rep != 10) %>%
  group_by(n_rep) %>%
  arrange(N_sites) %>%
  filter(med <= delta_mean) %>%
  slice(1) %>%           # first N where MDD < true effect
  ungroup() %>%
  mutate(facet_label = paste0(n_rep, " Microplots"))


# Combine and plot
summ_both <- bind_rows(summ_change, summ_mdd) %>%
  mutate(facet_label = paste0(n_rep, " Microplots"))

ggplot(summ_both  %>%
         filter(n_rep != 10),
       aes(x = N_sites, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med, linetype = metric), size = 1) +
  geom_hline(yintercept = delta_mean, lty = 2) +
  # geom_hline(yintercept = 4.0, lty = 3) +
  # vertical line at "achieves target power" N_sites
  geom_vline(
    data = crossings%>%
      filter(n_rep != 10),
    aes(xintercept = N_sites),
    linetype = 3
  ) +
  facet_grid(~ facet_label, scales = "free_x") +
  labs(
    x = "Number of plots",
    y = expression("Carbon stock change"~(Mg~C~ha^-1~10~yr^-1)),
    color = "Metric", fill = "Metric", linetype = "Metric"
  ) +
  theme_classic(base_size = 13) +
  coord_cartesian(ylim = c(0, 8)) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal"
  ) +
  scale_color_manual(
    values = c(
       "#1A1A1A",       # black/dark
       "#4F6D7A"        # muted blue-gray
    )
  ) +
  scale_fill_manual(
    values = c(
      "#1A1A1A20",     # transparent black tint for ribbon
      "#4F6D7A40"      # transparent blue-gray tint
    )
  ) +
  scale_linetype_manual(
    values = c(
      "solid",
      "dashed"
    )
  )


ggsave("Plots/04_MDD_simulation.png", width = 10, height = 5)



ggplot(summ_both, aes(x = N_sites, color = metric, fill = metric)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(aes(y = med), size = 1) +
  geom_hline(yintercept = delta_mean, lty = 2) +
  facet_wrap(~ n_rep, scales = "free_x") +
  labs(
    title = "Simulated change of +2.3 Mg C ha-1",
    subtitle = "Between site variation ", sigma_b, "; Within site variation = ", sigma_w ,
    x = "Number of sites",
    y = "Carbon Stock Change",
    color = "Metric", fill = "Metric"
  ) +
  theme_minimal(base_size = 13) + coord_cartesian(ylim = c(0, 10), xlim = c(0,600))




# I want to look at how plots, sites, changes variation and MDD -----------


results_summary <- results %>%
  group_by(N_sites, n_rep) %>%
  summarize(
    mean_within = mean(within_sd_pooled, na.rm = TRUE),
    mean_between = mean(between_sd_pooled, na.rm = TRUE),
    mean_SE = mean(SE_pop_rate, na.rm = TRUE),
    mean_MDD = mean(MDD_rate_pop_a_p, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(results_summary, aes(x = N_sites, y = mean_MDD, color = factor(n_rep))) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "Effect of Increasing Number of Sites on MDD",
    x = "Number of Sites",
    y = "Mean MDD",
    color = "Plots per Site"
  ) +
  theme_minimal(base_size = 14)


p1 <- ggplot(results_summary, aes(x = n_rep, y = mean_within, color = factor(N_sites))) +
  geom_line(linewidth = 1) +
  labs(
    title = "Within-Site SD vs. Plots per Site",
    x = "Plots per Site",
    y = "Within-Site SD",
    color = "N sites"
  ) +
  theme_minimal(base_size = 14)

p2 <- ggplot(results_summary, aes(x = N_sites, y = mean_between, color = factor(n_rep))) +
  geom_line(linewidth = 1) +
  labs(
    title = "Between-Site SD vs. Number of Sites",
    x = "Number of Sites",
    y = "Between-Site SD",
    color = "Plots per Site"
  ) +
  theme_minimal(base_size = 14)

library(patchwork)
p1 / p2


results_summary <- results %>%
  group_by(N_sites, n_rep) %>%
  summarize(
    SE = mean(SE_pop_rate, na.rm = TRUE),
    between = mean(between_sd_rate^2, na.rm = TRUE),
    meas = mean(meas_var_rate_bar, na.rm = TRUE),
    .groups = "drop"
  )

results_long <- results_summary %>%
  pivot_longer(cols = c(between, meas),
               names_to = "component",
               values_to = "variance")

ggplot(results_long, aes(x = N_sites, y = variance, color = component)) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~ n_rep, scales = "free_y") +
  labs(
    title = "Variance Components Informing SE_pop_rate",
    x = "Number of Sites",
    y = "Variance",
    color = "Component"
  ) +
  theme_minimal(base_size = 14)


ggplot(results_summary, aes(x = N_sites)) +
  geom_line(aes(y = SE, color = "SE"), linewidth = 1.3) +
  geom_line(aes(y = between / sqrt(N_sites), color = "Between / sqrt(m)"), linewidth = 1) +
  geom_line(aes(y = meas / sqrt(N_sites), color = "Meas / sqrt(m)"), linewidth = 1) +
  facet_wrap(~ n_rep, scales = "free_y") +
  labs(
    title = "How Variance Components Combine to Form SE_pop_rate",
    x = "Number of Sites",
    y = "Standard Error Contribution",
    color = "Component"
  ) +
  theme_minimal(base_size = 14)


# Increasing n_rep mainly reduces measurement variance

# Between-rate variance: true environmental heterogeneity in rate of change across your sites.
# var(true rate at site i)



results_summary <- results %>%
  group_by(N_sites, n_rep) %>%
  summarize(
    SE      = mean(SE_pop_rate, na.rm = TRUE),
    between = mean(between_sd_rate^2, na.rm = TRUE),
    meas    = mean(meas_var_rate_bar, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(results_summary,
       aes(x = N_sites, y = SE, color = factor(n_rep))) +
  geom_line(linewidth = 1) +
  geom_point() +
  labs(
    title = "SE of Population Mean Rate vs Number of Sites",
    x = "Number of Sites",
    y = "SE_pop_rate",
    color = "Plots per Site (n_rep)"
  ) +
  theme_minimal(base_size = 14)


results_long <- results_summary %>%
  select(N_sites, n_rep, between, meas) %>%
  pivot_longer(cols = c(between, meas),
               names_to = "component",
               values_to = "variance")

ggplot(results_long,
       aes(x = N_sites, y = variance, color = factor(n_rep))) +
  geom_line(linewidth = 1) +
  geom_point() +
  facet_wrap(~ component, scales = "free_y") +
  labs(
    title = "Variance Components Informing SE_pop_rate",
    x = "Number of Sites",
    y = "Variance",
    color = "Plots per Site (n_rep)"
  ) +
  theme_minimal(base_size = 14)


results_long <- results_summary %>%
  select(N_sites, n_rep, between, meas) %>%
  pivot_longer(cols = c(between, meas),
               names_to = "component",
               values_to = "variance")

ggplot(results_long,
       aes(x = N_sites, y = variance,
           color = component, linetype = factor(n_rep))) +
  geom_line(linewidth = 1) +
  labs(
    title = "Between-site vs Measurement Variance Components",
    x = "Number of Sites",
    y = "Variance (rate units)",
    color = "Variance Component",
    linetype = "Plots per Site"
  ) +
  theme_minimal(base_size = 14)


results_contrib <- results_summary %>%
  mutate(
    between_contrib = between / N_sites,
    meas_contrib = meas / N_sites
  ) %>%
  pivot_longer(cols = c(between_contrib, meas_contrib),
               names_to = "component",
               values_to = "SE_contrib")

ggplot(results_contrib,
       aes(x = N_sites, y = SE_contrib,
           color = component, linetype = factor(n_rep))) +
  geom_line(linewidth = 1) +
  labs(
    title = "Contributions to SE_pop_rate by Variance Component",
    x = "Number of Sites",
    y = "Component Contribution to SE^2",
    color = "Component",
    linetype = "Plots per Site"
  ) +
  theme_minimal(base_size = 14)



results_stack <- results_summary %>%
  mutate(total = between + meas) %>%
  pivot_longer(cols = c(between, meas),
               names_to = "component",
               values_to = "variance")

ggplot(results_stack,
       aes(x = N_sites, y = variance,
           fill = component)) +
  geom_area(alpha = 0.7, position = "stack") +
  facet_wrap(~ n_rep, scales = "fixed") +
  labs(
    title = "Stacked Variance Components",
    x = "Number of Sites",
    y = "Variance (rate units)",
    fill = "Component"
  ) +
  theme_minimal(base_size = 14)




ggplot(results_long,
       aes(x = N_sites, y = variance,
           color = component, linetype = factor(n_rep))) +
  geom_line(linewidth = 1) +
  labs(title = "Between vs Measurement Variance (rate units)",
       x = "Number of Sites",
       y = "Variance",
       color = "Component",
       linetype = "n_rep") +
  theme_minimal(base_size = 14)


# MDD (annual) vs N_sites, faceted by cores per site
results %>%
  ggplot(aes(x = as.factor(N_sites), y = MDD_rate_pop_a_p)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ n_rep, scales = "free_x") +
  labs(x = "Number of sites", y = "Annual MDD", title = "Annual MDD vs. Sites and Cores per Site")

# MDC over 10-year horizon vs N_sites (if horizon_years was set)
results %>%
  ggplot(aes(x = as.factor(N_sites), y = MDC_horizon_alpha)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ n_rep, scales = "free_x") +
  labs(x = "Number of sites", y = "10-year MDC", title = "10-year MDC vs. Sites and Cores per Site")

# Summaries by (N_sites, n_rep)
summary_grid <- results %>%
  group_by(N_sites, n_rep) %>%
  summarise(
    mean_MDD = mean(MDD_rate_pop_a_p, na.rm = TRUE),
    mean_MDC = mean(MDC_rate_pop_alpha, na.rm = TRUE),
    mean_SE  = mean(SE_pop_rate, na.rm = TRUE),
    .groups = "drop"
  )
print(summary_grid)







