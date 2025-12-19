# An exploration of variation by ecozone to test whether it reduces between site variation.
# A: Robert Buchkowski & Alex Polussa

# Load packages:
library(lubridate)
library(lme4)
library(tidyverse)

# Load in the necessary data:
source("Scripts/00_load_data.R")
source("Scripts/99_functions.R")

# Create a data frame of change:

gpremeas_calc_all = gpremeas %>%
  # Some sites have been sampled 3 times (i.e., 2003, 2013, 2023). Let's get rid of these data for now so we only have ~10 year time steps.
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
    Change_TC = (last_TC - first_TC) / years_diff,
    Change_perC = (last_perC - first_perC) / years_diff,
    Change_BD = (last_BD - first_BD) / years_diff,
    Change_CF = (last_CF - first_CF) / years_diff,
    .groups = "drop"
  ) %>%
  filter(years_diff > 0)

# Add in the ecozone data:
gpremeas_calc_all = gpremeas_calc_all %>%
  # add ecozone/site covars
  left_join(
    sitedata %>%
      dplyr::select(nfi_plot, ecozone, slope, aspect, land_cover, land_pos, veg_type) %>%
      dplyr::distinct(nfi_plot, .keep_all = TRUE),
    by = "nfi_plot"
  )

# Fair compare between and within site vairation:

# Get the ecozone vector:
ecozones <- sort(unique(gpremeas_calc_all$ecozone))

results <- vector(mode = "list", length(ecozones))

for(i in 1:length(results)){

  Nplots = gpremeas_calc_all %>%
    filter(ecozone == ecozones[i]) %>%
    pull(nfi_plot) %>% unique() %>% length()

  if(Nplots > 10){
    results[[i]] = gpremeas_calc_all %>%
      filter(Horizon == "Mineral") %>%
      filter(ecozone == ecozones[i]) %>% bwvr_alex(data = .,
                                                   vars = c(
                                                     "first_TC", "last_TC", "first_perC", "last_perC",
                                                     "first_BD", "last_BD", "first_CF", "last_CF",
                                                     "Change_TC","Change_perC","Change_BD","Change_CF"
                                                   ),
                                                   unit_col = "nfi_plot") %>%
      select(Variable, Overall_Mean_Weighted, contains("_Pooled")) %>%
      mutate(Horizon = "Mineral") %>%
      bind_rows(
        gpremeas_calc_all %>%
          filter(Horizon == "Organic") %>%
          filter(ecozone == ecozones[i]) %>% bwvr_alex(data = .,
                                                       vars = c(
                                                         "first_TC", "last_TC", "first_perC", "last_perC",
                                                         "first_BD", "last_BD", "first_CF", "last_CF",
                                                         "Change_TC","Change_perC","Change_BD","Change_CF"
                                                       ),
                                                       unit_col = "nfi_plot") %>%
          select(Variable, Overall_Mean_Weighted, contains("_Pooled")) %>%
          mutate(Horizon = "Organic")
      ) %>%
      separate(Variable, into = c("Time", "name"), sep = "_") %>%
      mutate(name = case_when(
        name == "TC"    ~ "Total carbon",
        name == "BD"    ~ "Bulk density",
        name == "perC"  ~ "Percent carbon",
        name == "CF"  ~ "Coarse fraction",
        TRUE            ~ name  # keep original if no match
      )) %>%
      mutate(Time = case_when(
        Time == "first" ~ "First",
        Time == "last" ~ "Second",
        TRUE ~ Time
      )) %>%
      rename(Type = name) %>%
      pivot_longer(contains("_")) %>%
      mutate(name = case_when(
        grepl("Mean", name) ~ "Mean",
        grepl("Between", name) ~ "Between",
        grepl("Within", name) ~ "Within"
      )) %>%
      mutate(Type = ifelse(
        Time == "Change", paste0("Change in ", Type), Type
      )) %>%
      mutate(Time = ifelse(
        Time == "Change", "First", Time
      )) %>%
      mutate(Type =
               factor(Type, levels =
                        c(
                          "Total carbon", "Change in Total carbon",
                          "Percent carbon", "Change in Percent carbon",
                          "Bulk density", "Change in Bulk density",
                          "Coarse fraction", "Change in Coarse fraction"
                        )
               )) %>%
      mutate(ecozone = ecozones[i],
             N = Nplots)
  }
}

results = do.call("rbind", results)

results

# Get the data for 100 random draws of comparable sample size not associated with ecozone to compare:

samsize = unique(results$N)

nfiplots = unique(gpremeas_calc_all$nfi_plot)

random_selection_overall <- vector("list", 100)

for(rso in 1:100){
  random_selection <- vector(mode = "list", length(samsize))

  for(i in 1:length(random_selection)){

    nfiplotscur = sample(nfiplots, samsize[i])

    random_selection[[i]] = gpremeas_calc_all %>%
      filter(Horizon == "Mineral") %>%
      filter(nfi_plot %in% nfiplotscur) %>% bwvr_alex(data = .,
                                                      vars = c(
                                                        "first_TC", "last_TC", "first_perC", "last_perC",
                                                        "first_BD", "last_BD", "first_CF", "last_CF",
                                                        "Change_TC","Change_perC","Change_BD","Change_CF"
                                                      ),
                                                      unit_col = "nfi_plot") %>%
      select(Variable, Overall_Mean_Weighted, contains("_Pooled")) %>%
      mutate(Horizon = "Mineral") %>%
      bind_rows(
        gpremeas_calc_all %>%
          filter(Horizon == "Organic") %>%
          filter(nfi_plot %in% nfiplotscur) %>% bwvr_alex(data = .,
                                                          vars = c(
                                                            "first_TC", "last_TC", "first_perC", "last_perC",
                                                            "first_BD", "last_BD", "first_CF", "last_CF",
                                                            "Change_TC","Change_perC","Change_BD","Change_CF"
                                                          ),
                                                          unit_col = "nfi_plot") %>%
          select(Variable, Overall_Mean_Weighted, contains("_Pooled")) %>%
          mutate(Horizon = "Organic")
      ) %>%
      separate(Variable, into = c("Time", "name"), sep = "_") %>%
      mutate(name = case_when(
        name == "TC"    ~ "Total carbon",
        name == "BD"    ~ "Bulk density",
        name == "perC"  ~ "Percent carbon",
        name == "CF"  ~ "Coarse fraction",
        TRUE            ~ name  # keep original if no match
      )) %>%
      mutate(Time = case_when(
        Time == "first" ~ "First",
        Time == "last" ~ "Second",
        TRUE ~ Time
      )) %>%
      rename(Type = name) %>%
      pivot_longer(contains("_")) %>%
      mutate(name = case_when(
        grepl("Mean", name) ~ "Mean",
        grepl("Between", name) ~ "Between",
        grepl("Within", name) ~ "Within"
      )) %>%
      mutate(Type = ifelse(
        Time == "Change", paste0("Change in ", Type), Type
      )) %>%
      mutate(Time = ifelse(
        Time == "Change", "First", Time
      )) %>%
      mutate(Type =
               factor(Type, levels =
                        c(
                          "Total carbon", "Change in Total carbon",
                          "Percent carbon", "Change in Percent carbon",
                          "Bulk density", "Change in Bulk density",
                          "Coarse fraction", "Change in Coarse fraction"
                        )
               )) %>%
      mutate(N = samsize[i])
  }

  random_selection_overall[[rso]] = do.call("rbind", random_selection)
  cat("Done", rso, "\n")
}

random_selection_overall_grouped = do.call("rbind", random_selection_overall)

forplot = random_selection_overall_grouped %>% group_by(Time, Type, Horizon, name, N) %>%
  summarize(mean = mean(value),
            sd = sd(value)) %>%
  ungroup() %>%
  full_join(
    results,by = join_by(Time, Type, Horizon, name, N)
  ) %>%

  # Start by only looking at total carbon stocks
  filter(Type == "Total carbon") %>%
  # Start by only looking at the first measurement
  filter(Time == "First") %>%
  filter(name != "Mean") %>%
  mutate(ymin = mean-sd,
         ymax = mean+sd) %>%
  mutate(FL = paste0(ecozone,"\n(N=", N, ")"))


# Order by sample size:
forplot = forplot %>%
  mutate(FL = factor(FL, levels = forplot %>%
                       select(N, FL) %>%
                       distinct() %>%
                       arrange(N) %>%
                       pull(FL)))

png("Plots/Figure5.png", width = 9, height = 4, units = "in", res = 600)
forplot %>%
  ggplot(aes(x = FL, y = value, fill = Horizon)) + geom_bar(stat = "identity") + theme_classic() + geom_pointrange(aes(y = mean, ymin = ymin, ymax = ymax)) +
  ylab("Standard deviation") + xlab("Ecozone") + facet_wrap(.~interaction(name, Horizon, sep = " "), scales = "free_y", nrow = 4, ncol = 2, dir = "h") + scale_fill_manual(name = "Horizon", values = c("grey", "brown"), guide = F)
dev.off()
