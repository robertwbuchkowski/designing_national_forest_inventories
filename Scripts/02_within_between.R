# A script to produce figure 2 showing within and between site variation
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
    Change_TC = (last_TC - first_TC) / years_diff*10,
    Change_perC = (last_perC - first_perC) / years_diff*10,
    Change_BD = (last_BD - first_BD) / years_diff*10,
    Change_CF = (last_CF - first_CF) / years_diff*10,
    .groups = "drop"
  ) %>%
  filter(years_diff > 0)


# Calculate the total numbers included:
gpremeas_calc_all %>% group_by(nfi_plot, Horizon) %>% count() %>% group_by(Horizon) %>% count()

gpremeas_calc_all %>% select(nfi_plot) %>% distinct() %>% nrow()
gpremeas %>% select(nfi_plot) %>% distinct() %>% nrow()

# Regression to the mean:

png("Plots/FigureS1.png", width = 8, height = 8, units = "in", res = 600)
cowplot::plot_grid(
  gpremeas_calc_all %>%
    group_by(nfi_plot, Horizon) %>%
    summarize(first_TC = mean(first_TC),
              Change_TC = mean(Change_TC)) %>%
    ggplot(aes(x = first_TC, y = Change_TC)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(color = Horizon)) +
    facet_wrap(.~Horizon) + theme_classic() + scale_color_manual(name = "Horizon", values = c("grey", "brown")) +
    stat_smooth(method = "lm", color = "black") + ylab(expression(Change~"in"~carbon~stock~(Mg~ha^-1~10~yr^-1))) + xlab(expression(t[0]~carbon~stock~(Mg~ha^-1))),

  gpremeas_calc_all %>%
    mutate(calc_last_TC = pmax(0,first_TC + Change_TC)) %>%
    group_by(nfi_plot, Horizon) %>%
    summarize(first_TC = mean(first_TC),
              last_TC = mean(calc_last_TC),
              Change_TC = mean(Change_TC)) %>%
    mutate(average_TC = (first_TC + last_TC)/2) %>%
    ggplot(aes(x = average_TC, y = Change_TC)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(aes(color = Horizon)) +
    facet_wrap(.~Horizon) + theme_classic() + scale_color_manual(name = "Horizon", values = c("grey", "brown")) +
    stat_smooth(method = "lm", color = "black") + ylab(expression(Change~"in"~carbon~stock~(Mg~ha^-1~10~yr^-1))) + xlab(expression(Average~t[0]~"&"~t[1]~carbon~stock~(Mg~ha^-1))),
  nrow = 2, labels = "AUTO"
)
dev.off()

# Produce a data frame of the within and between variations for plotting:

sumdata = bwvr_alex(data = gpremeas_calc_all %>%
                      filter(Horizon == "Mineral"),
                    vars = c(
                      "first_TC", "last_TC", "first_perC", "last_perC",
                      "first_BD", "last_BD", "first_CF", "last_CF",
                      "Change_TC","Change_perC","Change_BD","Change_CF"
                    ),
                    unit_col = "nfi_plot") %>%
  select(Variable, Overall_Mean_Weighted, contains("_Pooled")) %>%
  mutate(Horizon = "Mineral") %>%
  bind_rows(
    bwvr_alex(data = gpremeas_calc_all %>%
                filter(Horizon == "Organic"),
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
    name == "TC"    ~ "Carbon stock",
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
                      "Carbon stock", "Change in Carbon stock",
                      "Percent carbon", "Change in Percent carbon",
                      "Bulk density", "Change in Bulk density",
                      "Coarse fraction", "Change in Coarse fraction"
                    )
           ))

sumdata %>%
  pivot_wider() %>%
  mutate(CVbet = 100*Between/Mean,
         CVwit = 100*Within/Mean) %>%
  View()

write.csv(sumdata, "Data/sumdata.csv",  row.names = F)


# Get rate changes to over 10 year time period

sumdata %>% filter(str_detect(Type, "Change") & name == "Mean")

# Change Mineral: 2.3 over 10 years
# Change Organic: -0.867 over 10 years

png("Plots/Figure2.png", width = 12, height = 6, units = "in", res = 600)
cowplot::plot_grid(
  sumdata %>%
    filter(name != "Mean") %>%
    filter(!grepl("Change", Type)) %>%
    ggplot(aes(x = interaction(Time, name, sep = "\n"), y = value, fill = Horizon)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() +
    ylab("Standard deviation (bar) or Mean (line)") + xlab("") + facet_wrap(.~interaction(Horizon, Type, sep = "\n"), scales = "free_y", nrow = 4, ncol = 2, dir = "h") + scale_fill_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") + scale_color_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") + scale_linetype_manual(name = "Measurement", values = c(2,3), guide = "none") +
    geom_hline(
      aes(yintercept = value, linetype = Time),data = sumdata %>%
        filter(name == "Mean") %>%
        filter(!grepl("Change", Type)) %>%
        mutate(LocX = ifelse(Horizon == "Mineral", 1.5, 3.5)) %>%
        mutate(LocX = ifelse(!grepl("Change", Type),
                             ifelse(Time == "First", LocX - 0.02, LocX + 0.02), LocX))),

  sumdata %>%
    filter(name != "Mean") %>%
    filter(grepl("Change", Type)) %>%
    ggplot(aes(x = interaction(name, sep = "\n"), y = value, fill = Horizon)) + geom_bar(stat = "identity", position = "dodge") + theme_classic() +
    ylab("") + xlab("") + facet_wrap(.~interaction(Horizon, Type, sep = "\n"), scales = "free_y", nrow = 4, ncol = 2, dir = "h") + scale_fill_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") + scale_color_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") + scale_linetype_manual(name = "Measurement", values = c(2,3)) +
    geom_hline(
      aes(yintercept = value), linetype = 2, data = sumdata %>%
        filter(name == "Mean") %>%
        filter(grepl("Change", Type))), rel_widths = c(1,1)
)
dev.off()
