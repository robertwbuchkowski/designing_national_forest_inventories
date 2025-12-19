# A script to produce figure 7 showing the relationship between % change and total change
# A: Robert Buchkowski & Alex Polussa

# Load packages:
library(lubridate)
library(lme4)
library(tidyverse)

# Load in the necessary data:
source("Scripts/00_load_data.R")
source("Scripts/99_functions.R")

R2 = tibble(Horizon = c("Mineral", "Organic"),
            R2 = c(
              MuMIn::r.squaredGLMM(nlme::lme(TC~perC,random = ~1|nfi_plot, data = gpremeas %>%
                                               filter(meas_num < 3) %>% filter(Horizon == "Mineral")))[1,"R2m"],
              MuMIn::r.squaredGLMM(nlme::lme(TC~perC,random = ~1|nfi_plot, data = gpremeas %>%
                                               filter(meas_num < 3) %>% filter(Horizon == "Organic")))[1,"R2m"]
            )) %>%
  mutate(label = paste0("R2=",round(R2, 2)),
         x = c(2.5, 10),
         y = 300)

p2 = gpremeas %>%
  filter(meas_num < 3) %>%
  mutate(meas_num = paste(meas_num-1)) %>%
  ggplot() +
  geom_point(aes(x = perC, y = TC, color = Horizon, shape = meas_num)) +
  facet_wrap(.~Horizon, scales = "free", nrow = 2) +
  scale_color_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") +
  scale_shape_discrete(name = "Measurement")+ theme_classic() + xlab("Carbon content (%)") + ylab("Carbon stock (Mg/ha)") + theme(legend.position = "top") +
  geom_text(aes(x = x, y = y, label = label), data = R2)

png("Plots/FigureS2_PERvsTOT.png", width = 4, height = 6, units = "in", res = 600)
p2
dev.off()

# Change result:

change_result = gpremeas %>%
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

# Percent carbon ------

R2 = tibble(Horizon = c("Mineral", "Organic"),
            R2 = c(
              MuMIn::r.squaredGLMM(nlme::lme(Change_TC~Change_perC,random = ~1|nfi_plot, data = change_result %>% filter(Horizon == "Mineral")))[1,"R2m"],
              MuMIn::r.squaredGLMM(nlme::lme(Change_TC~Change_perC,random = ~1|nfi_plot, data = change_result%>% filter(Horizon == "Organic")))[1,"R2m"]
            ))%>%
  mutate(label = paste0("R2=",round(R2, 2)),
         x = c(-20, -70),
         y = c(300,250))

png("Plots/Figure9.png", width = 4, height = 6, units = "in", res = 600)
change_result %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(x = Change_perC, y = Change_TC, color = Horizon)) +
  facet_wrap(.~Horizon, scales = "free", nrow = 2) +
  scale_color_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") +
  scale_shape_discrete(name = "Measurement")+ theme_classic() + xlab(expression(Change ~ "in" ~ carbon ~ content ~ "(" * "%" ~ 10 ~ yr^-1 * ")")) + ylab(expression(Change ~ "in" ~ carbon ~ stock ~ "(" * "Mg"~ ha^-1 ~ 10 ~yr^-1 * ")")) + theme(legend.position = "top")  +
  geom_text(aes(x = x, y = y, label = label), data = R2)
dev.off()

# Bulk density ------

# Change result:

R2 = tibble(Horizon = c("Mineral", "Organic"),
            R2 = c(
              MuMIn::r.squaredGLMM(nlme::lme(Change_TC~Change_BD,random = ~1|nfi_plot, data = change_result %>% filter(Horizon == "Mineral")))[1,"R2m"],
              MuMIn::r.squaredGLMM(nlme::lme(Change_TC~Change_BD,random = ~1|nfi_plot, data = change_result%>% filter(Horizon == "Organic")))[1,"R2m"]
            ))%>%
  mutate(label = paste0("R2=",round(R2, 2)),
         x = c(-2, -1),
         y = c(300,250))

png("Plots/FigureS3_change_BD.png", width = 4, height = 6, units = "in", res = 600)
change_result %>%
  ggplot() +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(aes(x = Change_BD, y = Change_TC, color = Horizon)) +
  facet_wrap(.~Horizon, scales = "free", nrow = 2) +
  scale_color_manual(name = "Horizon", values = c("grey", "brown"), guide = "none") +
  scale_shape_discrete(name = "Measurement")+ theme_classic() + xlab(expression(Change ~ "in" ~ bulk ~ density ~ "(" * "g" ~ cm^-3 ~ 10 ~ yr^-1 * ")")) + ylab(expression(Change ~ "in" ~ carbon ~ stock ~ "(" * "Mg"~ ha^-1 ~ 10 ~ yr^-1 * ")")) + theme(legend.position = "top")  +
  geom_text(aes(x = x, y = y, label = label), data = R2)
dev.off()
