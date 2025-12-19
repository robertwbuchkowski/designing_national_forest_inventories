# A script to analyze the 4 and less than 4 microplots:

# Load packages:
library(lubridate)
library(lme4)
library(nlme)
library(tidyverse)
library(DHARMa)
library(glmmTMB)

# Load in the necessary data:
source("Scripts/00_load_data.R")

alldata = gpremeas %>%
  # Keep the first two measurement times for each plot:
  filter(meas_num < 3) %>%
  group_by(nfi_plot, meas_num, Horizon) %>%
  summarise(TC = mean(TC),
            perC = mean(perC),
            BD = mean(BD),
            CF = mean(CF),
            N = n(),
            .groups = "drop") %>%
  pivot_longer(TC|perC|BD|CF) %>%
  mutate(name = case_when(
    name == "TC"    ~ "Carbon stock",
    name == "BD"    ~ "Bulk density",
    name == "perC"  ~ "Percent carbon",
    name == "CF"  ~ "Coarse fraction",
    TRUE            ~ name  # keep original if no match
  ))

# Mineral horizon carbon ------------

# Parameter:
param = "Carbon stock"
horizon_select = "Mineral"

# Mineral horizon data:
curdata = alldata %>%
  filter(Horizon == horizon_select) %>%
  filter(name == param) %>%
  mutate(meas_num = paste(meas_num))

curdata %>%
  group_by(N) %>%
  summarize(vmean = mean(value),
            vmed = median(value))

# Bias:
curdata %>%
  group_by(N) %>%
  summarize(vmean = mean(value)) %>%
  mutate(gmean = mean(curdata$value)) %>%
  mutate(bias = vmean - gmean)

m1 <- glmmTMB(value ~ N + meas_num + (1|nfi_plot), data = curdata,
              dispformula = ~N + meas_num,
              family = Gamma(link = "log"))

cat("--------------------------\n",param, "for", horizon_select,"\n--------------------------")
summary(m1)

lag(predict(m1, newdata = data.frame(N = c(1,2,3,4), meas_num = 1), type = "response", re.form = NA)) - predict(m1, newdata = data.frame(N = c(1,2,3,4), meas_num = 1), type = "response", re.form = NA)


simres <- simulateResiduals(m1, n = 2000)

DHARMa::testSimulatedResiduals(simres)
DHARMa::testCategorical(simres, catPred = curdata$N)
DHARMa::testCategorical(simres, catPred = curdata$meas_num)

# Run the same model with randomly chosen microplots instead of the average to confirm the same coefficients:

samples <- vector("list", 100)

for(i in 1:length(samples)){
  alldata_subset = gpremeas %>%
    # Keep the first two measurement times for each plot:
    filter(meas_num < 3) %>%
    group_by(nfi_plot, meas_num, Horizon) %>%
    summarise(N = n(),
              .groups = "drop") %>%
    left_join(
      # Take 1 random sample for the plot value:
      gpremeas %>%
        # Keep the first two measurement times for each plot:
        filter(meas_num < 3) %>%
        group_by(nfi_plot, meas_num, Horizon) %>%
        slice_sample(n = 1) %>%
        select(nfi_plot, meas_num, Horizon, TC, perC, BD, CF)
    ) %>%
    pivot_longer(TC|perC|BD|CF) %>%
    mutate(name = case_when(
      name == "TC"    ~ "Carbon stock",
      name == "BD"    ~ "Bulk density",
      name == "perC"  ~ "Percent carbon",
      name == "CF"  ~ "Coarse fraction",
      TRUE            ~ name  # keep original if no match
    )) %>%
    filter(Horizon == horizon_select) %>%
    filter(name == param) %>%
    mutate(meas_num = paste(meas_num))

  m1 <- glmmTMB(value ~ N + meas_num + (1|nfi_plot), data = alldata_subset,
                dispformula = ~N + meas_num,
                family = Gamma(link = "log"))

  samples[[i]] = data.frame(summary(m1)$coefficients$cond) %>% rownames_to_column() %>%
    mutate(Rep = i)
}

do.call("rbind",samples) %>%
  tibble() %>%
  filter(rowname == "N") %>%
  ggplot(aes(x = Estimate)) + geom_histogram()

do.call("rbind",samples) %>%
  tibble() %>%
  filter(rowname == "N") %>%
  ggplot(aes(x = `Pr...z..`)) + geom_histogram()

do.call("rbind",samples) %>%
  tibble() %>%
  filter(rowname == "N") %>%
  filter(`Pr...z..` < 0.05) %>% dim()


# Organic horizon carbon ------------

# Parameter:
param = "Carbon stock"
horizon_select = "Organic"

# Mineral horizon data:
curdata = alldata %>%
  filter(Horizon == horizon_select) %>%
  filter(name == param) %>%
  mutate(meas_num = paste(meas_num))


m1 <- glmmTMB(value ~ N + meas_num + (1|nfi_plot), data = curdata,
              dispformula = ~N + meas_num,
              family = lognormal(link = "log"))

cat("--------------------------\n",param, "for", horizon_select,"\n--------------------------")
summary(m1)


simres <- simulateResiduals(m1, n = 2000)

DHARMa::testSimulatedResiduals(simres)
DHARMa::testCategorical(simres, catPred = curdata$N)
DHARMa::testCategorical(simres, catPred = curdata$meas_num)

# Not significant in most of the model versions.

# Plot the result -------

pp1 = alldata %>%
  filter(name == param) %>%
  ggplot(aes(fill = Horizon)) + facet_wrap(.~Horizon, scales = "free_y", nrow = 4, ncol = 2, dir = "v") + geom_boxplot(aes(x = paste0(N), y = value))  +
  stat_summary(aes(x = paste0(N), y = value), fun = mean, geom = "point",
               color = "black",fill = "cyan", size = 3, shape = 23) + theme_classic() + theme(legend.position = "none") + ylab(param) + xlab("Microplots") + scale_fill_manual(name = "Horizon", values = c("grey", "brown")) + geom_text(data = tibble(x = 1.5, y = c(150, 200), Horizon = c("Mineral", "Organic"), txt = c("Microplots*\nMeasurement*", "ns")), aes(label = txt, x=x,y=y)) + ylab(expression(Carbon~stock~(Mg~ha^-1)))

png("Plots/Figure7.png", width = 6, height = 4, units = "in", res = 600)
pp1
dev.off()
