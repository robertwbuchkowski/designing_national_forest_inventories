# Load Data

library(lubridate)
library(lme4)
library(tidyverse)

{
  Mgpdata = read_csv("Data/all_gp_soil_mineral_sample.csv") %>%
    mutate(meas_date = ymd(meas_date)) %>% # Make the measurement date a date
    filter(!(nfi_plot == 1135306 & meas_num == 1 & pit_num == "MP4" & sample_num == 2))

  # Load in the organic horizon data:

  Ogpdata = read_csv("Data/all_gp_for_flr_org_sample.csv") %>%
    mutate(meas_date = ymd(meas_date)) # Make the measurement date a date

  # Load in the site information:

  sitedata = read_csv("Data/all_gp_site_info.csv") %>%
    mutate(meas_date = ymd(meas_date)) # Make the measurement date a date

  # This function calculates the sample proportion in 0 to 15-cm
  proportion_in_target_range_vec <-
    function(sample_lower, sample_topper,
             target_lower = 0, target_topper = 15) {
      # Calculate overlap bounds
      overlap_lower <- pmax(sample_lower, target_lower)
      overlap_upper <- pmin(sample_topper, target_topper)

      # Calculate overlap length (0 if no overlap)
      overlap_length <- ifelse(overlap_upper > overlap_lower,
                               overlap_upper - overlap_lower, 0)

      # Calculate sample length
      sample_length <- sample_topper - sample_lower

      # Avoid division by zero
      proportion <- ifelse(sample_length > 0, overlap_length / sample_length, NA)

      return(proportion)
    }


  Mgpdata = Mgpdata %>%
    group_by(nfi_plot) %>% # by plot
    filter(!is.na(layer_cc)) %>% # get rid of non-data lines
    filter(pit_num %in% c("MP1", "MP2", "MP3", "MP4")) %>% # Get rid of soil pit records, keep only microplots

    # First, -7 or -1 means missing data,
    # while -9 means there were no disc rocks.
    # Just set these to zero mass as an approximation.
    mutate(mass_cobble = ifelse(mass_cobble < 0, 0, mass_cobble),
           mass_disc_rocks = ifelse(mass_disc_rocks < 0, 0, mass_disc_rocks)) %>%

    # Calculate the coarse fraction:
    mutate(coarse_fraction = (mass_gravel + mass_cobble + mass_disc_rocks)/2.65) %>%
    # Replace missing gravel values with 0
    mutate(coarse_fraction = ifelse(coarse_fraction < 0, 0, coarse_fraction)) %>%

    # Calculate the bulk density using NFI formula:
    mutate(bulk_density_2mm_custom = mass_2mm*(1-soil_moisture)/(volume- coarse_fraction)) %>%

    # Convert coarse fraction to CF bulk density:
    mutate(coarse_fraction = coarse_fraction*2.65/volume) %>%

    # Use the NFI bulk density if the raw data are missing:
    mutate(bulk_density_2mm_chosen = ifelse(bulk_density_2mm_custom < 0 | bulk_density_2mm_custom > 5, bulk_density_2mm, bulk_density_2mm_custom)) %>%

    # Get the layer % carbon:
    mutate(perC = ifelse(toc > 0, toc, tc)*0.1) %>% # choose the available data and convert from g/kg to %

    # Calculate the layer_cc using chosen bulk density:
    mutate(layer_ccc_gravadj = perC/10* # convert back to g/kg
             bulk_density_2mm_chosen*
             (sample_bottom-sample_upper)) %>%

    # CALCULATE the Adjusted carbon:

    mutate(layer_cc_gravadj =
             layer_cc *
             (volume/(volume -((mass_cobble + mass_disc_rocks)/2.65)))) %>%

    filter(layer_cc_gravadj >0 | layer_ccc_gravadj > 0) %>%

    # Get rid of 2 samples with probable error
    # when layer_cc_gravadj is negative because
    # the estimated volume of gravel and disc rocks
    # is larger than the sample volume:
    filter(layer_ccc_gravadj >= 0) %>%

    # Calculate the proportion of the sample layer in 0 to 15-cm:
    mutate(prop0to15 = proportion_in_target_range_vec(sample_upper,sample_bottom)) %>%
    # Remove sites missing the properties that we want:
    filter(layer_ccc_gravadj >=0 & perC >=0 & bulk_density_2mm_chosen >=0 & coarse_fraction >=0) %>%
    # Discard all the rows where the samples where the proportion is zero:
    filter(prop0to15 > 0) %>%
    # Sum up across all samples within the pit:
    group_by(nfi_plot, meas_date, pit_num) %>%
    # Calculate the thickness and CC in 0-15
    reframe(Thickness_0_15 = sum((sample_bottom - sample_upper)*prop0to15),
            CC0_15 = sum(layer_ccc_gravadj*prop0to15),
            CC0_15_old = sum(layer_cc_gravadj*prop0to15),
            top = min(sample_upper),
            perC_m = sum(perC*prop0to15/sum(prop0to15)),
            bulk_density_2mm_m = sum(bulk_density_2mm_chosen*prop0to15/sum(prop0to15)),
            coarse_fraction = sum(coarse_fraction*prop0to15/sum(prop0to15))) %>%
    filter(top == 0) %>% # Keep only plots with surface samples.
    # Calculate the adjusted carbon content for all layers,
    # but do not average across soil pits as done by NFI.
    # Also, do not correct for soil pit cobble and stone, because that is a 'site' level measurement
    mutate(CC0_15_ADJ = CC0_15/Thickness_0_15 * 15) %>%
    # Convert from kg/m2 to Mg/ha:
    mutate(CC0_15_ADJ = CC0_15_ADJ * 10) %>%
    group_by(nfi_plot) %>%
    mutate(meas_year = year(meas_date)) %>%
    arrange(meas_year, .by_group = TRUE) %>%
    mutate(
      meas_num = dense_rank(meas_year)  # Assign unique visit number per year
    ) %>%
    ungroup() %>%
    rename(layer_cc_m = perC_m) %>% # Keep the name the same
    # Clean up the data frame:
    select(nfi_plot, meas_date, meas_num, pit_num, CC0_15_ADJ,layer_cc_m,bulk_density_2mm_m, coarse_fraction)


  Ogpdata = Ogpdata %>%
    group_by(nfi_plot) %>% # by plot
    filter(!is.na(layer_cc_total)) %>% # get rid of non-data lines
    filter(pit_num %in% c("MP1", "MP2", "MP3", "MP4")) %>% # Get rid of soil pit records, keep only microplots
    filter(layer_cc_total >= 0) %>% # get rid of negative flags that NFI uses for missing data

    # Calculate cf:
    mutate(CF = ifelse(mass_gt8mm > 0, ifelse(mass_char_gt8mm > 0, mass_gt8mm + mass_char_gt8mm, mass_gt8mm), 0)) %>%

    # Calculate ff:
    mutate(FF = ifelse(mass_8mm > 0, mass_8mm, mass_total)) %>%

    # Calculate %C:
    mutate(TC = ifelse(tc_8mm > 0, tc_8mm, toc_8mm)) %>%

    # Sum up across all samples within the pit:
    group_by(nfi_plot, meas_date, pit_num) %>%
    # Calculate the thickness and CC in 0-15
    summarize(cc_total = sum(layer_cc_total),
              bd_average = sum(FF)/sum(volume),
              tc_average = mean(TC* 0.1),# Convert from g/kg to percent
              cf_average = sum(CF)/sum(volume),
              .groups = "drop") %>%
    # Convert from kg/m2 to Mg/ha:
    mutate(cc_total = cc_total * 10) %>%
    group_by(nfi_plot) %>%
    mutate(meas_year = year(meas_date)) %>%
    arrange(meas_year, .by_group = TRUE) %>%
    mutate(
      meas_num = dense_rank(meas_year)  # Assign unique visit number per year
    ) %>%
    ungroup() %>%
    # Clean up the data frame:
    select(nfi_plot, meas_date, meas_num, pit_num, cc_total,tc_average,bd_average,cf_average)

  # Remove duplicated rows:
  Mgpdata = Mgpdata %>% distinct()

  # Remove mineral samples with more than 17% carbon. These are NOT necessarily associated with organic soils and appear to be field sampling errors where the organic horizon was partly included in the mineral horizon.
  Mgpdata = Mgpdata %>% filter(layer_cc_m < 17)


  #### * Latest ####
  Mgpdata_latest = Mgpdata %>%
    group_by(nfi_plot) %>%
    filter(meas_num == 2)

  Ogpdata_latest = Ogpdata %>%
    group_by(nfi_plot) %>%
    filter(meas_num == 2)


  #### * across time #####

  Mgpremeas = Mgpdata %>%
    # Create a list of nfi_plot where at least 2 measurements were taken
    filter(meas_num > 1) %>%
    select(nfi_plot) %>%
    distinct() %>%
    left_join( # Join to keep only those nfi plots
      Mgpdata,by = join_by(nfi_plot)
    ) %>%
    rename(TC = CC0_15_ADJ) %>%
    rename(perC = layer_cc_m)%>%
    rename(BD = bulk_density_2mm_m) %>%
    rename(CF = coarse_fraction) %>%
    mutate(Horizon = "Mineral")

  Ogpremeas = Ogpdata %>%
    # Create a list of nfi_plot where at least 2 measurements were taken
    filter(meas_num > 1) %>%
    select(nfi_plot) %>%
    distinct() %>%
    left_join( # Join to keep only those nfi plots
      Ogpdata,by = join_by(nfi_plot)
    ) %>%
    rename(TC = cc_total) %>%
    rename(perC = tc_average)%>%
    rename(BD = bd_average) %>%
    rename(CF = cf_average) %>%
    mutate(Horizon = "Organic")

  # Combine the mineral and organic horizons:
  gpremeas = Mgpremeas %>%
    bind_rows(Ogpremeas)

  # Check for duplicate removal:
  if(gpremeas %>%
     group_by(nfi_plot, meas_num, pit_num, Horizon) %>%
     summarize(N = n(), .groups = "drop") %>% pull(N) %>% max() != 1){
    print("Duplicated plots found. Check data read.")
  }

}


print("NOTE: The meas_num variable is Horizon unique. So meas_num = 1 for mineral horizons may not be the same date as meas_num = 1 for organic horizons. This happens if only one horizon was successfully measured during a visit.")

print("Mgpdata is full mineral dataset")
print("Ogpdata is full Organic dataset")

print("Mgpremeas is Mineral remeasurement dataset")
print("Ogpremeas is Organic remeasurement  dataset")

print("Mgpdata_latest is the latest mineral measurement")

print("gpremeas is the combined organic and mineral remeasurement dataset")


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

print("gpremeas_calc_all is the combined organic and mineral remeasurement dataset subset to include only the plots with remeasured microplots.")
