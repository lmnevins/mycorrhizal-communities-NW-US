# -----------------------------------------------------------------------------#
# Process precipitation data from H.J. Andrews Forest meteorological stations
# Original Author: L. McKinley Nevins 
# July 18, 2025
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     ggplot2 v 3.5.1
#                     lubridate v 1.9.4
#                     
# -----------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
require(tidyverse); packageVersion("tidyverse")
require(dplyr); packageVersion("dplyr")
require(ggplot2); packageVersion("ggplot2")
require(lubridate); packageVersion("lubridate")

#################################################################################
#                               Main workflow                                   #
#  Use daily precipitation data from the meteorological stations in HJA to      #
#  calculate as close to 1991-2020 30 year normals of precipitation data that   #
#  I can.                                                                       #
#                                                                               #
#################################################################################

# Using data from here: https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=MS004

# Two meteorological stations: 
# 1. PRIMET, which is the primary station down by the offices, and the station
# that is closest to my low-elevation sites 

# 2. UPLMET, which the high-elevation station along on Upper Lookout 

# load in downloaded data:

precip <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/HJA_precip_data.csv")

# Format for dates and trim
precip_clean <- precip %>%
 mutate(DATE = mdy_hm(DATE)) %>%  # Convert to datetime
  mutate(YEAR = year(DATE),
         MONTH = lubridate::month(DATE, label = TRUE),
         STATION = SITECODE) %>%
  dplyr::select(DATE, YEAR, MONTH, STATION, PRECIP_TOT_DAY, PROBE_CODE)

# Check for NAs
sum(is.na(precip_clean$PRECIP_TOT_DAY))  # NA's are distributed across the probes and 
# years, so 285 isn't a biggie 

# There are two probes for the UPLMET, so want to compare to see if there's any reason 
# why I shouldn't just use PROBE01 


# Compare two UPLMET probes 
upl_compare <- precip_clean %>%
  filter(STATION == "UPLMET") %>%
  group_by(YEAR, PROBE_CODE) %>%
  summarise(annual_precip = sum(PRECIP_TOT_DAY, na.rm = TRUE)) %>%
  pivot_wider(names_from = PROBE_CODE, values_from = annual_precip)

# Plot differences
probes <- ggplot(upl_compare, aes(x = PPTUPL01, y = PPTUPL02)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed", color = "red") +
  labs(title = "UPLMET Annual Precip: Probe 01 vs Probe 02",
       x = "PPTUPL01 (mm)",
       y = "PPTUPL02 (mm)")

probes

# probe 02 has a few higher values so I'm just going to stick with probe 1 data 

precip_clean$PROBE_CODE <- as.factor(precip_clean$PROBE_CODE)

precip_clean <- precip_clean %>%
  filter(!(STATION == "UPLMET" & PROBE_CODE == "PPTUPL02"))


# calculate average normals 
annual_normals <- precip_clean %>%
  group_by(STATION, YEAR) %>%
  summarise(annual_total = sum(PRECIP_TOT_DAY, na.rm = TRUE)) %>%
  group_by(STATION) %>%
  summarise(mean_annual_precip_mm = mean(annual_total),
            sd_annual_precip_mm = sd(annual_total))

# Low elevation annual average is 2,163.0 mm 
# high elevation annual average is 2,551.4 mm 


## Also get summer precip normals 

# Make months numeric for filtering 
precip_clean <- precip_clean %>%
  mutate(MONTH_NUM = month(DATE))  

summer_precip <- precip_clean %>%
  filter(MONTH_NUM %in% c(6, 7, 8)) %>%
  group_by(STATION, YEAR) %>%
  summarise(summer_total_mm = sum(PRECIP_TOT_DAY, na.rm = TRUE)) %>%
  group_by(STATION) %>%
  summarise(mean_JJA_precip_mm = mean(summer_total_mm),
            sd_JJA_precip_mm = sd(summer_total_mm))

# Low elevation summer average is 104.7 mm 
# high elevation summer average is 169.3 mm 

