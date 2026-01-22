# ---------------------------------------------------------------------------------------------#
# Process 30 year Normals (1990-2020) climate data for study sites and sub sites 
# 
# Original Author: L. McKinley Nevins 
# January 5, 2026
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     cowplot v 1.2.0
#                     lubridate v 1.9.4
#                     
# ---------------------------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
require(lubridate); packageVersion("lubridate")

###############################################################################################
#                                       Main workflow                                         #
#  Use data of monthly 30 year normalized temperature and precipitation for each site,        #
#  collected from NOAA weather stations for precip and PRISM for temp. Calculate summary      #
#  values and plot variation over the year to reflect summer drought conditions.              #                                                  #                                                                  #
#                                                                                             #
###############################################################################################

################ --
# (1) DATA PREP
################ --

wd <- "~/Dropbox/WSU/Mycorrhizae_Project/"
setwd(wd)


### NOAA DATA ### -- 

# Load in monthly NOAA data 
brookings_clim <- read.csv("./Climate_data/brookings_NOAA.csv")
carson_clim<- read.csv("./Climate_data/carson_fish_hatchery_NOAA.csv")
darrington_clim<- read.csv("./Climate_data/darrington_NOAA.csv")
dia_lake_clim<- read.csv("./Climate_data/diamond_lake_NOAA.csv")
gasquet_clim<- read.csv("./Climate_data/gasquet_NOAA.csv")
glacier_clim<- read.csv("./Climate_data/glacier_NOAA.csv")
lemolo_lake_clim<- read.csv("./Climate_data/lemolo_lake_NOAA.csv")
stampede_pass_clim<- read.csv("./Climate_data/stampede_pass_NOAA.csv")
toketee_airstrip_clim<- read.csv("./Climate_data/toketee_airstrip_NOAA.csv")

# Status: All of these have precipitation data but most of them don't have temperature, 
# so I'm going to use the PRISM data that was previously compiled to add in the temperature
# data to each file 

## ALL PRECIP FROM NOAA, ALL TEMP FROM PRISM

# Subset the NOAA datasets to just a few of interest 
brookings_clim <- dplyr::select(brookings_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
brookings_clim$Field_ID <- c("South_05")


carson_clim <- dplyr::select(carson_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
carson_clim$Field_ID <- c("WFDP")


darrington_clim <- dplyr::select(darrington_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
darrington_clim$Field_ID <- c("North_02")


dia_lake_clim <- dplyr::select(dia_lake_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
dia_lake_clim$Field_ID <- c("South_03") # Also south 04, but doesn't need to be duplicated, can 
# just note that it represents both 


gasquet_clim <- dplyr::select(gasquet_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
gasquet_clim$Field_ID <- c("South_06") # Also south 07, but doesn't need to be duplicated, can 
# just note that it represents both. I also did a little shuffling of the Field_IDs for 
# Gasquet, but it won't 


glacier_clim <- dplyr::select(glacier_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
glacier_clim$Field_ID <- c("North_01")


lemolo_lake_clim <- dplyr::select(lemolo_lake_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
lemolo_lake_clim$Field_ID <- c("South_01") 


stampede_pass_clim <- dplyr::select(stampede_pass_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
stampede_pass_clim$Field_ID <- c("North_03") 


toketee_airstrip_clim <- dplyr::select(toketee_airstrip_clim, STATION, NAME, month, MLY.PRCP.NORMAL)
# Add Field_ID column so things can be merged 
toketee_airstrip_clim$Field_ID <- c("South_02") 


# Stack all of these and convert precip from in to mm
all_NOAA_clim <- rbind(brookings_clim, carson_clim, darrington_clim, dia_lake_clim, gasquet_clim, 
                       glacier_clim, lemolo_lake_clim, stampede_pass_clim, toketee_airstrip_clim)

all_NOAA_clim <- all_NOAA_clim %>%
  dplyr::mutate(precip_mm = MLY.PRCP.NORMAL * 25.4)


### PRISM DATA ### -- 

# Load in monthly PRISM data
prism_clim <- read.csv("./Climate_data/climate_data_30yr.csv")

# This file has each of the locations stacked, but the month column is the same as the NOAA data. 
# So can separate out each location and then left join the data by month to get the variables all together in one 
# data file. 

# Subset a few columns 
prism_clim <- dplyr::select(prism_clim, Field_ID, Location, month, min_temp, mean_temp, max_temp)


# Right now these are in F so need to convert to C 
prism_clim <- prism_clim %>%
  dplyr::mutate(
    temp_max_C   = (max_temp - 32) * 5/9,
    temp_min_C   = (min_temp - 32) * 5/9,
    temp_mean_C   = (mean_temp - 32) * 5/9)


# Pull out Andrews data so it can be merged separately

# Location at a factor 
prism_clim$Location <- as.factor(prism_clim$Location)

# Pull just Andrews location 
andrews_prism_clim <- prism_clim[prism_clim$Location=="Andrews Forest, OR", ]

andrews_prism_clim$MONTH <- andrews_prism_clim$month


# Merge all_NOAA_clim with prism_clim according to Field_ID and month to get the temp and precip for 
# each site all in one place 

temp_precip <- left_join(all_NOAA_clim, prism_clim, by = c("Field_ID", "month"))

# This worked and has everything together for each site, except for Andrews 


### ANDREWS DATA ### -- 

# Load in the data for Andrews Forest, which comes directly from the two weather stations there at low
# and high elevations to get precipitation info. Will still use PRISM for temperature 

# The HJA_precip_data_wrangle.R script calculated the overall annual normals, but I need to calculate 
# the monthly 30 year normals to keep the format consistent with the rest of the sites 


# Using data from here: https://andlter.forestry.oregonstate.edu/data/abstract.aspx?dbcode=MS004

# Two meteorological stations: 
# 1. PRIMET, which is the primary station down by the offices, and the station
# that is closest to my low-elevation sites 

# 2. UPLMET, which the high-elevation station along on Upper Lookout 

# load in downloaded data:

andrews_precip <- read.csv("~/Dropbox/WSU/Mycorrhizae_Project/Climate_data/HJA_precip_data.csv")

# Format for dates and trim
andrews_precip_clean <- andrews_precip %>%
  mutate(DATE = mdy_hm(DATE)) %>%  # Convert to datetime
  mutate(YEAR = year(DATE),
         MONTH = lubridate::month(DATE, label = TRUE),
         STATION = SITECODE) %>%
  dplyr::select(DATE, YEAR, MONTH, STATION, PRECIP_TOT_DAY, PROBE_CODE)

# Add numeric month code 
andrews_precip_clean <- andrews_precip_clean %>%
  mutate(MONTH_NUM = month(DATE)) 

# Check for NAs
sum(is.na(andrews_precip_clean$PRECIP_TOT_DAY))  # NA's are distributed across the probes and 
# years, so 285 isn't a biggie 

# There are two probes for the UPLMET, so want to compare to see if there's any reason 
# why I shouldn't just use PROBE01 


# Compare two UPLMET probes 
upl_compare <- andrews_precip_clean %>%
  filter(STATION == "UPLMET") %>%
  group_by(YEAR, MONTH, PROBE_CODE) %>%
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

andrews_precip_clean$PROBE_CODE <- as.factor(andrews_precip_clean$PROBE_CODE)

andrews_precip_clean <- andrews_precip_clean %>%
  filter(!(STATION == "UPLMET" & PROBE_CODE == "PPTUPL02"))


# calculate monthly totals by year 
monthly_totals <- andrews_precip_clean %>%
  group_by(STATION, YEAR, MONTH_NUM, MONTH) %>%
  summarise(
    monthly_precip_mm = sum(PRECIP_TOT_DAY, na.rm = TRUE),
    .groups = "drop"
  )


# Calculate the monthly normals across the whole time period 
monthly_normals <- monthly_totals %>%
  group_by(STATION, MONTH_NUM, MONTH) %>%
  summarise(
    mean_monthly_precip_mm = mean(monthly_precip_mm),
    sd_monthly_precip_mm = sd(monthly_precip_mm),
    .groups = "drop"
  )


# Split the datasets and make them compatible to merge with the PRISM data 

monthly_normals$STATION <- as.factor(monthly_normals$STATION)

# Get high elevation site 
and_high_precip <- monthly_normals[monthly_normals$STATION=="UPLMET", ]

# Get low elevation site 
and_low_precip <- monthly_normals[monthly_normals$STATION=="PRIMET", ]

# Add in Field_ID column that matches the PRISM data 
and_high_precip$Field_ID <- c("Andrews_02") 
and_low_precip$Field_ID <- c("Andrews_01") 


# set month formats before merging 
and_high_precip <- and_high_precip %>%
  mutate(MONTH = as.integer(MONTH_NUM))

and_low_precip <- and_low_precip %>%
  mutate(MONTH = as.integer(MONTH_NUM))

andrews_prism_clim <- andrews_prism_clim %>%
  mutate(MONTH = as.integer(MONTH))


# Merge with the andrews_prism_clim dataset according to the Field_ID column 
and_high_final <- left_join(and_high_precip, andrews_prism_clim, by = c("Field_ID", "MONTH"))

and_low_final <- left_join(and_low_precip, andrews_prism_clim, by = c("Field_ID", "MONTH"))


###################################### ---- 

# Andrews datasets are ready to go. For the others, pull them individually from the temp_precip df 
# and then format them 

## Organizing these in order from north to south so it's easier to keep track of them 


# Glacier - North_01
only_glacier <- temp_precip[temp_precip$Field_ID=="North_01", ]

# Darrington - North_02
only_darrington <- temp_precip[temp_precip$Field_ID=="North_02", ]

# Stampede Pass - North_03 
only_stampede_pass <- temp_precip[temp_precip$Field_ID=="North_03", ]

# Carson - WFDP
only_carson <- temp_precip[temp_precip$Field_ID=="WFDP", ]

## Andrews ### 

# Lemolo Lake - South_01
only_lemolo_lake <- temp_precip[temp_precip$Field_ID=="South_01", ]

# Toketee Air Strip - South_02
only_toketee_air_strip <- temp_precip[temp_precip$Field_ID=="South_02", ]

# Diamond Lake - South_03, South_04
only_diamond_lake <- temp_precip[temp_precip$Field_ID=="South_03", ]

# Brookings - South_05
only_brookings <- temp_precip[temp_precip$Field_ID=="South_05", ]

# Gasquet - South_06, South_07
only_gasquet <- temp_precip[temp_precip$Field_ID=="South_06", ]


###############################################################################################

############################# --
# (2) PLOT MONTHLY NORMALS 
############################# --

# Growing season dates were determined from USDA hardiness zones for each site. The values entered as 
# grow_start and grow_end are bumped a little bit in the actual plotting. 


#Plotting without legends to aid in spacing 

#### North_01 - Glacier ###

# growing season - early April-late October 
grow_start <- 5
grow_end   <- 10

north_01 <- ggplot(only_glacier, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "North_01 - Canyon Creek, WA") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

north_01



#### North_02 - Darrington ###

# growing season - early April-late October 
grow_start <- 5
grow_end   <- 10

north_02 <- ggplot(only_darrington, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "North_02 - Sauk River, WA") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

north_02



#### North_03 - Stampede Pass ###

# growing season - early April-late October 
grow_start <- 5
grow_end   <- 10

north_03 <- ggplot(only_stampede_pass, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "North_03 - Stampede Pass, WA") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

north_03


#### WFDP - Carson ###

# growing season - early April-late October 
grow_start <- 5
grow_end   <- 10

WFDP <- ggplot(only_carson, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "WFDP - Wind River Forest Dynamics Plot, WA") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

WFDP


#### Andrews_01 - H.J. Andrews Experimental Forest, OR, Low Elevation ###

# growing season - late April-early November
grow_start <- 5
grow_end   <- 11

Andrews_01 <- ggplot(and_low_final, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.25, xmax = grow_end + 0.25,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = mean_monthly_precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.25, grow_end + 0.25), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "Andrews_01 - H.J. Andrews Experimental\n Forest, OR - Low Elevation") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

Andrews_01


#### Andrews_02 - H.J. Andrews Experimental Forest, OR, High Elevation ###

# growing season - late April-early November
grow_start <- 5
grow_end   <- 11

Andrews_02 <- ggplot(and_high_final, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.25, xmax = grow_end + 0.25,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = mean_monthly_precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.25, grow_end + 0.25), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "Andrews_01 - H.J. Andrews Experimental\n Forest, OR - High Elevation") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

Andrews_02



#### South_01 -Toketee Falls, OR ###

# growing season - early April- late October
grow_start <- 5
grow_end   <- 10

South_01 <- ggplot(only_lemolo_lake, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "South_01 - Toketee Falls, OR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

South_01



#### South_02 -Toketee Falls, OR ###

# growing season - early April-late October
grow_start <- 5
grow_end   <- 10

South_02 <- ggplot(only_toketee_air_strip, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "South_02 - Toketee Falls, OR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

South_02



#### South_03 -Toketee Falls, OR, Prospect, OR###

# growing season - early April-late October
grow_start <- 5
grow_end   <- 10

South_03 <- ggplot(only_diamond_lake, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.75, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.75, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "South_03 - Toketee Falls, OR; South_04 - Prospect, OR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

South_03


#### South_05 - Brookings, OR ###

# growing season - Mid March - early December
grow_start <- 4
grow_end   <- 12

South_05 <- ggplot(only_brookings, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.5, xmax = grow_end + 0.25,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.5, grow_end + 0.25), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "South_05 - Brookings, OR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

South_05


#### South_06, South_07 - Gasquet, CA ###

# growing season - early April-late October
grow_start <- 5
grow_end   <- 10

South_06 <- ggplot(only_gasquet, aes(x = month)) +
  annotate(
    "rect", xmin = grow_start - 0.25, xmax = grow_end + 0.75,
    ymin = -Inf, ymax = Inf, fill = "darkolivegreen3", alpha = 0.15) +
  geom_line(aes(y = precip_mm, color = "Precipitation"), linewidth = 1) +
  geom_line(aes(y = temp_mean_C * 5, color = "Mean Temperature"), linewidth = 1) +
  geom_line(aes(y = temp_max_C * 5, color = "Max Temperature"), linetype = "dashed") +
  geom_line(aes(y = temp_min_C * 5, color = "Min Temperature"),linetype = "dotted") +
  geom_vline(
    xintercept = c(grow_start - 0.25, grow_end + 0.75), linetype = "dashed",
    color = "darkolivegreen4", linewidth = 0.7) +
  scale_x_continuous(breaks = 1:12,
                     labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")) + 
  scale_y_continuous(name = "Precipitation (mm)", sec.axis = sec_axis(~ . / 5, name = "Temperature (°C)")) +
  scale_color_manual(
    values = c("Precipitation" = "deepskyblue3", "Mean Temperature" = "darkorange2",
               "Max Temperature" = "firebrick", "Min Temperature" = "goldenrod3")) +
  labs(x = "Month", color = "", title = "South_06; South_07 - Gasquet, CA") +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text  = element_text(size = 11),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title   = element_text(size = 11, colour = "black"))

South_06


## can rerun this plot with the legend at the bottom to get it for compiling the figures 


###############################################################################################

######################## --
# (3) ORGANIZE PLOTS 
######################## --

# Put plots in order from N to S

# north_01, north_02, north_03, WFDP, Andrews_01, Andrews_02, South_01, South_02, South_03, 
# South_05, South_06


# Too many to fit together, so splitting out to WA sites and OR and CA sites 
WA_climate_plots <- plot_grid(north_01, north_02, north_03, WFDP,
                               ncol = 2, nrow = 2, labels = c('(a)', '(b)', '(c)', '(d)'), 
                               align = "hv", hjust = -0.1)

WA_climate_plots


# Break into two to keep sizing the same 
ORCA_climate_plots_1 <- plot_grid(Andrews_01, Andrews_02, South_01, South_02,
                              ncol = 2, nrow = 3, labels = c('(a)', '(b)', '(c)', '(d)'), 
                              align = "hv", hjust = -0.1)

ORCA_climate_plots_1



ORCA_climate_plots_2 <- plot_grid(South_03, South_05, South_06,
                                  ncol = 2, nrow = 3, labels = c('(e)', '(f)', '(g'), 
                                  align = "hv", hjust = -0.1)

ORCA_climate_plots_2



# -- END -- # 
