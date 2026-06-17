## -----------------------------------------------------------------------------#
# Project: "Taxonomic and functional composition of mycorrhizal communities respond 
# differently to host identity and environment"
#
# Process 30 year Normals (1990-2020) climate data for study sites and sub sites 
# 
# Original Author: L. McKinley Nevins 
#
# January 5, 2026
# 
# Software versions:  R v 4.4.1
#                     tidyverse v 2.0.0
#                     dplyr v 1.1.4
#                     vegan 2.6.10
#                     ggplot2 v 3.5.1
#                     cowplot v 1.2.0
#                     lubridate v 1.9.4
#                     ggord v 1.1.8
#                     ggordiplots v 0.4.3
#                     ggvegan v 0.2.1
#                     plotly v 4.11.0
#                     ggrepel v 0.9.7
#                     
# ---------------------------------------------------------------------------------------------#

# PACKAGES, SCRIPTS, AND SETUP ####
library(tidyverse); packageVersion("tidyverse")
library(dplyr); packageVersion("dplyr")
library(vegan); packageVersion("vegan")
library(ggplot2); packageVersion("ggplot2")
library(cowplot); packageVersion("cowplot")
library(lubridate); packageVersion("lubridate")
library(ggord); packageVersion("ggord")
library(ggordiplots); packageVersion("ggordiplots")
library(ggvegan); packageVersion("ggvegan")
library(plotly); packageVersion("plotly")
library(ggrepel); packageVersion("ggrepel")

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



###############################################################################################

############################################### --
# (4) GET MAIN AXES OF SITE ENVIRO VARIATION  
############################################### --

## Do PCAs for both AM and ECM hosts to get main axes of environmental variation across sites 

# Can then use the axis values for each tree in models of the effects of site environment and 
# host identity on the fungal community composition


# Load in the sample data file containing the subset of environmental data 
AM_env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/AM_enviro_all_2025.csv")

ECM_env <- read.csv(file = "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/FINAL/EM_enviro_no_THPL.csv")


# Subset environmental datasets to remove some of the soil micronutrient data that is not of interest right now 

AM_env <- dplyr::select(AM_env, Sample_ID, Site, Host_ID, Field_ID, Location, elev, mean_precip_mm, 
                        mean_summer_precip_mm, MAT, pct_N, pct_C, org_matter, ph, Sand, Silt, Clay, avg_July_SPEI, 
                        count_mod_dry, apr1_SWE)


ECM_env <- dplyr::select(ECM_env, Sample_ID, Site, Host_ID, Field_ID, Location, elev, mean_precip_mm, 
                        mean_summer_precip_mm, MAT, pct_N, pct_C, org_matter, ph, Sand, Silt, Clay, avg_July_SPEI, 
                        count_mod_dry, apr1_SWE)


# Stack these to get all locations 
all_env <- rbind(AM_env, ECM_env)


all_env <- dplyr::select(all_env, Site, Field_ID, Location, elev, mean_precip_mm, mean_summer_precip_mm, MAT, pct_N, 
                         pct_C, org_matter, ph, Sand, Silt, Clay, avg_July_SPEI, count_mod_dry, apr1_SWE)


all_env <- all_env %>% distinct(Field_ID, .keep_all = TRUE)



################ --
## ALL SITES ## 
############### -- 

#PCA of environmental variation across sites 
pca_all_enviro = prcomp(all_env[4:17], center = T, scale = T)

sd.pca_all_enviro = pca_all_enviro$sdev
loadings.pca_all_enviro = pca_all_enviro$rotation
names.pca_all_enviro = colnames(all_env[4:17])
scores.pca_all_enviro = as.data.frame(pca_all_enviro$x)
scores.pca_all_enviro$Site = all_env$Site
scores.pca_all_enviro$Field_ID = all_env$Field_ID
scores.pca_all_enviro$Location = all_env$Location
summary(pca_all_enviro)



# PCA scores are 'scores.pca_all_enviro' with column for different grouping variables 

loadings.pca_all_enviro <- as.data.frame(loadings.pca_all_enviro)

# get proportion of variance explained to add to each axis label 
pca_var <- pca_all_enviro$sdev^2  # Eigenvalues (variance of each PC)
pca_var_explained <- pca_var / sum(pca_var) * 100  # Convert to percentage


#set colors for sites
palette <- c("#580E70", "#47E5BB", "#F8BD4B", "#9D072C")

# Shapes for Field_IDs

# North_01  North_02  North_03  WFDP  Andrews_01  Andrews_02  South_01  South_02  South_03  South_04 
# South_05  South_06  South_07 
Field_ID_code <- c(15, 0, 7, 16, 17, 2, 23, 5, 9, 11, 13, 4, 8) 


# Change loadings names to something cleaner 
new_loadings <- c("Elev", "MAP", "MSP", "MAT", "Pct_N", "Pct_C", "Org_Matter", "pH", "Sand", 
                  "Silt", "Clay", "SPEI", "Mod_Dry", "SWE")
rownames(loadings.pca_all_enviro) <- new_loadings


# Plot the Results by site alone
PCA_plot_all_enviro <- ggplot(scores.pca_all_enviro, aes(x = PC1, y = PC2, color = Site)) +
  geom_point(size = 3, aes(shape = Field_ID)) +
  geom_segment(data = loadings.pca_all_enviro, aes(x = 0, y = 0, xend = PC1 * 10, yend = PC2 * 10),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") + 
  geom_text_repel(data = loadings.pca_all_enviro, aes(x = PC1 * 11, y = PC2 * 11, label = rownames(loadings.pca_all_enviro)),
                  color = "black", size = 5, max.overlaps = 10) +
  theme_minimal(base_size = 14) +
  scale_colour_manual(values=palette, 
                      name="Site",
                      breaks=c("Northern", "WFDP", "Andrews", "Southern"),
                      labels=c("Northern", "WFDP", "Andrews", "Southern")) +
  scale_shape_manual(values=Field_ID_code, 
                     name="Field ID",
                     breaks=c("North_01", "North_02", "North_03", "WFDP", "Andrews_01", "Andrews_02", "South_01",
                              "South_02", "South_03", "South_04", "South_05", "South_06", "South_07"),
                     labels=c("North_01", "North_02", "North_03", "WFDP", "Andrews_01", "Andrews_02", "South_01",
                              "South_02", "South_03", "South_04", "South_05", "South_06", "South_07")) +
  labs(x = paste0("PC1 (", round(pca_var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(pca_var_explained[2], 1), "%)")) +
  theme(axis.line = element_line(color = "black", linewidth = 0.75, linetype = "solid")) +
  theme(legend.title = element_text(colour="black", size=14, face="bold")) +
  theme(legend.text = element_text(colour="black", size = 14)) +
  theme(axis.text.x = element_text(colour="black", size = 14),
        axis.text.y = element_text(colour="black", size = 14)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)) + 
  theme(legend.position = "right")

PCA_plot_all_enviro



##Broken-Stick test for the significance of the loadings
print(pca_all_enviro)

plot(pca_all_enviro, type = "l")


ev = pca_all_enviro$sdev^2

evplot = function(ev) {
  # Broken stick model (MacArthur 1957)
  n = length(ev)
  bsm = data.frame(j=seq(1:n), p=0)
  bsm$p[1] = 1/n
  for (i in 2:n) bsm$p[i] = bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p = 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op = par(mfrow=c(2,1),omi=c(0.1,0.3,0.1,0.1), mar=c(1, 1, 1, 1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

evplot(ev)

## Screenshots of these outputs were saved 



summary(pca_all_enviro)

#can safely retain first two - PC1, PC2, and maybe the third one too 

# figure out which environmental factors were most strongly associated with the first two PC's

# Get top 4 variables for PC1
top_PC1 <- loadings.pca_all_enviro[order(abs(loadings.pca_all_enviro$PC1), decreasing = TRUE), ][1:4, ]

# Get top 4 traits for PC2
top_PC2 <- loadings.pca_all_enviro[order(abs(loadings.pca_all_enviro$PC2), decreasing = TRUE), ][1:4, ]

# Get top 4 traits for PC3
top_PC3 <- loadings.pca_all_enviro[order(abs(loadings.pca_all_enviro$PC3), decreasing = TRUE), ][1:4, ]


# PC1 is showing a spread according to Mod_dry, MAP, Sand, and Silt. Essentially precipitation and then
# soil texture. This spread is really mainly between a few of the southern sites.  

# PC2 is showing a spread according to Pct_N, Pct_C, Org_Matter and SPEI. Essentially soil properties and drought conditions. 
# This one is very strongly driven by WFDP with low values, which was the site with the highest soil C and organic matter 
# content, and then South_06 which was the ALRU riverbed site that had essentially no organic matter in the soil. 

# PC3 is showing a spread according to SWE, MSP, Clay, and MAT. So climate, snow conditions, and then clay content 
# that can also relate to water holding in the soil. 

# Grab PC1, PC2, and PC3 for each Field_ID
all_enviro_PCs <- dplyr::select(scores.pca_all_enviro, PC1, PC2, PC3, Site, Field_ID, Location)

# Save this file for later 
write.csv(all_enviro_PCs, "~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/PCA_Outputs/all_enviro_PCs.csv")


###############################################################################################

################################ --
# (5) ORGANIZE ENVIRO PCA PLOTS
################################ --

# Save all environmental factors PCA
ggsave("~/Dropbox/WSU/Mycorrhizae_Project/Community_Analyses/PCA_Outputs/all_enviro_PCA.png", 
       plot = PCA_plot_all_enviro, width = 7.5, height = 6.5, units = "in", dpi = 300)



# -- END -- # 
