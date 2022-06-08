#0. setup ####
# setwd('C:/Users/Alice Carter/git/nhc_50yl/')
setwd('C:/Users/alice.carter/git/nhc_50yl')
source('src/helpers.R')
library(lubridate)
library(tidyverse)
library(zoo)
library(xts)
library(dygraphs)
# library(devtools)
# install_github('streampulse/StreamPULSE', dependencies=TRUE)
library(StreamPULSE)

  then_col = "brown3"
  now_col = "gray"
  fall_col = "brown3"

sites <- read_csv("data/siteData/NHCsite_metadata.csv") %>%
  slice(c(1:5,7))

# downloads all sensor data for every site in "sites" and saves it in raw
#1. download and prep data files ####
# a. download
source('src/download_format_data/download_rawSPdata.R')

# b. correct level data
source("src/download_format_data/correct_level_data.R")

# c. gap fill level data at NHC and UNHC
source("src/download_format_data/gap_fill_level.R")

# d. Build rating curves for NHC and UNHC
#    Calculate discharge, interpolate between sites based on watershed area
source("src/download_format_data/update_rating_curves.R")

# e. Calculate coefficients for depthxQ and velocityxQ relationships
source("src/download_format_data/calc_avg_channel_crosssections.R")

# f. get light data from StreamLight package:
source("src/download_format_data/get_light_from_StreamLight.R")

# g. Prepare datafiles for modeling metabolism
source("src/download_format_data/prep_raw_SPdata.R")

# these files are now ready to be run for SM metabolism or Hall metabolism


#2. Run Metabolism with stream metabolizer:####
# src/metabolism/run_streamMetabolizer_raymond_K600.R
# when running, make sure to manually select the sites you are running at the bottom.


#3. Compare Metabolism across datasets ####
# a. Load and compile metabolism runs from SM
  source("src/metabolism/compile_sm_models.R")

# b. Compile all model results and convert to units of C
 # Generates summary data table for paper and calculates CVs for results
 # Also plots cumulative metabolism and within year met CV across scales
  source("src/metabolism/compile_all_model_results.R")
# c. Compare metabolism estimates run with adjusted median depth to match those
  # reported in 1969:
  source('src/metabolism/compare_cbp_met_with_different_depths.R')

# c. Metabolism Plots
sites <- sites[c(1:5,7),]
  # plot daily metabolism estimates from stream metabolizer for all
  # sites in 2019 and NHC and UNHC in all years:
  source("src/plot_data/plot_met_with_Q.R")
  # generate bootstrapped CIs for direct calc estimates and plot
  # source("NHC_2019_metabolism/src/analyze_data/bootstrap_metabolism_comparison.R")
  source("src/analyze_data/bootstrap_metabolism_comparison_streamMetabolizer.R") #**
  # plot variation in Met and P/R across scales
    # updated version is for comparing Hall 72 data to modern stream metabolizer calcs
    # this file also builds the multipanel comparison across sites and years
  source("src/plot_data/plot_direct_calc_met_results.R")
  # plot temp x ER and GPP relationships
  # this script generates figures 3, 4
  source("src/analyze_data/calculate_Q10_SM_plot.R")
  # collect geomorphic driver variables, pair with metabolism, plot regressions
  # calculates slopes, r2, pvals for manuscript
    # plots from metabolism summary table - percent differences across scales
    source('src/plot_data/plot_met_summary_dat.R')


#4. Process Hall 1970 data: ####
# Calculations and data processing:
  # Calculate discharge based on rating curve and daily stage plots, pair with temperature
  source('src/analyze_data/calculate_hall_Q.R')
  # Compile discharge, depth, and K600 data from Hall 1970. Snap daily stages
    # and temperatures to match peak flow dates from table 5
    # Plot distributions of Q, K, D then and now for comparison
  source('src/analyze_data/nhc_temp_Q_k_depth_comparison.R')


################################################################################
# Other plots:
  # plot temperature and discharge comparison
  source("src/analyze_data/nhc_temp_Q_comparison.R")
  # plot land cover and terrestrial metabolism
  # source("src/plot_data/plot_riparian_NPP_LAI.R") # not this one
  # source("src/plots/plot_watershed_landcover.R")

  # summarize temperature, discharge, and DO across all metabolism site years
    # quantify hypoxia, compare drivers to met at annual timescale
    source('src/analyze_data/summarize_physical_met_drivers.R')

  # summarize hall and todays nutrient chemistry
    source('src/analyze_data/munge_nutrient_concentrations.R')
# plot multipanel climate figure
  source('src/plot_data/multi_panel_watershed_change_fig.R')


#2. Run Hall Metabolism ####
# # this section is no longer in use - see old git repo to run this
# # a. calculate k based on Churchill 1962 equation
# source("hall_50yl/code/src/analysis/calculate_k_from_hall_eq.R")
#
# # b. model metabolism following Hall 1970 method
# source("hall_50yl/code/src/analysis/model_metabolism_hall_method.R")
