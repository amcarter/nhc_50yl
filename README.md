# Metabolism in New Hope Creek 50 years after the first measurement

This Project contains the data and code to model the metabolism across 3 years at 6 sites in New Hope Creek, NC, the site of the first annual metabolism study. There is also code for downloading climate drivers, quantifying change over time, and visualizing the modern and historical data, as well as building a time series model of metabolism to hindcast predictions to 1969.

Code and data in this repository were generated by Alice M Carter, Michael J Vlah, and Emily S Bernhardt.


**Contents**
  
1. [Data Sets](#data-sets-description)  
    - [Site Data](#site-data)  
    - [Longitudinal Sampling](#longitudinal-sampling)  
    - [Watershed Data](#watershed-data)  
    - [Hall 1970 Data](#hall-1970-data)  
    - [Timeseries Model Fits](#timeseries-model-fits) 
2. [Workflow](#workflow)  
    - [Download and Format Data](#download-and-format-data)  
    - [Metabolism Model](#metabolism-model)  
    - [Data Analysis](#data-analysis)  
    - [Plot Scripts](#plot-scripts)
    - [Timeseries Model](#timeseries-model) 

<!-- Data Sets description -->
## Data Sets Description

<!-- Site Data -->
## Site Data

**1. NHC_site_metadata.csv** 
    - metadata for sites included in the modern study including latitude, longitude, sensor offsets from the stream bed, upstream geomorphic characteristics, watershed areas (ws_area.km2), percent impervious surface cover (isc.percent), comid and slope from the national hydrography dataset. 

**2. all_nhc_ysi_data.csv** 
    - Manually collected data from sensor visits 

**3. NOAA_airpres.csv**
    - Air pressure data downloaded from the National Oceanic and Atmospheric Administration based on the latitude and longitude of the study reach outlet. 

<!-- Rating Curves-->
## Rating Curves

**1. discharge_measurements.csv**
    - Manually measured discharge and velocities from all sites.

**2. discharge_measurements_NHC_UNHC.csv**
    - data from streampulse paired with discharge measurements for the upstream and downstream sites

**3. NHC_UNHC_Q.csv**
    - discharge measurements for the upstream and downstream sites

**4. flow_dates_filter.csv**
    - dates on which the flow was above the rating curves or that changed by too large a magnitude within a single day and violates the metabolism model assumptions and therefore that day needs to be excluded.

**5. interpolatedQ_allsites.csv**
    - daily discharge at all sites interpolated from UNHC and NHC based on watershed area and stream length
 
<!-- Longitudinal Sampling -->
## Longitudinal Sampling
NOTE: most of this data not used in this study.

**1. nhc_channel_crosssections.csv**
    - Compiled data from geomorphic surveys in the 1 km upstream of the upstream and downstream sampling sites

**2. nhc_habitat_survey_20190225.csv**
    - Survey of habitat changes along the entire study reach on 25 Feb 2019.

**3. nhc_habitat_survey_20190307.csv**
    - Survey of habitat changes along the entire study reach on 7 March 2019.

**4. NHCLongitudinalDO_20190308.csv**
    - Longitudinal samples of water chemistry and widths and thalweg depths along the entire study reach on 8 March 2019.

**5. NHCLongitudinalDO_20191009.csv**
    - Longitudinal samples of water chemistry and widths and thalweg depths along the entire study reach on 9 Oct 2019.

<!-- Watershed Data-->
## Watershed Data

**1. nhcblands_usgs_stats.csv**
    - data and calculated metrics (Richards Baker Index (RBI), autoregressive coefficient (ar_1), and 5th quantile of discharge (q05)) for discharge at the Eno River USGS gauge, from a basin neighboring NHC.

**2. nlcd_1992-2016_summary.csv**
    - Data from the National Land Cover Dataset on total land area (CellTally) in the NHC watershed in each land use category during each survey year.

**3. nldas.csv**
    - Data on precipitation from NLDAS summarized over the watershed area

**4. noaa_air_temp.csv**
    - Daily air temperature data from the National Oceanic and Atmospheric Administration based on the latitude and longitude of the study reach outlet.

**5. prism_raw.csv**
    - Precipitation data from PRISM.

<!-- Hall 1970 data -->
## Hall 1970 Data

Data from tables and digitized from plots in Hall 1970 and 1972.

**1. hall_discharge_temp_daily.csv**
    - Temperature and discharge compiled from digitized figures in Hall 1970 and 1972.

**2. hall_table_5_flood_stages_cbp.csv** 
    - Data table 5 from Hall 1970

**3. hall_table_13_p.csv**
    - Data table 13 from Hall 1970

**4. hall_table_14_nitrogen.csv**
    - Data table 14 from Hall 1970

**5. hall_table_15.csv**
    - Data table 15 from Hall 1970

**6. hall_table_26.csv**
    - Data table 26 from Hall 1970

**7. hall_tableA2_k_morphology.csv**
    - Data table A2 from Hall 1970

**8. hall_tableA2_k_morphology_extended.csv** 
    - Data table A2 from Hall 1970 extended to include calculations of gas exchange in other units

**9. hall_figure5_digitized_ratingcurve.csv**
    - Data digitized from figure 5 from Hall 1970

**10. hall_figure26_digitized_dailystage.csv**
    - Data digitized from figure 26 from Hall 1970

**11. hall_figure27_digitized_mean_daily_temp.csv**
    - Data digitized from figure 27 from Hall 1970

<!-- Timeseries Model Fits -->
## Time Series Model Fits

**1. Arima_hindcast_models_coefficients.csv**
    - coefficients from Arima models fit on GPP and ER
    
**2. BRMS_hindcast_models_coefficients.csv**
    - coefficients from BRMS models fit on GPP and ER


<!-- Workflow -->
## Workflow

**1. master.R**
    - single script outlining the steps and with code to run some of the files to work through the data prep and analysis in this project
    
**2. helpers.R**
    - Script with custom functions used throughout this project.


<!-- download and format data -->
### Download and Format Data

**1. download_rawSPdata.R**  
    - Download raw oxygen data from publicly available StreamPulse portal.
    
**2. correct_level_data.R**
    - Downloads raw air pressure data from NOAA based on the latitude and longitude of sites. Uses NOAA air pressure data paired with manual depth measurements collected at each sensor visit (all_nhc_ysi_data.csv) to convert water pressure data into stream level at each of the sites. Sensor data must be downloaded into the appropriate folder first by running the above script. Note, this file must first be run with the sites dataframe containing only NHC (row 1) to generate the corrected NHC level file which is then used to correct the other sites, by re running with all sites in the sites dataframe.

**3. gap_fill_level.R**
    - gap fills level data at the upstream and downstream sites based on data from the other sensor when one has missing data. 

**4. update_rating_curves.R**
    - Calculates discharge rating curves for the upstream and downstream sites based on field measurements. 

**5. calc_avg_channel_crosssections.R** 
    - Calculates average reach widths and depths in the 1km reach upstream of each site based on data from geomorphic surveys and on longitudinal transects of the entire study reach. Data referenced are contained in the longitudinal_sampling data folder.

**6. calc_discharge_from_crosssection.R**
    - Calculates a discharge based on widths, depths, and velocities collected using hand profiling with an electromagnetic sensor.

**7. prep_raw_SPdata.R**
    - Formats raw sensor data to be ready to run using StreamMetabolizer. Calculates average depths, incident light, discharge, and formats datetimes as solar times. 

**8. compare_q_to_usgs_site.R**
    - downloads the full USGS discharge record from the site on lower NHC to compare flows to those measured at our upstream sites.
    
**9. compare_sensor_to_handheld_DO.R**
    - plot the deviations of sensor DO and temperature to measurements taken by hand during site visits with the YSI.
    
**10. get_light_from_StreamLight.R**
    - download light data from NLDAS and APPEARS databases and process using the streamlight package to estimate stream surface light at all sites.

<!-- metabolism model -->
### Metabolism Model

**1. run_streamMetabolizer_raymond_K600.R**
    - Runs stream metabolizer models on the prepared datafiles using a Bayesian partially pooled model with both observation and process error. Models are run on individual site years of data.
 
**2. inspect_model_fits.R**
    - Functions for evaluating, reformatting, and plotting stream metabolizer output objects. This script can be sourced at the beginning of a file and doesn’t need to be opened.

**3. compile_sm_models.R**
    - Compiles and filters all model output objects removing bad data fits, out of bounds solutions and days with flow conditions that violate model assumptions. Plots model fit information for the supplementary materials.

**4. compile_all_model_results.R**
    - Compiles metabolism estimates from stream metabolizer with reported values from Hall 1970 and 1972. Converts all data into units of g C from g O2. Calculates summary metrics including peak metabolic windows and cumulative annual metabolism. Generates data included in Table 1.

**5. compare_cbp_met_with_different_depths.R**
    - Compares the metabolism estimates from site NHC_5 calculated with modern depths to those re-estimated using adjusted depths to match those reported in Hall 1970. 

**6. compare_pooled_vs_individual_nhc_met_years.R**
    - look at differences in model runs on individual site years vs pooling across all years at a single site. 
    - plots a comparison figure for the SI
    
**7. run_metab_script.sh**
    - shell script to run metabolism models on an HPC
    
<!-- data analysis -->
### data analysis

**1. bootstrap_metabolism_comparison_streamMetabolizer.R**
    - Calculates bootstrapped 95% confidence intervals for seasonally weighted data from this study and from Hall 1972. Plots seasonal distribution of sample dates, seasonal mean rates, and bootstrapped means.
 
**2. calculate_Q10_SM_plot.R**
    - Calculates relationships between metabolic rates and drivers in the modern dataset. Plots figures 3, 4 and generates statistics for results section.

**3. calculate_hall_Q.R**
    - uses digitized data from figures in Hall 1970 to reconstruct daily discharge and temperature from 1969.

**4. nhc_temp_Q_k_depth_comparison.R**
    - Compares distributions of measured depth, discharge and estimated K600 from this study and from Hall 1970, 1972. Plots supplementary figure 2.

**5. nhc_temp_Q_comparison.R**
    - Compares daily records of temperature and discharge at NHC from 1969 and 2019. Generates figure 6.

**6. summarize_physical_met_drivers.R**
    - Compiles and summarizes drivers of metabolic patterns in the modern dataset. 

**7. munge_nutrient_concentrations.R**
    - Compares nutrients from Hall 1970, 1972 to today’s measured concentrations. Calculates regressions for DOC and nitrate as drivers of metabolism.

**8. get_spatial_predictors.R**
    - access GEE, delineates watersheds - may not be functional and isn't used in this project anymore.

<!-- Plot Scripts -->
### Plot Scripts

**1. plot_met_with_Q.R**
    - Plots daily metabolism with discharge for all site years in the modern dataset. Generates figure 2

**2. plot_direct_calc_met_results.R**
    - plots monthly averaged metabolism data for the modern and Hall datasets as well as comparisons of metabolic patterns across sites and years. Generates figure 8. 

**3. plot_met_summary_dat.R**
    - Compares geomorphic driver variables with metabolism estimates. Calculates statistics for manuscript and plots regressions.

**4. multi_panel_watershed_change_fig.R**
    - plots data for watershed and climate change. Generates figure 5.
    
**5. plot_hall_now_metab_excessER.R**
    - plots the modern and historical metabolism in a plot that matches Hall 1972 for the CB site.


    
<!-- Timeseries Model -->
### Timeseries Model

**1. ARIMA_BRMS_models_of_met.R**
    - timeseries models using 1) the Arima function in the forecast package then 2)BRMS models to hindcast metabolism and make plots.

**2. trial_mods_met_of_temp.R**
    - trial timeseries models to hindcast.
    
**3. fit_stan_models.R**
    - prep data and run stan code to fit time series models and make a forecast 
    
**4. ar1_ER.stan**
    - Stan code to run an ar1 timeseries model of respiration

**5. ar1_GPP.stan**
    - Stan code to run an ar1 timeseries model of GPP.

**6. lme_met_of_temp_Q.R** 
    - early modeling attempt. No longer in use.
    