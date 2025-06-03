# Synthetic-Forecast-v1-FIRO-DISES
Synthetic forecast model to support FIRO work under DISES funding. Version 1 is model developed in WRR manuscript, Brodeur et al. (2023) 'Synthetic forecast ensembles for evaluating Forecast Informed Reservoir Operations'.  This method is referred to as 'syn-M1' in the manuscript entitled 'Synthetic ensemble forecasts: operations-based evaluation and inter-model comparison for reservoir systems across California' published in Water Resources Research in 2025.

The model is setup generically to run synthetic forecasts for any main hindcast location with an arbitrary number of associated sites. The main location is specified by the 'loc' variable at the top of all main workflow scripts, which will point to the location specific subdirectory (i.e. ./my_directory/_main_hindcast_location_) in the 'data' sub-repo once the data have been downloaded from the referenced Hydroshare resources. The user must also specify the 'keysite' variable as the primary site to condition the sampling procedure. Typically, this would be the main reservoir inflow point for a smaller system. For a system with multiple reservoir inflow points, it is up to user discretion, but one strategy is to use the largest (by annual inflow magnitude) inflow point. Recommended settings for these primary user defined variables are indicated below. The scripts will create and write to an 'out' sub-repo with a similar structure to the 'data' sub-repo. The scripts are setup, by default, to fit the model to all available hindcast data with some leftout water years for validation and generate a user-defined number of samples across the entire observational record. These time data ranges can be modified by the user.

---
Setup for forecast generation at Prado dam system (ADO), including main reservoir inflow (ADOC1). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/b6788237717c41e0bcc69bcaa851694f/). Starting user-defined settings are:
  
loc = 'ADO'   
keysite = 'ADOC1'   
n_samp = 100   

Brodeur, Z. P. (2025). ADO Synthetic Forecast Generation Dataset, HydroShare, https://doi.org/10.4211/hs.b6788237717c41e0bcc69bcaa851694f   

---
Setup for forecast generation at Lake Mendocino system (LAM), including reservoir inflow (LAMC1) and downstream local flows at Ukiah and Hopland (UKAC1, HOPC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/e51d9821c8d84682b642eb0818ac3137/). Starting user-defined settings are:
  
loc = 'LAM'   
keysite = 'LAMC1'   
n_samp = 100   

Brodeur, Z. P. (2025). LAM Synthetic Forecast Generation Dataset, HydroShare, https://doi.org/10.4211/hs.e51d9821c8d84682b642eb0818ac3137    

---
Setup for synthetic forecast generation at New Hogan Lake system (NHG), including reservoir inflow (NHGC1) and downstream Mud Slough site (MSGC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/dfa02b83bbde4ae3888ffafeb4446a5b/). Starting user-defined settings are:
  
loc = 'NHG'   
keysite = 'NHGC1'   
n_samp = 100   

Brodeur, Z. P. (2025). NHG Synthetic Forecast generation dataset, HydroShare, https://doi.org/10.4211/hs.dfa02b83bbde4ae3888ffafeb4446a5b   

---
Setup for forecast generation at selected sites of the Yuba-Feather system (YRS), including reservoir inflow at Lake Oroville (ORDC1) and New Bullards Bar (NBBC1) and downstream local flows at Marysville junction (MRYC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/29a7c696ee4e4766883078ca0d681884/). Starting user-defined settings are:
  
loc = 'YRS'   
keysite = 'ORDC1'  and 'NBBC1' 
n_samp = 100   

Brodeur, Z. P. (2025). YRS Synthetic Forecast Generation Dataset, HydroShare, https://doi.org/10.4211/hs.29a7c696ee4e4766883078ca0d681884    

---
#### Note: After downloading and extracting data from Hydroshare resources above, ensure local directory path for HEFS data is configured: './Synthetic-Forecast-v2-FIRO-DISES/data/_main_hindcast_location_/...', where '...' are the site specific sub-repos defined in 'Data' section below. Unzipping the files can result in duplication in the data path and this must be corrected for the code to function.

## Dependencies
- R package 'BigVAR'
- R package 'stringr'
- R package 'lubridate'

   
Information below describes setup and execution of the model:   
## Data
Data that is downloaded from the Hydroshare resources above is already processed to the requirements below, assuming the 'Note' at the end is heeded. These datasets contain scripts used to process the raw HEFS hindcast files, which may be useful for pre-processing data for a new site to match the specifications detailed below.

In the ./data/_main_hindcast_location_ folder there are two required sets of files. 

The first is a .csv file called 'observed_flows.csv' that contains the observed flows for all sites of interest for the entire period for which observations are available across all locations. The requirements for 'observed_flows.csv' are as follows:
1) The observed flow matrix represents daily flows
2) The first column is named "Date" and has dates formatted as yyyy-mm-dd
3) The remaining columns each have a different site, and are named using the site ID (e.g., NHGC1)
4) The units of flow are kcfs

The second set of files are located in the directory ./data/_main_hindcast_location_/HEFS/, and must conform to the following structure: 
1) There should be a separate folder under ./data/_main_hindcast_location_/HEFS for each site, and the site name should be somewhere in the title of that folder
2) Within each site folder, there should be a set of .csv files, one for each day that a hindcast is available
3) The date should be somewhere in the name of each file, in the format yyyymmdd (standard for HEFS output)
4) we assume all forecasts are provided hourly, and are issued at 12 GMT
5) the units of flow in the forecasts is kcfs
6) the first column includes the date, and all other columns include forecasts for different ensemble members

Note: Most locations contain a smaller (2 month) set of HEFS data specific to a February 1986 flood event. This must be configured exactly as above with the directory structure: ./data/_main_hindcast_location_/HEFS86/   
Note2: The data files for each location contain a 6 year subset of validation years that is, by default, left out of the synthetic forecast fitting data/procedure. These left out years are generated by the './src/val_yrs_select.R' script to attempt maximize annual flow distributional similarity between the CAL and VAL data for each location/keysite.

## Workflow
### 1. Data Processing 
Note: Specification of the 'loc' variable in these scripts will process all associated sites
Processes the forecast and observation data to standard R data structures for Period of Record (POR) hindcast:
1) ./src/data_processing.R
   - all HEFS and observations data files for POR
   
### 2. Fitting
Note: Must specify the 'loc' variable for the overall location and the 'keysite_name' variable for the site used to condition the kNN sampling. Other user defined parameters are set to the default used in the study, including '5fold-test' for the fully out-of-sample 5-fold procedure outlined in the manuscript.  
Fits the multivariate statistical model to the POR hindcasts and observations:
1) ./src/fit_synthetic_forecasts.R
   - main fitting script; calls multiple functions from the ./src/functions sub-repo

### 3. Generation
Note: Must specify the 'loc' variable for the overall location and the 'keysite_name' variable for the site used to condition the kNN sampling. Other user defined parameters are set to the default used in the study, including '5fold-test' for the fully out-of-sample 5-fold procedure 
The user needs to run the following scripts in this order for the model to produce the synthetic forecasts. _Note that the ./src/generate_synthetic_forecasts.R script calls the function ./src/functions/syn_gen.R, which holds the actual synthetic forecast model_:
1) ./src/generate_synthetic_forecasts.R

The output of the first two steps is an R array that is saved as an R data structure file (.rds). In order to further post-process data for transfer to other models, languages, etc, there are two output options:   

2) ./src/5fold-val_collate.R
   - combines the folds to a single array for the 5-fold generation schemes ('5fold', '5fold-test')
3) ./src/data_writeout.R
   - writes individual .csv files in the same format as the input HEFS .csv files for each generated sample
4) ./src/data_writeout_ncdf.R
   - writes both HEFS and synthetic HEFS files to a netCDF file
5) ./src/slice_plot-ens.R
   - slices a 10x sample subset from the generated synthetic forecast array for plotting; the raw arrays for large sample runs (e.g. 100 samples) require too much RAM for typical personal computers

All scripts create and output metadata to the ./out/_main_hindcast_location_/ subdirectory. For sites with separate 1986 data, there are separate scripts with a '_86.R' suffix to process those specific data subsets.  

## 4. Verification

Detailed forecast verification for the synthetic forecasts is located at this repo: [https://github.com/zpb4/Synthetic-Forecast_Verification](https://github.com/zpb4/Synthetic-Forecast_Verification). This repo is designed to integrate seamlessly with the Synthetic Forecast generation output if located on the same root directory.


## Contact

Zach Brodeur: zpb4@cornell.edu
