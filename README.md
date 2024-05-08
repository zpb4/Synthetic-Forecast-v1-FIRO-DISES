# Synthetic-Forecast-v1-FIRO-DISES
Synthetic forecast model to support FIRO work under DISES funding. Version 1 is model developed in WRR manuscript, Brodeur et al. (2024) 'Synthetic forecast ensembles for evaluating Forecast Informed Reservoir Operations'
   
---
Setup for forecast generation at Prado dam system (ADO), including main reservoir inflow (ADOC1). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/b6788237717c41e0bcc69bcaa851694f/). Starting user-defined settings are:
  
loc = 'ADO'   
keysite = 'ADOC1'   
n_samp = 10   

---
Setup for forecast generation at Lake Mendocino system (LAM), including reservoir inflow (LAMC1) and downstream local flows at Ukiah and Hopland (UKAC1, HOPC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/e51d9821c8d84682b642eb0818ac3137/). Starting user-defined settings are:
  
loc = 'LAM'   
keysite = 'LAMC1'   
n_samp = 10   

---
Setup for synthetic forecast generation at New Hogan Lake system (NHG), including reservoir inflow (NHGC1) and downstream Mud Slough site (MSGC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/dfa02b83bbde4ae3888ffafeb4446a5b/). Starting user-defined settings are:
  
loc = 'NHG'   
keysite = 'NHGC1'   
n_samp = 10   

 

---
Setup for forecast generation at selected sites of the Yuba-Feather system (YRS), including reservoir inflow at Lake Oroville (ORDC1) and New Bullards Bar (NBBC1) and downstream local flows at Marysville junction (MRYC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/29a7c696ee4e4766883078ca0d681884/). Starting user-defined settings are:
  
loc = 'YRS'   
keysite = 'ORDC1'   
n_samp = 10   

---
#### Note: After downloading and extracting data from Hydroshare resources above, ensure local directory path for HEFS data is configured: './Synthetic-Forecast-v2-FIRO-DISES/data/_main_hindcast_location_/...', where '...' are the site specific sub-repos defined in 'Data' section below. Unzipping the files can result in duplication in the data path and this must be corrected for the code to function.
   
Information below describes setup and execution of the model:   
# data

In the ./data folder there are two required sets of files. 

The first is a .csv file called 'observed_flows.csv' that contains the observed flows for all sites of interest for the entire period for which observations are available across all locations. The requirements for 'observed_flows.csv' are as follows:
1) The observed flow matrix represents daily flows
2) The first column is named "Date" and has dates formatted as yyyy-mm-dd
3) The remaining columns each have a different site, and are named using the site ID (e.g., ADOC1)
4) The units of flow are kcfs

The second set of files are located in the directory .data/HEFS/, and must conform to the following structure: 
1) There should be a separate folder under ./data/HEFS for each site, and the site name should be somewhere in the title of that folder
2) Within each site folder, there should be a set of .csv files, one for each day that a hindcast is available
3) The date should be somewhere in the name of each file, in the format yyyymmdd (standard for HEFS output)
4) we assume all forecasts are provided hourly, and are issued at 12 GMT
5) the units of flow in the forecasts is kcfs
6) the first column includes the date, and all other columns include forecasts for different ensemble members

# workflow

The user needs to run the following scripts in this order for the model to produce the synthetic forecasts:
1) ./scr/data_processing.R
2) .scr/create_synthetic_forecasts.R
3) .scr/data_writeout.R

When running create_synthetic_forecasts.R, there are a few arguments the user can specify, including the number of synthetic forecast samples to create, as well as the time period over which to create them. This script calls the function .scr/syn_gen.R, which holds the actual synthetic forecast model. 

Finally, there is a plotting script ./scr/plot_ensembles.R which can be used to visualize the results. 

