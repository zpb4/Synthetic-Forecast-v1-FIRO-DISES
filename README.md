# Synthetic-Forecasts-v1-FIRO-DISES
Synthetic forecast model to support FIRO work under DISES funding. Version 1 is model developed in WRR manuscript, Brodeur et al. (2023) 'Synthetic forecast ensembles for evaluating Forecast Informed Reservoir Operations'

- 'main' branch holds the basic model including observed data for the New Hogan site. HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/f63ead2d62414940a7d90acdc234a5d1/) and must be extracted to 'data' sub-repo to run the model.  
   
- 'NHG' branch is full implementation of model for New Hogan reservoir inflow (NHGC1) and downstream local flows at Mud Slough site (MSGC1L). HEFS data is stored on a zip file [here](https://www.hydroshare.org/resource/f63ead2d62414940a7d90acdc234a5d1/) and must be extracted to 'data' sub-repo to run the model.   

   
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

