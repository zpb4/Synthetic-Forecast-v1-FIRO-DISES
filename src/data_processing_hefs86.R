#Script to process csv HEFS data from RFC

library(lubridate)
library(stringr)

print(paste('datapro start',Sys.time()))

#input main location ID
loc <- 'YRS'  #current options: 'NHG' 'YRS' 'LAM'

#---------------------Get the daily observations for each site ----------------------------
#REQUIREMENTS FOR 'observed_flows.csv'
#1) The observed flow matrix represents daily flows
#2) The first column is named "Date" and has dates formatted as yyyy-mm-dd
#3) The remaining columns each have a different site, and are named using the site ID (e.g., ADOC1)
#4) The units of flow are kcfs

obs <- read.csv(paste('./data/',loc,'/observed_flows.csv',sep=''),header=TRUE)
obs <- data.frame(obs)
ixx_obs <- as.POSIXlt(ymd(obs[,1]),tz = "UTC")  
site_names <- sort(colnames(obs)[-1])
n_sites <- length(site_names)
n_obs <- nrow(obs)

#drop the date column while keeping as a data frame
col_names <- colnames(obs)
col_names <- col_names[col_names!='Date']
col_names <- sort(col_names) #rearrange alphabetically to jive w/'folder_list' below
obs <- subset(obs,select=col_names)

#create a directory to store daily forecast .csv files for HEFS
for(sites in site_names){
  if (!dir.exists(paste('./out/',loc,'/HEFS86/',sites,sep=''))) {
    dir.create(paste('./out/',loc,'/HEFS86/',sites,sep=''),recursive=T)
  }
}
#------------------------------------------------------------------------------------------




#-------------------Process the hourly HEFS forecasts for each site -----------------------
#REQUIREMENTS FOR HEFS hindcasts
#1) All HEFS hindcasts should be stored under ./data/HEFS
#2) There should be a separate folder under ./data/HEFS for each site, and the site name should be somewhere in the title of that folder
#3) Within each site folder, there should be a set of .csv files, one for each day that a hindcast is available
#4) The date should be somewhere in the name of each file, in the format yyyymmdd (standard for HEFS output)
#5) we assume all forecasts are provided hourly, and are issued at 12 GMT
#6) the units of flow in the forecasts is kcfs
#7) the first column includes the date, and all other columns include forecasts for different ensemble members

#names of HEFS folders
folder_list <- list.files(paste('./data/',loc,'/HEFS86',sep=''))

#get metadata for the forecasts
all_files <- list.files(paste('./data/',loc,'/HEFS86/',folder_list[1],'/',sep=""))
toskip <- 0
temp <- read.csv(paste('./data/',loc,'/HEFS86/',folder_list[1],'/',all_files[1],sep=""),header=TRUE,skip=toskip)
if(!is.numeric(temp[1,2])) {
  toskip <- toskip + 1
  temp <- read.csv(paste('./data/',loc,'/HEFS86/',folder_list[1],'/',all_files[1],sep=""),header=TRUE,skip=toskip)
}
#the first row is for the current hour (12 GMT), so we drop that one
if(str_split(temp$X[1],pattern=' ')[[1]][2]=='12:00:00'){
  temp <- temp[-1,]}

#number of ensemble members        
n_ens <- ncol(temp)-1

#get the saved HEFS dates
HEFS_dates <- sort(as.numeric(gsub("\\D", "", all_files)))
start_date <- ymd(substr(HEFS_dates[1],1,8))
end_date <- ymd(substr(HEFS_dates[length(HEFS_dates)],1,8))
ixx_hefs <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date),'day'))
ixx_hefs_out <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date+15),'day'))
n_hefs <- length(ixx_hefs)

#get the number of leads (one row is for the current hour, so we drop that one)
leads <- nrow(temp)/24
leads <- min(leads,15) #only want max of 15d dynamical forecasts (some have 30 leads)

#function to convert hourly forecasts to daily average forecasts
agg_fun_dly <- function(x){out <- apply(matrix(as.numeric(x),nrow=24,byrow = F),2,mean); return(out)}


#loop through sites and create hefs_forward
hefs_forward <- array(NA,c(n_sites,n_ens,n_hefs,leads))
for(i in 1:n_sites){
  
  #get files for current site
  all_files <- list.files(paste('./data/',loc,'/HEFS86/',folder_list[i],'/',sep=""))
  
  for (j in 1:n_hefs) {
    cur_file <- grep(substr(HEFS_dates[j],1,8),all_files)
    temp <- read.csv(paste('./data/',loc,'/HEFS86/',folder_list[i],'/',all_files[cur_file],sep=""),header=TRUE,skip=toskip)
    #the first row is for the current hour (12 GMT), so we drop that one
    if(str_split(temp$X[1],pattern=' ')[[1]][2]=='12:00:00'){
      temp <- temp[-1,]}
    #convert to daily forecasts 
    hefs_inp = t(apply(temp[,-1],2,agg_fun_dly))
    hefs_forward[i,,j,] <- hefs_inp[,1:leads]
    #write and save a set of daily HEFS forecasts from hourly inputs
    hefs_dly_out<-t(hefs_forward[i,,j,])
    dly_dates<-as.character(ixx_hefs_out[(j+1):(j+dim(hefs_dly_out)[1])])
    hefs_dly_out<-cbind(dly_dates,hefs_dly_out)
    colnames(hefs_dly_out)<-c('Date',paste(col_names[i],c('',paste('.',1:(dim(hefs_dly_out)[2]-2),sep='')),sep=''))
    filename<-str_replace(all_files[cur_file],'hefs_hourly','daily')
    if(str_split(filename,'_')[[1]][2]!=col_names[i]){filename<-str_replace(filename,str_split(filename,'_')[[1]][2],col_names[i])}
    write.csv(hefs_dly_out,file=paste('./out/',loc,'/HEFS86/',col_names[i],'/',filename,sep=''),row.names = F)
  }  
}

saveRDS(hefs_forward,paste('./out/',loc,'/hefs_forward_86.rds',sep=''))

save.image(paste('./out/',loc,'/data_prep_rdata86.RData',sep=''))

print(paste('datapro end',Sys.time()))

rm(list=ls());gc()

###############################################END################################################