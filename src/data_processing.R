#Script to process csv HEFS data from RFC

library(lubridate)
library(stringr)

print(paste('datapro start',Sys.time()))
#---------------------Get the daily observations for each site ----------------------------
#REQUIREMENTS FOR 'observed_flows.csv'
#1) The observed flow matrix represents daily flows
#2) The first column is named "Date" and has dates formatted as yyyy-mm-dd
#3) The remaining columns each have a different site, and are named using the site ID (e.g., ADOC1)
#4) The units of flow are kcfs

obs <- read.csv('./data/observed_flows.csv',header=TRUE)
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

#create an 'out' repository to store metadata
if (!dir.exists('./out')) {
  dir.create('./out')
}

if (!dir.exists('./out/HEFS')) {
  dir.create('./out/HEFS')
}

#create a directory to store daily forecast .csv files for HEFS
for(sites in site_names){
  if (!dir.exists(paste('./out/HEFS/',sites,sep=''))) {
    dir.create(paste('./out/HEFS/',sites,sep=''))
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
folder_list <- list.files('./data/HEFS')

#get metadata for the forecasts
all_files <- list.files(paste('./data/HEFS/',folder_list[1],'/',sep=""))
toskip <- 0
temp <- read.csv(paste('./data/HEFS/',folder_list[1],'/',all_files[1],sep=""),header=TRUE,skip=toskip)
if(!is.numeric(temp[1,2])) {
  toskip <- toskip + 1
  temp <- read.csv(paste('./data/HEFS/',folder_list[1],'/',all_files[1],sep=""),header=TRUE,skip=toskip)
}
#the first row is for the current hour (12 GMT), so we drop that one
temp <- temp[-1,]

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

#function to convert hourly forecasts to daily average forecasts
agg_fun_dly <- function(x){out <- apply(matrix(as.numeric(x),nrow=24,byrow = F),2,mean); return(out)}


#loop through sites and create hefs_forward
hefs_forward <- array(NA,c(n_sites,n_ens,n_hefs,leads))
for(i in 1:n_sites){
  
  #get files for current site
  all_files <- list.files(paste('./data/HEFS/',folder_list[i],'/',sep=""))
  
  for (j in 1:n_hefs) {
    cur_file <- grep(substr(HEFS_dates[j],1,8),all_files)
    temp <- read.csv(paste('./data/HEFS/',folder_list[i],'/',all_files[cur_file],sep=""),header=TRUE,skip=toskip)
    #the first row is for the current hour (12 GMT), so we drop that one
    temp <- temp[-1,]
    #convert to daily forecasts 
    hefs_forward[i,,j,] <- t(apply(temp[,-1],2,agg_fun_dly))
    #write and save a set of daily HEFS forecasts from hourly inputs
    hefs_dly_out<-t(hefs_forward[i,,j,])
    dly_dates<-as.character(ixx_hefs_out[(j+1):(j+dim(hefs_dly_out)[1])])
    hefs_dly_out<-cbind(dly_dates,hefs_dly_out)
    colnames(hefs_dly_out)<-c('Date',paste(col_names[i],c('',paste('.',1:(dim(hefs_dly_out)[2]-2),sep='')),sep=''))
    filename<-str_replace(all_files[cur_file],'hefs_hourly','daily')
    if(str_split(filename,'_')[[1]][2]!=col_names[i]){filename<-str_replace(filename,str_split(filename,'_')[[1]][2],col_names[i])}
    write.csv(hefs_dly_out,file=paste('./out/HEFS/',col_names[i],'/',filename,sep=''),row.names = F)
  }  
}

saveRDS(hefs_forward,'./out/hefs_forward.rds')

#----------------------Cumulative stats for HEFS------------------------------------

#get forward looking cumulative totals from the hindcast
hefs_forward_cumul <- apply(hefs_forward,c(1,2,3),FUN=sum)

#also calculate the fractional values for the lead times
hefs_forward_frac <- aperm(apply(hefs_forward,c(1,2,3),function(x){x/sum(x)}),c(2,3,4,1))

#take ensemble average of cumulative totals
hefs_forward_cumul_ens_avg <- apply(hefs_forward_cumul,c(1,3),FUN=mean)

###define the residuals across ensemble members#
hefs_forward_cumul_ens_resid <- hefs_forward_cumul
for (j in 1:n_sites) {
  for(e in 1:n_ens) {
    hefs_forward_cumul_ens_resid[j,e,] <- hefs_forward_cumul[j,e,] - hefs_forward_cumul_ens_avg[j,]
  }  
}
#--------------------------------------------------------------------------------------



#----------------------Cumulative stats for obs----------------------------------------

#####prepare observed data######
#calculate cumulative observed flow totals
#these are forward looking total, not including the current day
#therefore, we cannot have a value for the first time entry. 
#we only retain these values for dates where there are a full 15 days that follow
obs_forward_all_leads <- array(NA,c(n_sites,length(1:(n_obs-leads)),leads))

for (j in 1:n_sites) {
  for(i in 1:(n_obs-leads)) {
    obs_forward_all_leads[j,i,] <- obs[(i+1):(i+leads),j]
  }
}

#these are the dates that PRECEDE the 15 day cumulative totals, i.e., 
#on date t, here is the 15 day total over the NEXT 15 days
ixx_obs_forward <- ixx_obs[1:(n_obs-leads)]
n_obs_forward <- length(ixx_obs_forward)

#subset cumul obs from to hindcast period
obs_forward_all_leads_hind <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_hefs,,drop=FALSE]
#################################


save.image("out/data_prep_rdata.RData")

print(paste('datapro end',Sys.time()))

rm(list=ls());gc()

###############################################END################################################