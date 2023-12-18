
rm(list=ls())
library(lubridate)

print(paste('data-write start',Sys.time()))
#----------------------------------------

load("out/data_prep_rdata.RData")
syn_hefs_forward <- readRDS('out/syn_hefs_forward.rds')
ixx_sim <- readRDS('out/ixx_sim.rds') 
n_samp <- readRDS('out/n_samp.rds') 

#get rid of all previously saved text files
#CAUTION: dont do this unless you're sure
all_out <- list.files('./out')
folders_to_delete <- all_out[grep('syn_samp_out',all_out)]
setwd('./out')
for (h in 1:length(folders_to_delete)) {
  unlink(folders_to_delete[h],recursive=TRUE)
}
setwd('../')

#write out the data by synthetic sample and site
for (i in 1:n_samp) {
  #create folder for each sample if needed
  cur_dir <- paste("./out/syn_samp_out",i,sep='')
  if (!dir.exists(cur_dir)) {
    dir.create(cur_dir)
  }
  
  for (j in 1:n_sites) {
    #create folder for each site under each sample if needed
    cur_dir_site <- paste(cur_dir,"/",site_names[j],sep='')
    if (!dir.exists(cur_dir_site)) {
      dir.create(cur_dir_site)
    }
    
    #loop through dates
    for (t in 1:length(ixx_sim)) {
      cur_syn_hefs_forward <- t(syn_hefs_forward[i,j,,t,])
      tt <- paste(gsub("-","",ymd(ixx_sim[t])),site_names[j],"daily",sep="_")
      tt <- paste(cur_dir_site,"/",tt,'.csv',sep='')
      forecast_dates <- seq((ixx_sim[t]+60*60*24),(ixx_sim[t]+60*60*24*leads),by="day") #add number of seconds in a day to get to the next day
      col_names <- c("Date",rep(site_names[j],ncol(cur_syn_hefs_forward)))
      cur_syn_hefs_forward_final <- data.frame(forecast_dates,cur_syn_hefs_forward)
      names(cur_syn_hefs_forward_final) <- col_names
      write.csv(cur_syn_hefs_forward_final,file=tt,row.names = FALSE)
    }
  }  
}

print(paste('data-write end',Sys.time()))
