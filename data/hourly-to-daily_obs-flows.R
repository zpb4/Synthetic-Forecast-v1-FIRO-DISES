#script to convert hourly flows to appropriately configured daily 'observed_flow.csv'
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

obs_in<-read.csv('../data/2021100112_calaveras_historical_flow.csv')
date_vec<-obs_in$GMT

#remove hourly dates and infer daily date vector
dates<-strftime(date_vec[date_vec!=''],format="%Y-%m-%d")
dates_daily<-unique(dates)

agg_fun_dly <- function(x){out <- apply(matrix(as.numeric(x),nrow=24,byrow = F),2,mean); return(out)}

n_idx<-which(obs_in$GMT!='')
site_idx<-which(colnames(obs_in)!='GMT')

data_mat<-apply(obs_in[n_idx,site_idx],2,as.numeric)
data_out<-apply(data_mat,2,agg_fun_dly)

obs_flow_out<-cbind(dates_daily[-c(1)],data_out)
cnames<-colnames(obs_flow_out)
cnames[1]<-'Date'
colnames(obs_flow_out)<-cnames

write.csv(obs_flow_out,'../data/observed_flows.csv',row.names = F)

rm(list=ls());gc()

##################################END#################################################

