
rm(list=ls())
library(lubridate)
library(ncdf4)


#----------------------------------------
loc <- 'YRS'
site <- 'NBBC1'
parm <- 'a'
cal_val_setup <- 'cal'

load(paste('./out/',loc,'/data_prep_rdata.RData',sep=''))

#remove previous files if existing (netcdf will not overwrite files)
unlink(paste('./out/',loc,'/Qf-hefs.nc',sep=''),recursive=TRUE)
unlink(paste('./out/',loc,'/Qf-syn',parm,'_',cal_val_setup,'.nc',sep=''),recursive=TRUE)

syn_hefs_forward <- readRDS(paste('./out/',loc,'/syn_hefs_forward-',parm,'_',site,'_',cal_val_setup,'.rds',sep=''))
ixx_gen <- readRDS(paste('./out/',loc,'/ixx_gen.rds',sep='')) 
n_samp <- readRDS(paste('./out/',loc,'/n_samp.rds',sep='')) 
#syn_hefs_forward <- syn_hefs_forward[1:10,,,,]

#add a single entry dimension to match synthetic forecasts
hefs_fwd<-array(NA,c(1,dim(hefs_forward)))
hefs_fwd[1,,,,]<-hefs_forward

#rearrange dimensions to match NHG model in Python
hefs_out<-aperm(hefs_fwd,c(1,2,5,3,4))

#---------------------------------------------------
#hefs
#define the dimensions
ens_dim<-ncdim_def('ensemble','',0:(dim(hefs_fwd)[1]-1))
site_dim<-ncdim_def('site','',0:(dim(hefs_fwd)[2]-1))
trace_dim<-ncdim_def('trace','',0:(dim(hefs_fwd)[3]-1))
date_dim<-ncdim_def('date',paste('days since',as.character(ixx_hefs[1]),'00:00:00'),0:(dim(hefs_fwd)[4]-1),calendar = 'proleptic_gregorian')
ld_dim<-ncdim_def('lead','',0:(dim(hefs_fwd)[5]-1))

#write the variable to the netcdf file and save
hefs_var<-ncvar_def('hefs','kcfs',dim=list(ens_dim,site_dim,ld_dim,trace_dim,date_dim))
hefs_nc<-nc_create(paste('./out/',loc,'/Qf-hefs.nc',sep=''),hefs_var,force_v4 = F)
ncvar_put(hefs_nc,hefs_var,hefs_out)
nc_close(hefs_nc)

#-------------------------------------------------------------
#syn-hefs
shefs_out<-aperm(syn_hefs_forward,c(1,2,5,3,4))
#define the dimensions
ens_dim<-ncdim_def('ensemble','',0:(dim(syn_hefs_forward)[1]-1))
site_dim<-ncdim_def('site','',0:(dim(syn_hefs_forward)[2]-1))
trace_dim<-ncdim_def('trace','',0:(dim(syn_hefs_forward)[3]-1))
date_dim<-ncdim_def('date',paste('days since',as.character(ixx_gen[1]),'00:00:00'),0:(dim(syn_hefs_forward)[4]-1),calendar = 'proleptic_gregorian')
ld_dim<-ncdim_def('lead','',0:(dim(syn_hefs_forward)[5]-1))

#write the variable to the netcdf file and save
shefs_var<-ncvar_def('syn','kcfs',dim=list(ens_dim,site_dim,ld_dim,trace_dim,date_dim))
shefs_nc<-nc_create(paste('./out/',loc,'/Qf-syn',parm,'_',site,'_',cal_val_setup,'.nc',sep=''),shefs_var,force_v4 = F)
ncvar_put(shefs_nc,shefs_var,shefs_out)
nc_close(shefs_nc)

rm(list=ls());gc()

#################################END#############################################

