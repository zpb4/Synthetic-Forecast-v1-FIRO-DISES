args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

#////////////////////////////////////////////////////////////////////
loc = 'YRS'
keysite = 'NBBC1'
vers = 'bvar1-sged'
parm = 'a'
n_samp<-100
fit_gen_strategy = 'all'
cal_val_setup = 'cal'

##gen_start<-   '1979-10-02'      #  '1979-10-02'     
##gen_end<-     '2021-10-01'      #  '2021-10-01

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

library(future)
library(future.apply)
library(fGarch)
plan(multicore,workers=35)

if(cal_val_setup != '5fold'){
  load(paste('./out/',loc,'/model-fit_',vers,'_',keysite,'_',cal_val_setup,'.RData',sep=''))}
if(cal_val_setup == '5fold'){
  load(paste('./out/',loc,'/model-fit_',vers,'_',keysite,'_',cal_val_setup,'-',idx,'.RData',sep=''))}

source('./src/functions/syn_gen.R')

if(fit_gen_strategy=='all'){
  gen_start<-ixx_obs[1]
  gen_end<-ixx_obs[length(ixx_obs)]
}

nsamp_lst = as.list(1:n_samp)

#7) Generate synthetic forecasts
print(paste('start syngen',Sys.time()))

gen_start=as.POSIXlt(gen_start,tz='UTC')
gen_end=as.POSIXlt(gen_end,tz='UTC')
  
ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC")
#update VAR or AR model indicator for synthetic generation
var_ar<-'var'
if(var_mod=='ar'){var_ar<-'ar'}
  
#7a. Function inputs
#syn_gen<-function(n_samp,obs,knn_samps,obs_idx,at_mat,sep_par,cexp_fit,norm_fit,var_coefs,lag,cumul_samp_span,use_mean,seasonal,fit_start,fit_end,gen_start,gen_end,parallel)
syn_gen_out <- function(x){syn_gen(x,obs=obs,knn_samps=knn_samps,obs_idx=ixx_obs,at_mat,sep_par=sep_param,distn,cexp_fit=loess_fit,norm_fit,sd_arr,lreg_coefs=lreg_coef,
                                   zero_mod=zero_mod_fit,var_coefs,var_ar=var_ar,lag=var_lag,cumul_samp_span,use_mean=var_mean,seasonal,fit_start,fit_end,gen_start,gen_end,leave_out_years,use_mpi=TRUE)}

fut_vec <- future_lapply(nsamp_lst,syn_gen_out,future.seed=TRUE)
#fut_vec <- syn_gen_out(1)

syn_hefs_fwd<-array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

for(m in 1:n_samp){
  #rearrange flow forecast to forward looking format
  for(s in 1:n_sites){
    syn_hefs_fwd[m,s,,,]<-rearrange_to_fwd_forecast(fut_vec[[m]][s,,,])
  }
  print(paste('samp',m,Sys.time()))
}

if(cal_val_setup != '5fold'){ 
  saveRDS(syn_hefs_fwd,file=paste('./out/',loc,'/syn_hefs_forward-',parm,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))}
if(cal_val_setup == '5fold'){ 
  saveRDS(syn_hefs_fwd,file=paste('./out/',loc,'/syn_hefs_forward-',parm,'_',keysite_name,'_',cal_val_setup,'-',idx,'.rds',sep=''))}


saveRDS(ixx_gen,file=paste('./out/',loc,'/ixx_gen.rds',sep=''))
saveRDS(n_samp,file=paste('./out/',loc,'/n_samp.rds',sep=''))
  
print(paste('end syngen',Sys.time()))

rm(list=ls());gc()

##########################################END########################################