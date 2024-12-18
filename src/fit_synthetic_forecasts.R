args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

#/////////////////////////////////////////
#Primary user defined settings

loc = 'YRS'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'
n_samp = 100            #number of samples to generate
keysite_name = 'NBBC1'   #specific site ID for 'keysite' which conditions the kNN sampling
fit_gen_strategy = 'all'  #'all' fits to all available fit data and generate across all available observations
vers<-'bvar1-sged'
cal_val_setup = 'cal'

if(cal_val_setup=='cal'){
  leave_out_years = c()
}

if(cal_val_setup=='val'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('./data/',loc,'/opt_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

if(cal_val_setup=='5fold'){
  wy = 90:119
  wy_arr = array(NA,c(5,6))
  set.seed(1)
  for(i in 1:5){
    samp = sample(wy,6,replace=F)
    wy_arr[i,] = samp + 1900
    wy = wy[!(wy%in%samp)]
  }
  leave_out_years = wy_arr[idx,]
}

print(leave_out_years)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
print(paste('start task',Sys.time()))

hefs_forward<-readRDS(paste('./out/',loc,'/hefs_forward.rds',sep=''))
load(paste('./out/',loc,'/data_prep_rdata.RData',sep=''))
source('./src/functions/cbias-hsked_functions.R')
source('./src/functions/var-fit_function.R')
source('./src/functions/sep_function.R')
source('./src/functions/knn_samp.R')

#specify start parameters
#either input 'all' to fit to all available fit data and generate across all available observations
#or input 'specify' and delineate specific dates in 'fit-start/end' and 'gen-start/end'
fit_gen_strategy<- 'all'

#date start and end for fitting; date doesn't matter if 'all' selected above
fit_start<-   '1989-10-01'       #  '1989-10-01'    
fit_end<-     '2019-09-30'       #  '2019-09-30'    
#date start for generation; date doesn't matter if 'all' selected above
gen_start<-   '1979-10-02'      #  '1979-10-02'     
gen_end<-     '2021-10-01'      #  '2021-10-01'

seasonal<-'monthly'  # 'monthly' for monthly fits; 'annual' for annual fits; o/w specify list w/monthly subsets in each list entry
cumul_samp_span<-15 #span of leads for knn sampling, set to # of leads generally

#other parameters that can be changed for experimentation
ctr<-'median'
min_quantile<-0.5
var_lag<-1
var_mean<-F
var_mod<-'bvar' # 'bvar' for BigVAR lasso penalized, 'svar' for Bigtime sparseVAR model, 'varx' for rmgarch robust least squares model
use_sigest<-T
distn<-'sged'
parallel<-T
use_mpi<-T

#index for selected keysite
keysite <- which(site_names==keysite_name)

#date/time manipulation
if(fit_gen_strategy=='all'){
  fit_rng<-ixx_hefs[ixx_hefs%in%ixx_obs]
  fit_start<-fit_rng[1]
  fit_end<-fit_rng[length(fit_rng)]
  gen_start<-ixx_obs[1]
  gen_end<-ixx_obs[length(ixx_obs)]
}

#correct fit dates to UTC convention
fit_start=as.POSIXlt(fit_start,tz='UTC')
fit_end=as.POSIXlt(fit_end,tz='UTC')

#rearrange fwd looking forecast to obs-synch
hefs<-array(NA,dim(hefs_forward))

for(i in 1:n_sites){
  hefs[i,,,]<-fwd_forecast_rearrange(hefs_forward[i,,,])
}

#synch obs and hefs to fit period
obs_fit <- obs[which(ixx_obs==fit_start):which(ixx_obs==fit_end),]
hefs_fit <- hefs[,,which(ixx_hefs==fit_start):which(ixx_hefs==fit_end),,drop=FALSE]
#-----------------------------------------
#1) kNN sample based on 'keysite' observations
print(paste('start knn',Sys.time()))

#1a. Function inputs
##knn_samp<-function(n_samp,obs,cumul_samp_span,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi)

knn<-knn_samp(n_samp,obs[,keysite],obs_idx=ixx_obs,cumul_samp_span,seasonal,fit_start,fit_end,gen_start,gen_end,leave_out_years,use_mpi)

#1b. Function outputs
knn_samps<-knn   #fitted loess model per season and K

#saveRDS(knn,paste('./out/',loc,'/knn.rds',sep=''))

print(paste('end knn',Sys.time()))
#-----------------------------------------
#2) Fit conditional expectation model
print(paste('start cexpfit',Sys.time()))
#2a. Function inputs
##fit_cexp<-function(obs,forecast,ctr,seasonal,start_date,end_date)

cexp_fit<-fit_cexp(obs_fit,hefs_fit,ctr,seasonal,fit_start,fit_end,leave_out_years,noise_reg = T)

#2b. Function outputs
loess_fit<-cexp_fit[[1]]   #fitted loess model per season and K
cexp_mat<-cexp_fit[[2]]    #nxK matrix of conditional expectation estimates
err_mat<-cexp_fit[[3]]     #exnxK matrix of errors calculated as: cexp - forecast

#saveRDS(cexp_fit,paste('./out/',loc,'/cexp-fit.rds',sep=''))

print(paste('end cexpfit',Sys.time()))

#---------------------------------------------
#3) Fit logistic regression model for forecast intermittency
print(paste('start zeromodfit',Sys.time()))
#3a. Function inputs
##fit_logreg(obs_fit,hefs_fit,seasonal,fit_start,fit_end,leave_out_years)

#lreg_fit<-fit_logreg(obs_fit,hefs_fit,seasonal,fit_start,fit_end,leave_out_years)
zero_mod_fit<-fit_zero_flow_mod(obs_fit,hefs_fit,seasonal,fit_start,fit_end,leave_out_years)

#3b. Function outputs
#lreg_coef<-lreg_fit #fitted logistic regression model coefficients per season and K

#saveRDS(zero_mod_fit,paste('./out/',loc,'/zero-mod-fit.rds',sep=''))

print(paste('end zeromodfit',Sys.time()))
#-----------------------------------------
#4) Fit linear heteroscedastic model
print(paste('start hskedfit',Sys.time()))
#4a. Function inputs
##linear_hsked_fit<-function(err_mat,cexp_mat,obs,scale_ref,seasonal,min_quantile,start_date,end_date)

hsked_fit<-linear_hsked_fit(err_mat,cexp_mat,obs_fit,scale_ref = 'cexp',seasonal,min_quantile,fit_start,fit_end,leave_out_years)

#4b. Function outputs
norm_fit<-hsked_fit[[1]]   #linear normalization model coefficients
sd_arr<-hsked_fit[[2]]     #standard deviation estimate as a backup
norm_err<-hsked_fit[[3]]   #matrix of normalized errors

#saveRDS(hsked_fit,paste('./out/',loc,'/hsked-fit.rds',sep=''))

print(paste('end hskedfit',Sys.time()))
#-----------------------------------------
#5) Fit VAR model to normalized errors
print(paste('start varfit',Sys.time()))
#5a. Funciton inputs
#var_fit<-function(err_mat,lag,use_mean,seasonal,start_date,end_date)

var_est<-var_fit(norm_err,lag=var_lag,var_mod=var_mod,use_mean = var_mean,seasonal,fit_start,fit_end,leave_out_years,parallel=T,use_mpi)

#5b. Function outputs
at_mat<-var_est[[1]]       #nxK matrix of independent residuals, a_t
var_coefs<-var_est[[2]]    #matrix of VAR coefficients by ensemble and season

saveRDS(var_est,paste('./out/',loc,'/var-fit_',vers,'.rds',sep=''))

print(paste('end varfit',Sys.time()))
#---------------------------------------------
#6) Fit SEP distribution to VAR residuals by ensemble, season, and K
print(paste('start sepfit',Sys.time()))
#6a. Function inputs
#forc_err_SEPfit<-function(err_mat,use_sigest,seasonal,start_date,end_date,parallel)

sep_est<-forc_err_SEPfit(at_mat,use_sigest,distn,seasonal,fit_start,fit_end,leave_out_years,parallel=T,use_mpi)

sep_param<-sep_est         #matrix of SEP parameter estimates by ensemble, season, and K

saveRDS(sep_est,paste('./out/',loc,'/sep-fit_',vers,'.rds',sep=''))

print(paste('end fit',Sys.time()))

if(cal_val_setup != '5fold'){
  save.image(paste('./out/',loc,'/model-fit_',vers,'_',keysite_name,'_',cal_val_setup,'.RData',sep=''))}

if(cal_val_setup == '5fold'){
  save.image(paste('./out/',loc,'/model-fit_',vers,'_',keysite_name,'_',cal_val_setup,'-',idx,'.RData',sep=''))}

#---------------------------------------------

print(paste('end task',Sys.time()))

rm(list=ls());gc()
####################################END#########################################