
print(paste('start task',Sys.time()))

hefs_forward<-readRDS('./out/hefs_forward.rds')
load("./out/data_prep_rdata.RData")
source('./src/functions/cbias-hsked_functions.R')
source('./src/functions/var-fit_function.R')
source('./src/functions/sep_function.R')
source('./src/functions/syn_gen.R')
source('./src/functions/knn_samp.R')

#specify start parameters
#either input 'all' to fit to all available fit data and generate across all available observations
#or input 'specify' and delineate specific dates in 'fit-start/end' and 'gen-start/end'
fit_gen_strategy<- 'all'
generate_syn_forecasts<-T  #if 'T', generates forecasts, 'F' saves image for generation in separate routine (may be useful for large generation runs)

#date start and end for fitting; date doesn't matter if 'all' selected above
fit_start<-   '1989-10-01'       #  '1989-10-01'    
fit_end<-     '2019-09-30'       #  '2019-09-30'    
#date start for generation; date doesn't matter if 'all' selected above
gen_start<-   '1979-10-02'      #  '1979-10-02'     
fit_end<-     '2021-10-01'      #  '2021-10-01'

leave_out_years <- c(1995,2000,2005,2010,2015)

keysite<-'NHGC1'  #select site for kNN sampling
seasonal<-'monthly'  # 'monthly' for monthly fits; 'annual' for annual fits; o/w specify list w/monthly subsets in each list entry
n_samp<-100
vers<-'bvar1-sged'
cumul_samp_span<-15 #span of leads for knn sampling, set to # of leads generally

#other parameters that can be changed for experimentation
ctr<-'median'
min_quantile<-0.1
var_lag<-1
var_mean<-F
var_mod<-'bvar' # 'bvar' for BigVAR lasso penalized, 'svar' for Bigtime sparseVAR model, 'varx' for rmgarch robust least squares model
use_sigest<-T
distn<-'sged'
parallel<-T
use_mpi<-T

#date/time manipulation
if(fit_gen_strategy=='all'){
  fit_start<-ixx_hefs[1]
  fit_end<-ixx_hefs[length(ixx_hefs)]
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
hefs_fit <- hefs[,,which(ixx_hefs==fit_start):which(ixx_hefs==fit_end),]

#-----------------------------------------
#1) kNN sample based on 'keysite' observations
print(paste('start knn',Sys.time()))

#1a. Function inputs
##knn_samp<-function(n_samp,obs,cumul_samp_span,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi)

knn<-knn_samp(n_samp,obs[,which(site_names==keysite)],obs_idx=ixx_obs,cumul_samp_span,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi)

#1b. Function outputs
knn_samps<-knn   #fitted loess model per season and K

#saveRDS(cexp_fit,'./out/cexp-fit.rds')

print(paste('end knn',Sys.time()))
#-----------------------------------------
#2) Fit conditional expectation model
print(paste('start cexpfit',Sys.time()))
#2a. Function inputs
##fit_cexp<-function(obs,forecast,ctr,seasonal,start_date,end_date)

cexp_fit<-fit_cexp(obs_fit,hefs_fit,ctr,seasonal,fit_start,fit_end,leave_out_years)

#2b. Function outputs
loess_fit<-cexp_fit[[1]]   #fitted loess model per season and K
cexp_mat<-cexp_fit[[2]]    #nxK matrix of conditional expectation estimates
err_mat<-cexp_fit[[3]]     #exnxK matrix of errors calculated as: cexp - forecast

#saveRDS(cexp_fit,'./out/cexp-fit.rds')

print(paste('end cexpfit',Sys.time()))

#---------------------------------------------
#3) Fit logistic regression model for forecast intermittency
print(paste('start lregfit',Sys.time()))
#3a. Function inputs
##fit_cexp<-function(obs,forecast,ctr,seasonal,start_date,end_date)

lreg_fit<-fit_logreg(obs_fit,hefs_fit,seasonal,fit_start,fit_end,leave_out_years)

#3b. Function outputs
lreg_coef<-lreg_fit #fitted logistic regression model coefficients per season and K

#saveRDS(cexp_fit,'./out/cexp-fit.rds')

print(paste('end lregfit',Sys.time()))
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

#saveRDS(hsked_fit,'./out/hsked-fit.rds')

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

saveRDS(var_est,paste('./out/var-fit_',vers,'.rds',sep=''))

print(paste('end varfit',Sys.time()))
#---------------------------------------------
#6) Fit SEP distribution to VAR residuals by ensemble, season, and K
print(paste('start sepfit',Sys.time()))
#6a. Function inputs
#forc_err_SEPfit<-function(err_mat,use_sigest,seasonal,start_date,end_date,parallel)

sep_est<-forc_err_SEPfit(at_mat,use_sigest,distn,seasonal,fit_start,fit_end,leave_out_years,parallel=T,use_mpi)

sep_param<-sep_est         #matrix of SEP parameter estimates by ensemble, season, and K

saveRDS(sep_est,paste('./out/sep-fit_',vers,'.rds',sep=''))

print(paste('end fit',Sys.time()))

save.image("./out/model-fit_rdata.RData")

#load("./out/model-fit_rdata.RData")
#---------------------------------------------
if(generate_syn_forecasts==TRUE){
#7) Generate synthetic forecasts
print(paste('start syngen',Sys.time()))

ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC")
#update VAR or AR model indicator for synthetic generation
var_ar<-'var'
if(var_mod=='ar'){var_ar<-'ar'}

#7a. Function inputs
#syn_gen<-function(n_samp,obs,knn_samps,obs_idx,at_mat,sep_par,cexp_fit,norm_fit,var_coefs,lag,cumul_samp_span,use_mean,seasonal,fit_start,fit_end,gen_start,gen_end,parallel)

syn_gen_out <- syn_gen(n_samp,obs=obs,knn_samps=knn_samps,obs_idx=ixx_obs,at_mat,sep_par=sep_param,distn,cexp_fit=loess_fit,norm_fit,sd_arr,lreg_coefs=lreg_coef,var_coefs,var_ar=var_ar,lag=var_lag,cumul_samp_span,use_mean=var_mean,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi)

syn_hefs_fwd <- syn_gen_out[[1]]
syn_hefs_err <- syn_gen_out[[2]]

saveRDS(syn_hefs_fwd,file=paste('./out/syn_hefs_forward_',vers,'.rds',sep=''))
saveRDS(syn_hefs_err,file=paste('./out/syn_hefs_err_',vers,'.rds',sep=''))
saveRDS(ixx_gen,file='./out/ixx_gen.rds')
saveRDS(n_samp,file='./out/n_samp.rds')

print(paste('end syngen',Sys.time()))
}

print(paste('end task',Sys.time()))

rm(list=ls());gc()
####################################END#########################################