

#function to rearrange forward looking forecast to obs-synched forecast
rearrange_to_fwd_forecast<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,1:(dim(forecast)[2]-i),i]<-forecast[,(i+1):dim(forecast)[2],i]
  }
  return(forecast_out)
}

logreg_pred<-function(x,coeffs){
  p <- exp(coeffs[1] +  x * coeffs[2]) / (1 + exp(coeffs[1] +  x * coeffs[2]))
  prd<-c()
  for(i in 1:length(x)){
    prd[i]<-rbinom(1,1,p)
  }
  return(prd)
}

syn_gen<-function(n_samp,obs,knn_samps,obs_idx,at_mat,sep_par,distn,cexp_fit,norm_fit,sd_arr,lreg_coefs,var_coefs,var_ar,lag,cumul_samp_span,use_mean,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi){

source('./src/functions/hybrid-loess-fit_fun.R')
source('./src/functions/sep_function.R')
  
#parallelization code
if(use_mpi==FALSE){
  library(doParallel)
  parallel::detectCores()
  n.cores <- parallel::detectCores()
  my.cluster<-try(getDefaultCluster())
  if(class(my.cluster)=='try-error'|class(my.cluster)=='NULL'){
    my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')}
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()}
if(use_mpi==TRUE){
  library(doMPI)
  cl <- startMPIcluster()
  registerDoMPI(cl)}

#model specifications
ar_lag <- lag #lag order
ens_num <- dim(at_mat)[1]
Kdim <- dim(at_mat)[3]
if(is.null(dim(obs)[2])==TRUE){n_sites <- 1}else{n_sites <- dim(obs)[2]}
leads <- Kdim/n_sites
n <- n_samp #number of simulations
res_env<-10 #envelope for generated residuals (must be less than 'res_env' times absolute value of empirical)
obs_env<-3 #envelope for synthetic forecasts during simulation, ensure simulations less than 'obs_env' times max obs
span<-cumul_samp_span #

#date time indices
#fitted period
ix_fit<-seq(as.Date(fit_start),as.Date(fit_end),'day')
ixx_fit<-as.POSIXlt(ix_fit)

#generation period
ix_gen<-seq(as.Date(gen_start),as.Date(gen_end),'day') 
ixx_gen<-as.POSIXlt(ix_gen)

#seasonal lists
season_list<-seasonal

if(seasonal=='monthly'){
  season_list<-vector('list',12)
  for(i in 1:12){
    season_list[[i]]<-i
  }
}

if(seasonal=='annual'){
  season_list<-list(1:12)
}

#observation arrays that match forecast dimension nxK
#process obs and forecasts when only 1 site
if(is.null(dim(obs)[2])==TRUE){
  obs_mat<-matrix(rep(obs,dim(at_mat)[3]),ncol=dim(at_mat)[3],byrow=F)
}

#concatenate across lead time dimension if more than one site
if(dim(obs)[2]>1){
  ncols=dim(at_mat)[3]/dim(obs)[2]
  obs_mat<-matrix(rep(obs[,1],ncols),ncol=ncols,byrow=F)
  for(i in 2:dim(obs)[2]){
    obs_sset<-matrix(rep(obs[,i],ncols),ncol=ncols,byrow=F)
    obs_mat<-cbind(obs_mat,obs_sset)}
  } 
  
obs_mat_fit<-obs_mat[which(obs_idx==fit_start):which(obs_idx==fit_end),]
obs_mat_gen<-obs_mat[which(obs_idx==gen_start):which(obs_idx==gen_end),]

#generate conditional expectation matrix
cexp_gen<-array(NA,dim(obs_mat_gen))

for(i in 1:length(season_list)){
  seas<-which(ixx_gen$mon%in%(season_list[[i]]-1))
  for(k in 1:dim(at_mat)[3]){
    hy_fit<-cexp_fit[[i]][[k]]
    cexp<-hyb_loess_out(hy_fit[[1]],hy_fit[[2]],hy_fit[[3]],hy_fit[[4]],obs_mat_gen[seas,k])
    cexp[cexp<0]<-0
    bad_vals<-which(is.na(cexp)==T|cexp==Inf|cexp==-Inf)
    cexp[bad_vals]<-obs_mat_gen[seas,k][bad_vals]
    cexp_gen[seas,k]<-cexp
  }
}

#1b) create sequential date/time index for sequential VAR simulation (continuity between months)
if(seasonal=='monthly'){
  yr_idx<-ixx_gen$year
  yr_idx_lst<-vector('list',12)

  for(i in 1:12){
    seas<-which(ixx_gen$mon==(i-1))
    yr_idx_lst[[i]]<-yr_idx[seas]
  }
  
  yrs_uniq<-unique(ixx_gen$year)
  
  yr_seq<-c()
  mo_seq<-c()
  for(i in 1:length(yrs_uniq)){
    yr_seq<-c(yr_seq,rep(yrs_uniq[i],length(unique(ixx_gen$mon[ixx_gen$year==yrs_uniq[i]]))))
    mo_seq<-c(mo_seq,unique(ixx_gen$mon[ixx_gen$year==yrs_uniq[i]]))
  }
  mo_seq<-mo_seq+1 #set to 1 indexed vector
}

#2) Synthetic forecast generations

#parallel generation across samples
synflow_out <- foreach(m = 1:n,.combine='c',.inorder=F,.packages=c('fGarch'),.export=c('sep_samp','sep_quant','lin_mod','logreg_pred')) %dopar% {

  #arrays to store each run
  syn_flow<-array(NA,c(dim(at_mat)[1],length(ixx_gen),dim(at_mat)[3]))
  syn_err<-array(NA,c(dim(syn_flow)))
  
  knn_lst<-knn_samps[[m]]
  
  for(e in 1:dim(at_mat)[1]){
    
    #generate empirical copula via Schaake shuffle 
    syn_empcop_knn<-vector('list',length(season_list))

    for(i in 1:length(season_list)){
      seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
      seas_gen<-which(ixx_gen$mon%in%(season_list[[i]]-1))
      inp_mat<-at_mat[e,seas_fit,]
      empcop<-apply(inp_mat,2,function(x){rank(x,ties.method = 'random')})
      syn_empcop_mat<-array(NA,c(dim(inp_mat)[1],dim(inp_mat)[2]))
      for(k in 1:dim(inp_mat)[2]){
        #generate new ats
        if(distn=='sep'){
          syn_at<-sep_samp(dim(inp_mat)[1],theta=sep_par[e,i,k,1],sigma=sep_par[e,i,k,2],beta=sep_par[e,i,k,3],alpha=sep_par[e,i,k,4])}
        if(distn=='sged'){
          syn_at<-rsged(dim(inp_mat)[1],mean=sep_par[e,i,k,1],sd=sep_par[e,i,k,2],nu=sep_par[e,i,k,3],xi=sep_par[e,i,k,4])}
        #ensure absolute value of ats within specified envelop
        extr_syn_at_idx<-which(abs(syn_at)>(res_env*max(abs(inp_mat[,k]))))  
        syn_at[extr_syn_at_idx]<-rnorm(length(extr_syn_at_idx)) 
        #reorder new ats to match empirical data
        rnk_syn_at<-rank(syn_at,ties.method = 'random')
        for(j in 1:length(rnk_syn_at)){
        syn_empcop_mat[j,k]<-syn_at[which(rnk_syn_at==empcop[j,k])]
        }
      }
      id<-knn_lst[[i]]
      syn_empcop_knn_samp<-array(NA,c(length(ixx_gen),dim(inp_mat)[2]))
      syn_empcop_knn_samp<-syn_empcop_mat[id,]
      syn_empcop_knn[[i]]<-syn_empcop_knn_samp
    }

    #specify starting residuals for primary generations
    app_mat<-syn_empcop_knn[[mo_seq[1]]][1:ar_lag,]
    
    #generation proceeds sequentially along months to maintain VAR continuity
    for(i in 1:length(yr_seq)){
      #load ats and VAR coefficients
      syn_at_mat<-syn_empcop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
      var_coef<-var_coefs[e,mo_seq[i],,]
      seas_gen<-which(ixx_gen$mon==(mo_seq[i]-1) & ixx_gen$year==yr_seq[i])
      cexp<-cexp_gen[seas_gen,]
      actual_obs<-obs_mat_gen[seas_gen,]
      if(is.null(dim(syn_at_mat))==TRUE){
        syn_at_mat<-t(matrix(syn_at_mat))
        cexp<-t(matrix(cexp))
        actual_obs<-t(matrix(actual_obs))
      }
      
      #establish maximum (historical) observation for generated envelope calculations
      seas_mth_sset<-which(obs_idx$mon==(mo_seq[i]-1))
      actual_obs_seas_sset<-obs_mat[seas_mth_sset,]
      max_vec<-apply(actual_obs_seas_sset,2,max)#/apply(aob_seas,2,mean)
      
      #arrays to store simulated data
      #matrix for VAR and 'denormalized' residuals
      syn_resid_mat<-matrix(0,ncol=dim(syn_at_mat)[2],nrow=(dim(syn_at_mat)[1]+ar_lag))
      #matrix to store VAR residuals separately (VAR is fitted to normalized residuals)
      syn_var_mat<-matrix(0,ncol=dim(syn_at_mat)[2],nrow=(dim(syn_at_mat)[1]+ar_lag))
      #maintain previously generated (out to 'ar' lag) residuals from preceding month
      syn_var_mat[1:ar_lag,]<-app_mat
      
      #monthly VAR residual and denormalized residual (ie synthetic forecast errors) generation
      for(j in (ar_lag+1):(dim(syn_at_mat)[1]+ar_lag)){
        #simulate rows of new normalized residuals from VAR coefficients
        if(var_ar=='var'){
          if(use_mean==FALSE){
            raw_var_res<-c(t(syn_var_mat[(j-1):(j-ar_lag),])) %*% t(var_coef) + syn_at_mat[(j-ar_lag),]}
          if(use_mean==TRUE){
            raw_var_res<-c(t(syn_var_mat[(j-1):(j-ar_lag),]),1) %*% t(var_coef) + syn_at_mat[(j-ar_lag),]}
          var_res<-raw_var_res[1,]
        }
        #run this if AR(p) model
        if(var_ar=='ar'){
          if(use_mean==FALSE){
            raw_var_res<-apply(syn_var_mat[(j-1):(j-ar_lag),] * t(var_coef),2,sum) + syn_at_mat[(j-ar_lag),]}
          if(use_mean==TRUE){
            raw_var_res<-apply(cbind(syn_var_mat[(j-1):(j-ar_lag),],rep(1,dim(syn_at_mat)[2])) * t(var_coef),2,sum) + syn_at_mat[(j-ar_lag),]}
          var_res<-raw_var_res
        }
        
        #corresponding row of conditional expectation matrix
        zero_sim_ref<-cexp[(j-ar_lag),]
        act_obs<-actual_obs[(j-ar_lag),]
        
        #replace extreme VAR residuals if would create exceedingly large simulation values (> specified X times max obs)
        extr_idx<-which(((zero_sim_ref-var_res)-obs_env*max_vec)>0)
        var_res[extr_idx]<-rnorm(length(extr_idx))
        
        #zero-truncate VAR residuals to prevent auto-correlation in negative simulation space
        ##sub_zero_idx<-which((zero_sim_ref-var_res)<0)
        ##var_res[sub_zero_idx]<-zero_sim_ref[sub_zero_idx]
        ##syn_var_mat[j,]<-var_res
        
        #denormalize (scale) residuals
        hsked_var<-cexp[(j-ar_lag),]
        norm_val<-c()
        #normalization value calculated per time step, per site
        for(k in 1:dim(syn_at_mat)[2]){
          norm_val[k]<-lin_mod(norm_fit[[mo_seq[i]]][[k]][1],norm_fit[[mo_seq[i]]][[k]][2],hsked_var[k])
          if(norm_val[k]<0){norm_val[k]<-sd_arr[mo_seq[i],k]}
        }
        raw_err<-var_res*norm_val
        
        #replace extreme raw errors if would create exceedingly large simulation values (> specified X times max obs)
        extr_idx<-which(((zero_sim_ref-raw_err)-obs_env*max_vec)>0)
        raw_err[extr_idx]<-rnorm(length(extr_idx))*mean(sd_arr[mo_seq[i],])
        
        #apply intermittency model
        for(k in 1:dim(syn_at_mat)[2]){
          if((zero_sim_ref[k]-raw_err[k])<=0){
            if(distn=='sep'){
              samp_vec<-sep_samp(100,theta=sep_par[e,mo_seq[i],k,1],sigma=sep_par[e,mo_seq[i],k,2],beta=sep_par[e,mo_seq[i],k,3],alpha=sep_par[e,mo_seq[i],k,4])*norm_val[k]}
            if(distn=='sged'){
              samp_vec<-rsged(100,mean=sep_par[e,mo_seq[i],k,1],sd=sep_par[e,mo_seq[i],k,2],nu=sep_par[e,mo_seq[i],k,3],xi=sep_par[e,mo_seq[i],k,4])*norm_val[k]}
            int_pred<-try(logreg_pred(act_obs[k],lreg_coefs[[mo_seq[i]]][[k]]),silent=TRUE)
            if(class(int_pred)=='try-error'|class(int_pred)=='integer'){int_pred<-0}
            if(int_pred==0){raw_err[k]<-zero_sim_ref[k]}
            if(int_pred==1){
              new_res<-try(sample(samp_vec[samp_vec<zero_sim_ref[k]],1),silent=TRUE)
              if(class(new_res)=='try-error'){new_res<-runif(1)*zero_sim_ref[k]}
              raw_err[k]<-new_res}
            var_res[k]<-raw_err[k]/norm_val[k]
            if(is.na(var_res[k])==T|var_res[k]==Inf|var_res[k]==-Inf){var_res[k]<-raw_err[k]}
          }
        }
        syn_resid_mat[j,]<-raw_err
        syn_var_mat[j,]<-var_res
      }
      
    syn_err[e,seas_gen,]<-syn_resid_mat[(ar_lag+1):j,]
    syn_flow[e,seas_gen,]<-cexp_gen[seas_gen,] - syn_resid_mat[(ar_lag+1):j,]
    app_mat<-syn_var_mat[(j-ar_lag+1):j,]
    }
  }
return(list(syn_flow,syn_err))
  
}

#store concatenated output to single RDS array
syn_flow_comb<-array(NA,c(n,n_sites,dim(at_mat)[1],length(ixx_gen),leads))
syn_err_comb<-array(NA,c(n,n_sites,dim(at_mat)[1],length(ixx_gen),leads))

for(m in 1:n){
  #rearrange flow forecast to forward looking format
  for(s in 1:n_sites){
    col_idx<-(leads*(s-1)+1):(leads*s)
    syn_flow_comb[m,s,,,]<-rearrange_to_fwd_forecast(synflow_out[[m*2-1]][,,col_idx])
    syn_err_comb[m,s,,,]<-synflow_out[[m*2]][,,col_idx]
  }
}

return(list(syn_flow_comb,syn_err_comb))

if(use_mpi==FALSE){
  stopCluster(cl = my.cluster)}
if(use_mpi==TRUE){
  closeCluster(cl)
  Rmpi::mpi.quit()}
}



###################################END######################################

