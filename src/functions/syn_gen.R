
library(doParallel)
library(abind)

#function to rearrange forward looking forecast to obs-synched forecast
rearrange_to_fwd_forecast<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,1:(dim(forecast)[2]-i),i]<-forecast[,(i+1):dim(forecast)[2],i]
  }
  return(forecast_out)
}

syn_gen<-function(n_samp,obs,obs_idx,at_mat,sep_par,distn,cexp_fit,norm_fit,var_coefs,var_ar,lag,cumul_samp_span,use_mean,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi){

source('./src/functions/hybrid-loess-fit_fun.R')
source('./src/functions/sep_function.R')
  
#parallelization code
if(use_mpi==FALSE){
  parallel::detectCores()
  n.cores <- parallel::detectCores()
  my.cluster<-try(getDefaultCluster())
  if(class(my.cluster)=='try-error'){
    my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')}
  print(my.cluster)
  doParallel::registerDoParallel(cl = my.cluster)
  foreach::getDoParRegistered()}
if(use_mpi==TRUE){
  cl <- startMPIcluster()
  registerDoMPI(cl)}

#model specifications
ar <- lag #lag order
ens_num <- dim(at_mat)[1]
leads <- dim(at_mat)[3]
no_sites <- n_sites #NHG and MSG
n<- n_samp #number of simulations
res_env<-10 #envelope for generated residuals (must be less than X times absolute value of empirical)
obs_env<-3 #envelope for residuals during simulation (ie after VAR application), ensure simulations less than X times max obs
span<-cumul_samp_span

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

#observations for fit and generation periods
obs_fit<-obs[which(obs_idx==fit_start):which(obs_idx==fit_end)]
obs_gen<-obs[which(obs_idx==gen_start):which(obs_idx==gen_end)]

#observation arrays that match forecast dimension nxK
obs_mat<-matrix(rep(obs,dim(at_mat)[3]),ncol=dim(at_mat)[3])
obs_mat_fit<-obs_mat[which(obs_idx==fit_start):which(obs_idx==fit_end),]
obs_mat_gen<-obs_mat[which(obs_idx==gen_start):which(obs_idx==gen_end),]

#create cumulative sampling vector 
cumul_samp_fun<-function(x,span){
  out<-c()
  for(i in 1:length(x)){
    idx<-i:min(length(x),(i+span))
    out[i]<-sum(x[idx])
  }
  return(out)
}

obs_fit_samp<-cumul_samp_fun(obs_fit,span=span)
obs_gen_samp<-cumul_samp_fun(obs_gen,span=span)

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

#1c) kNN and kernel weighted sampling specification
if(seasonal=='monthly'){
  knn<-round(sqrt(length(which(ixx_gen$mon==1))))
}
tot<-sum(rep(1,knn) / 1:knn)
wts<-rep(1,knn) / 1:knn / rep(tot,knn)

#2) Synthetic forecast generations


#parallel generation across samples
synflow_out <- foreach(m = 1:n,.combine='c',.inorder=F,.packages=c('fGarch'),.export=c('sep_samp','sep_quant','lin_mod')) %dopar% {

  #arrays to store each run
  syn_flow<-array(NA,c(dim(at_mat)[1],length(ixx_gen),dim(at_mat)[3]))
  syn_err<-array(NA,c(dim(syn_flow)))
  
  #kNN sampling procedure (simulated from fitted period)
  #kNN index common across ensemble members
  knn_lst<-vector('list',12)
  
  for(i in 1:length(season_list)){
    knn_vec<-c()
    seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
    seas_gen<-which(ixx_gen$mon%in%(season_list[[i]]-1))
    gen_obs<-obs_gen_samp[seas_gen]
    fit_obs<-obs_fit_samp[seas_fit]
    for(j in 1:length(gen_obs)){
        knn_samp<-abs(fit_obs-gen_obs[j])
        #ob_smp<-matrix(rep(ob_smp,length(seas_fit)),ncol=no_sites,byrow=T)
        #y<-sqrt(apply((ob_smp - ob)^2,1,sum)) #find NEP closest by Euclidean distance
        knn_sort<-sort(knn_samp,index.return=T)
        id<-sample(knn_sort$ix[1:knn],1,prob=wts)
        knn_vec[j]<-id
    }
  knn_lst[[i]]<-knn_vec
  }
  
  
  for(e in 1:dim(at_mat)[1]){
    
    #generate empirical copula via Schaake shuffle 
    syn_empcop_knn<-vector('list',length(season_list))

    for(i in 1:length(season_list)){
      seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
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
      syn_empcop_knn_samp<-array(NA,dim(syn_empcop_mat))
      syn_empcop_knn_samp<-syn_empcop_mat[id,]
      syn_empcop_knn[[i]]<-syn_empcop_knn_samp
    }

    #specify starting residuals for primary generations
    app_mat<-syn_empcop_knn[[mo_seq[1]]][1:ar,]
    
    #generation proceeds sequentially along months to maintain VAR continuity
    for(i in 1:length(yr_seq)){
      #load ats and VAR coefficients
      syn_at_mat<-syn_empcop_knn[[mo_seq[i]]][which(yr_idx_lst[[mo_seq[i]]]==yr_seq[i]),]
      var_coef<-var_coefs[e,mo_seq[i],,]
      seas_gen<-which(ixx_gen$mon==(mo_seq[i]-1) & ixx_gen$year==yr_seq[i])
      cexp<-cexp_gen[seas_gen,]
      if(is.null(dim(syn_at_mat))==TRUE){
        syn_at_mat<-t(matrix(syn_at_mat))
        cexp<-t(matrix(cexp))
      }
      actual_obs<-obs_mat_gen[seas_gen,]
      
      #establish maximum (historical) observation for generated envelope calculations
      seas_mth_sset<-which(obs_idx$mon==(mo_seq[i]-1))
      actual_obs_seas_sset<-obs_mat[seas_mth_sset,]
      max_vec<-apply(actual_obs_seas_sset,2,max)#/apply(aob_seas,2,mean)
      
      #arrays to store simulated data
      #matrix for VAR and 'denormalized' residuals
      syn_resid_mat<-matrix(0,ncol=dim(syn_at_mat)[2],nrow=(dim(syn_at_mat)[1]+ar))
      #matrix to store VAR residuals separately (VAR is fitted to normalized residuals)
      syn_var_mat<-matrix(0,ncol=dim(syn_at_mat)[2],nrow=(dim(syn_at_mat)[1]+ar))
      #maintain previously generated (out to 'ar' lag) residuals from preceding month
      syn_var_mat[1:ar,]<-app_mat
      
      #monthly VAR residual and denormalized residual (ie synthetic forecast errors) generation
      for(j in (ar+1):(dim(syn_at_mat)[1]+ar)){
        #simulate rows of new normalized residuals from VAR coefficients
        if(var_ar=='var'){
          if(use_mean==FALSE){
            raw_var_res<-c(t(syn_var_mat[(j-1):(j-ar),])) %*% t(var_coef) + syn_at_mat[(j-ar),]}
          if(use_mean==TRUE){
            raw_var_res<-c(t(syn_var_mat[(j-1):(j-ar),]),1) %*% t(var_coef) + syn_at_mat[(j-ar),]}
          var_res<-raw_var_res[1,]
        }
        #run this if AR(p) model
        if(var_ar=='ar'){
          if(use_mean==FALSE){
            raw_var_res<-apply(syn_var_mat[(j-1):(j-ar),] * t(var_coef),2,sum) + syn_at_mat[(j-ar),]}
          if(use_mean==TRUE){
            raw_var_res<-apply(cbind(syn_var_mat[(j-1):(j-ar),],rep(1,dim(syn_at_mat)[2])) * t(var_coef),2,sum) + syn_at_mat[(j-ar),]}
          var_res<-raw_var_res
        }
        
        #corresponding row of conditional expectation matrix
        zero_sim_ref<-cexp[(j-ar),]
        
        #replace extreme VAR residuals if would create exceedingly large simulation values (> specified X times max obs)
        extr_idx<-which(((zero_sim_ref-var_res)-obs_env*max_vec)>0)
        var_res[extr_idx]<-rnorm(length(extr_idx))
        
        #zero-truncate VAR residuals to prevent auto-correlation in negative simulation space
        sub_zero_idx<-which((zero_sim_ref-var_res)<0)
        var_res[sub_zero_idx]<-zero_sim_ref[sub_zero_idx]
        syn_var_mat[j,]<-var_res
        
        #denormalize (scale) residuals
        hsked_var<-cexp[(j-ar),]
        norm_val<-c()
        #normalization value calculated per time step, per site
        for(k in 1:dim(syn_at_mat)[2]){
          norm_val[k]<-lin_mod(norm_fit[[mo_seq[i]]][[k]][1],norm_fit[[mo_seq[i]]][[k]][2],hsked_var[k])
          if(norm_val[k]<0){norm_val[k]<-sd_arr[mo_seq[i],k]}
        }
        raw_err<-var_res*norm_val
        
        #zero truncate denormalized residuals (syn forecast errors)
        sub_zero_idx<-which((zero_sim_ref-raw_err)<0)
        raw_err[sub_zero_idx]<-zero_sim_ref[sub_zero_idx]
        syn_resid_mat[j,]<-raw_err
      }
      
    syn_err[e,seas_gen,]<-syn_resid_mat[(ar+1):j,]
    syn_flow[e,seas_gen,]<-cexp_gen[seas_gen,] - syn_resid_mat[(ar+1):j,]
    app_mat<-syn_var_mat[(j-ar+1):j,]
    }
  }
return(list(syn_flow,syn_err))
  
}

#store concatenated output to single RDS array
syn_flow_comb<-array(NA,c(n,dim(at_mat)[1],length(ixx_gen),dim(at_mat)[3]))
syn_err_comb<-array(NA,c(n,dim(at_mat)[1],length(ixx_gen),dim(at_mat)[3]))

for(m in 1:n){
  #rearrange flow forecast to forward looking format
  syn_flow_comb[m,,,]<-rearrange_to_fwd_forecast(synflow_out[[m*2-1]])
  syn_err_comb[m,,,]<-synflow_out[[m*2]]
}

return(list(syn_flow_comb,syn_err_comb))

if(use_mpi==FALSE){
  stopCluster(cl = my.cluster)}
if(use_mpi==TRUE){
  closeCluster(cl)
  Rmpi::mpi.quit()}
}



###################################END######################################

