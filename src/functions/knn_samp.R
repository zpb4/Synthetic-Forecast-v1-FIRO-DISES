


knn_samp<-function(n_samp,obs,obs_idx,cumul_samp_span,seasonal,fit_start,fit_end,gen_start,gen_end,use_mpi){
  
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
  n <- n_samp #number of simulations
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
  
  #observations for fit and generation periods
  obs_fit<-obs[which(obs_idx==fit_start):which(obs_idx==fit_end)]
  obs_gen<-obs[which(obs_idx==gen_start):which(obs_idx==gen_end)]
  
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
  
  #1c) kNN and kernel weighted sampling specification
  knn<-round(sqrt(length(ixx_fit))) #annual knn based on sqrt(n) where n is length of fit period
  
  if(seasonal=='monthly'){ #monthly knn based on sqrt(n) where n is length of a monthly subset of fit period
    knn<-round(sqrt(length(which(ixx_fit$mon==1))))
  }
  if(seasonal!='monthly'&seasonal!='annual'){ #multimonth knn based on sqrt(n) where n is length of a multi-monthly subset of fit period
    knn<-round(sqrt(length(which(ixx_fit$mon%in%(season_list[[i]]-1)))))
  }
  
  tot<-sum(rep(1,knn) / 1:knn)
  wts<-rep(1,knn) / 1:knn / rep(tot,knn)
  
  #2) Synthetic forecast generations
  
  #kNN sampling procedure (simulated from fitted period)
  #kNN index common across ensemble members
  knn_out<-foreach(m = 1:n,.inorder=F)%dopar%{
    knn_lst<-vector('list',12)

    for(i in 1:length(season_list)){
      knn_vec<-c()
      seas_fit<-which(ixx_fit$mon%in%(season_list[[i]]-1))
      seas_gen<-which(ixx_gen$mon%in%(season_list[[i]]-1))
      gen_obs<-obs_gen_samp[seas_gen]
      fit_obs<-obs_fit_samp[seas_fit]
      for(j in 1:length(gen_obs)){
        knn_samp<-abs(fit_obs-gen_obs[j])
        knn_sort<-sort(knn_samp,index.return=T)
        id<-sample(knn_sort$ix[1:knn],1,prob=wts)
        knn_vec[j]<-id
      }
      knn_lst[[i]]<-knn_vec
    }
    return(knn_lst)
  }

  return(knn_out) #return an nsamps long list of kNN sampled knn_lsts
}
