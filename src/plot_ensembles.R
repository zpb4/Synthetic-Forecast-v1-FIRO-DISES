

rm(list=ls())
library(fields)

#----------------------------------------
vers<-'bvar3-sep'

load("out/data_prep_rdata.RData")
syn_hefs_forward <- readRDS(paste('out/syn_hefs_forward_',vers,'.rds',sep=''))
ixx_sim <- readRDS('out/ixx_sim.rds') 
n_samp <- readRDS('out/n_samp.rds') 

#################################################################


##############Ensemble Density Plots###############

site<-'NHGC1'
cur_site <- which(col_names==site)


obs_rank <- 4 #pick which obs event to plot (1 largest, 2 second largest, etc)
#get date for plot
obs_date_loc <- which(obs[,cur_site]==sort(obs[,cur_site],decreasing=TRUE)[obs_rank])  #index for maximum observation
obs_date <- ixx_obs[obs_date_loc]
hefs_date_loc <- which(ixx_hefs==obs_date)
syn_hefs_date_loc <- which(ixx_sim==obs_date)

num_plot_rows <- min(n_samp+1,4)
par(mfcol=c(num_plot_rows,4),mar=c(3,3,1,1),mgp=c(3,.4,0))
#ymax <- .5*max(obs[,cur_site],hefs_forward[cur_site,,,],syn_hefs_forward[,cur_site,,,])

for (ld in c(1,4,7,10)) {
  obs_idx <- obs_date_loc-ld #ref index for ensemble plots
  hefs_idx <- hefs_date_loc-ld
  syn_idx <- syn_hefs_date_loc-ld
  
  ymax <- 1.25*max(obs[obs_idx:(obs_idx+15),cur_site])
  #1) Forward looking forecast plots
  #hefs
  plot(0:15,obs[obs_idx:(obs_idx+15),cur_site],type='l',lwd=3,main="",ylim=c(0,ymax),
       xlab='',ylab='')
  abline(v=ld,col='darkorange3',lty=2)
  text(ld+4.5,ymax-1,obs_date,col='darkorange3')
  if (length(hefs_date_loc)>0) {
    cur_sum <- round(sum(apply(hefs_forward[cur_site,,hefs_idx,],2,FUN=mean)))
    tt <- paste('HEFS ld',ld,'sum',cur_sum)
    title(tt)
    for(i in 1:n_ens){
      lines(0:15,c(obs[obs_idx,cur_site],hefs_forward[cur_site,i,hefs_idx,]),col=i)
    }
    lines(0:15,obs[obs_idx:(obs_idx+15),cur_site],lwd=3)
    mtext('Flow (kcfs)',side=2,line=1.8)
    mtext('Days from fore. init.',side=1,line=1.8)
  }
    
  #synthetic
  for (samp in 1:(num_plot_rows-1)) {
    plot(0:15,obs[obs_idx:(obs_idx+15),cur_site],type='l',lwd=3,main="",ylim=c(0,ymax),
         xlab='',ylab='')
    abline(v=ld,col='darkorange3',lty=2)
    text(ld+4.5,ymax-1,obs_date,col='darkorange3')
    if (length(syn_hefs_date_loc)>0) {
      cur_sum <- round(sum(apply(syn_hefs_forward[samp,cur_site,,syn_idx,],2,FUN=mean)))
      tt <- paste('syn',samp,'ld',ld,'sum',cur_sum)
      title(tt)
      for(i in 1:n_ens){
        lines(0:15,c(obs[obs_idx,cur_site],syn_hefs_forward[samp,cur_site,i,syn_idx,]),col=i)
      }
      lines(0:15,obs[obs_idx:(obs_idx+15),cur_site],lwd=3)
      mtext('Flow (kcfs)',side=2,line=1.8)
      mtext('Days from fore. init.',side=1,line=1.8)
    }
  }
}


########################################################



#################Aggregate validation#######################################
#cur_site <- 1
cur_samp <- 1
cur_e <- 1
#check autocorrelation for specific lead times (for one sample, ensemble member and site)
par(mfcol=c(2,4),mar=c(2,2,1,1),oma=c(3,3,1,1))
for (ld in c(1,3,5,7)) {
  acf(hefs_forward[cur_site,cur_e,,ld],xlim=c(1,10))
  title(paste("HEFS, lead",ld))
  if (ld==1) {mtext('HEFS',side=2,line=2)}
  acf(syn_hefs_forward[cur_samp,cur_site,cur_e,,ld],xlim=c(1,10),main='syn-HEFS')
  title(paste("syn-HEFS, lead",ld))
  if (ld==1) {mtext('syn-HEFS',side=2,line=2)}
}

#check cross correlation across lead times
par(mfcol=c(1,2),mar=c(2,2,1,1),oma=c(3,3,1,1))
image.plot(1:leads,1:leads,cor(hefs_forward[cur_site,cur_e,,]),zlim=c(0,1))
image.plot(1:leads,1:leads,cor(syn_hefs_forward[cur_samp,cur_site,cur_e,,]),zlim=c(0,1))

#check correlation across ensemble members
par(mfrow=c(4,2),mar=c(2,2,1,1),oma=c(3,3,1,1),mgp=c(3,.4,0))
for (ld in c(1,4,7,10)) {
  image.plot(1:n_ens,1:n_ens,cor(t(hefs_forward[cur_site,,,ld])),zlim=c(0,1))
  if (ld==1) {mtext('HEFS',side=3,line=0)}
  mtext(paste('Lead',ld),side=2,line=1.5)
  image.plot(1:n_ens,1:n_ens,cor(t(syn_hefs_forward[cur_samp,cur_site,,,ld])),zlim=c(0,1))
  if (ld==1) {mtext('syn-HEFS',side=3,line=0)}
}

########################################################



##################Event based statistics################################

if (n_samp<2) {
  {stop("cant make these plots without 2 samples or more")}
} else {
  #pick which obs event to plot (1 largest, 2 second largest, etc)
  #cur_site <- 1
  all_leads <- c(1,4,7,10)
  num_events <- 10
  syn_HEFS_cumul_mean <- array(NA,c(n_samp,num_events,length(all_leads)))
  HEFS_cumul_mean <- array(NA,c(num_events,length(all_leads)))
  obs_cumul_mean <- array(NA,c(num_events,length(all_leads)))
  syn_HEFS_cumul_sd <- array(NA,c(n_samp,num_events,length(all_leads)))
  HEFS_cumul_sd <- array(NA,c(num_events,length(all_leads)))
  
  for (i in 1:num_events) {   
    #get date for plot
    keep <- ixx_obs%in%ixx_hefs
    obs_date_loc <- which(obs[keep,cur_site]==sort(obs[keep,cur_site],decreasing=TRUE)[i])  #index for maximum observation
    obs_date <- ixx_obs[keep][obs_date_loc]
    obs_forward_date_loc <- which(ixx_obs_forward==obs_date)
    hefs_date_loc <- which(ixx_hefs==obs_date)
    syn_hefs_date_loc <- which(ixx_sim==obs_date)
  
    for (ld_indx in 1:length(all_leads)) {
      ld <- all_leads[ld_indx]
      obs_idx <- obs_date_loc-ld #ref index for ensemble plots
      obs_forward_idx <- obs_forward_date_loc-ld
      hefs_idx <- hefs_date_loc-ld
      syn_idx <- syn_hefs_date_loc-ld
    
      #ensemble mean
      HEFS_cumul_mean[i,ld_indx] <- mean(hefs_forward[cur_site,,hefs_idx,])
      syn_HEFS_cumul_mean[,i,ld_indx] <- apply(syn_hefs_forward[,cur_site,,syn_idx,],1,FUN=mean)
      obs_cumul_mean[i,ld_indx] <- mean(obs_forward_all_leads[cur_site,obs_forward_idx,])
      #ensemble standard deviation
      HEFS_cumul_sd[i,ld_indx] <- sd(apply(hefs_forward[cur_site,,hefs_idx,],FUN=mean,1))
      syn_HEFS_cumul_sd[,i,ld_indx] <- apply(apply(syn_hefs_forward[,cur_site,,syn_idx,],c(1,2),FUN=mean),1,FUN=sd)
  
    }
  }
  
  par(mfcol=c(2,4),mar=c(2,2,1,1))
  for (ld_indx in 1:length(all_leads)) {
    ymin <- min(syn_HEFS_cumul_mean[,,ld_indx],HEFS_cumul_mean[,ld_indx],obs_cumul_mean[,ld_indx])
    ymax <- max(syn_HEFS_cumul_mean[,,ld_indx],HEFS_cumul_mean[,ld_indx],obs_cumul_mean[,ld_indx])
    boxplot(syn_HEFS_cumul_mean[,,ld_indx],main=paste("Ens. Mean, Lead",all_leads[ld_indx]),ylim=c(ymin,ymax))
    points(HEFS_cumul_mean[,ld_indx],col="red",pch=16)
    points(obs_cumul_mean[,ld_indx],col="blue",pch=17)
    ymin <- 0.2 #min(syn_HEFS_cumul_sd[,,ld_indx],HEFS_cumul_sd[,ld_indx])
    ymax <- 1.0 #max(syn_HEFS_cumul_sd[,,ld_indx],HEFS_cumul_sd[,ld_indx])
    boxplot(syn_HEFS_cumul_sd[,,ld_indx],main=paste("Ens. SD, Lead",all_leads[ld_indx]),ylim=c(ymin,ymax))
    points(HEFS_cumul_sd[,ld_indx],col="red",pch=16)
  }
}
########################################################################
