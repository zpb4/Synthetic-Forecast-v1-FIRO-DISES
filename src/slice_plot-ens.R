#script to slice out smaller subset of output array for plotting on local machine

vers<-'bvar1-sged'

syn_hefs_fwd <- readRDS(paste('out/syn_hefs_forward_',vers,'.rds',sep=''))
#syn_hefs_fwd <- readRDS('out/syn_hefs_forward.rds')
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,]

saveRDS(syn_hefs_forward,paste('out/syn_hefs_forward_',vers,'_plot-ens.rds',sep=''))

rm(list=ls());gc()


#########################################END##########################