#Forecast Verification Functions
fwd_forecast_rearrange<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,(i+1):dim(forecast)[2],i]<-forecast[,1:(dim(forecast)[2]-i),i]
  }
  return(forecast_out)
}

ens_rank<-function(ens_forc,obs){
  f_sort<-sort(ens_forc)
  dif<-f_sort - obs
  if(all(dif<0)==T){rnk<-length(dif)+1} else {y<-min(dif[dif>=0]);rnk<-which(dif==y)}
  if(length(rnk)>1){rnk<-sample(rnk,1)}
  return(rnk)
}

roll_sum<-function(x){out<-c();out[1]<-x[1];
for(i in 1:(length(x)-1)){out[i+1]<-out[i]+x[i+1]};return(out)}

chi2<-function(freqs){
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  chi2<-((m)/n)*sum((freqs - (n / m))^2)
  prob<-pgamma(chi2,shape=((m-1)/2),scale=2,lower.tail=F)
  return(list(chi2,prob))
}

chi2_comp<-function(freq_ref,freq_tst){
  freq_ref[freq_ref<1]<-1
  m<-length(freq_ref) 
  n<-sum(freq_ref)
  chi2<-sum(((freq_tst - freq_ref)^2)/freq_ref)
  prob<-pgamma(chi2,shape=((m-1)/2),scale=2,lower.tail=F)
  return(list(chi2,prob))
}

RI<-function(freqs){
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  ri<-(1/n)*sum(abs(freqs - (n / (m+1))))
  return(ri)
}

Ent1<-function(freqs){
  freqs<-freqs[freqs>0]
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  ent<-((-1)/log(m+1))*sum((freqs/n)*log(freqs/n))
  return(ent)
}

Ent<-function(freqs){
  freqs[freqs<1]<-1
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  ent<-((-1)/log(m+1))*sum((freqs/n)*log(freqs/n))
  return(ent)
}

mn_ens<-function(ens){
  x_bar_t<-(1/length(ens))*sum(ens)
  s2<-(1/length(ens))*sum((ens - x_bar_t)^2)
  return(x_bar_t)
}

st2_ens<-function(ens){
  x_bar_t<-(1/length(ens))*sum(ens)
  s2<-(1/length(ens))*sum((ens - x_bar_t)^2)
  return(s2)
}

bin_spread<-function(ens_vars,ens_num){
  sprd<-(ens_num/(ens_num-1))*(1/length(ens_vars))*sum(ens_vars)
  return(sqrt(sprd))
}

bin_mse<-function(ens_mns,obs,ens_num){
  mse<-(ens_num/(ens_num+1))*(1/length(ens_mns))*sum((ens_mns - obs)^2)
  return(sqrt(mse))
}

eCRPS<-function(ens,obs){
  m<-length(ens)
  t1<-(1/m)*sum(abs(ens-obs))
  t2<-array(0,c(m,m))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      t2[i,j]<-abs(ens[i]-ens[j])
    }
  }
  t2_res<-(1 / (m*(m-1)))*sum(t2)
  return(t1 - t2_res)
}

eDSS<-function(var_ens,mn_ens,obs){
  ss<-log(var_ens) + (((obs - mn_ens)^2)/var_ens^2)
  return(ss)
}

##################################END###################################################