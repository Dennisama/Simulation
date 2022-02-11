#Simulation
library(astsa)
library(stats)
library(arfima)
library(fracdiff)
require(svMisc)

setwd("H:/TimeSeriesTest/LiAndBiao/Table8")

#Size of the sample
n0 = 512

#The number of beta.k, we have beta.1, beta.2, ...
v = 30
k0 = c(1:v)

#The Fourier Frequencies
m = n0/2 - 1
j = c(1:m)
mu.j = 2*pi*j/n0

###1/3###
# mu.j.filter.3 = mu.j[seq(2,m,3)]
# sub.length.3 = m/3
# n0.sub.3 = 2*(sub.length.3 + 1)

###1/5###
# mu.j.filter.5 = mu.j[seq(3,m,5)]
# sub.length.5 = m/5
# n0.sub.5 = 2*(sub.length.5 + 1)

#The number of the simulation
N = 5000

#beta.0.ar1 is used to store all beta.0.
# beta.k.ar1 is used to store all beta.k, this is a 6-row, 5000-column matrix.
beta.0.pow = rep(0,N)
beta.k.pow = matrix(0,nrow = v,ncol = N)
beta.0.emp = rep(0,N)
beta.k.emp = matrix(0,nrow = v,ncol = N)
###1/3###
# beta.0.pow.3 = rep(0,N)
# beta.k.pow.3 = matrix(0,nrow = v,ncol = N)
###1/5###
# beta.0.pow.5 = rep(0,N)
# beta.k.pow.5 = matrix(0,nrow = v,ncol = N)
#The following two are used to calculate the empirical values.
###1/3###
# beta.0.emp.3 = rep(0,N)
# beta.k.emp.3 = matrix(0,nrow = v,ncol = N)
###1/5###
# beta.0.emp.5 = rep(0,N)
# beta.k.emp.5 = matrix(0,nrow = v,ncol = N)

#Spectral Density Function
#@parameter w.unit: this is our fourier frequency
#@parameter phi.1: this is the coefficient of AR(1), the true value is 0.8
#@parameter sig2_f: the variance of noises, notated as sigma square.
####Table 5####
sp.ar1.f.est = function(w.unit,phi.1,sig2_f){(sig2_f/((1 + phi.1^2) - 2*phi.1*cos(w.unit)))/(2*pi)}
####Table 6####
# sp.long.f.est = function(w.loc,phi.loc,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}
####Table 7####
# sp.ar2.f.est = function(w.loc,phi.loc,theta.loc,sig2_f){sig2_f/(2*pi)*(1 + theta.loc^2 + 2*theta.loc*cos(w.loc))/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}
####Table 8####
# sp.long.f.est = function(w.loc,d.loc,sig2_f){sig2_f*(4*(sin(w.loc/2))^2)^(-d.loc)/(2*pi)}


#Beta.0 Function
#It is used to calculate beta.0
#@parameter I: The Periodogram
#@parameter sp: The Spectral Density
beta.0.f = function(I,sp){
  sum(I/sp)*2*sqrt(pi)/n0
}

#Beta.k Function
#It is used to calculate beta.k
#@parameter I: The Periodogram
#@parameter sp: The Spectral Density
beta.k.f = function(k,I,sp,mu.j.in){
  sum(cos(k*mu.j.in)*I/sp)*2*sqrt(2*pi)/n0
}

#r.rate.1.thr is to store the rejection rate at 5%
#r,rate.2.thr is to store the rejection rate at 10%
r.rate.1 = rep(0,(v + 1))
r.rate.2 = rep(0,(v + 1))
r.rate.1.emp = rep(0,(v + 1))
r.rate.2.emp = rep(0,(v + 1))
###1/3###
# r.rate.1.thr.third = rep(0,(v + 1))
# r.rate.2.thr.third = rep(0,(v + 1))
###1/5###
# r.rate.1.thr.fifth = rep(0,(v + 1))
# r.rate.2.thr.fifth = rep(0,(v + 1))
#The following two are for comparing empirical crictical values.
###1/3###
# r.rate.1.emp.third = rep(0,(v + 1))
# r.rate.2.emp.third = rep(0,(v + 1))
###1/5###
# r.rate.1.emp.fifth = rep(0,(v + 1))
# r.rate.2.emp.fifth = rep(0,(v + 1))

phi.est.v.pow = rep(0,N)
# theta.est.v.pow = rep(0,N)
# d.est.v.pow = rep(0,N)
sig2.est.v.pow = rep(0,N)

phi.est.v = rep(0,N)
# theta.est.v = rep(0,N)
# d.est.v = rep(0,N)
sig2.est.v = rep(0,N)


###############
##Simulations##
###############



for(sim in c(1:N)){
  #######Table 5######
  # simulation of AR(2)
  x.sim.ar1.pow = arima.sim(list(order = c(2,0,0),ar = c(.8,-.15)),n = n0)
  prd.pow = Mod(fft(x.sim.ar1.pow- mean(x.sim.ar1.pow)))^2/(2*pi*n0)
  x.prd.ar1.pow = prd.pow[2:(n0/2)]
  # simulation of AR(1)
  x.sim.ar1 = arima.sim(list(order = c(1,0,0),ar = .8),n = n0)
  prd = Mod(fft(x.sim.ar1 - mean(x.sim.ar1)))^2/(2*pi*n0)
  x.prd.ar1 = prd[2:(n0/2)]

  ##MLE##
  x.sim.pow = arima0(x.sim.ar1.pow,order = c(1,0,0),method = "ML",include.mean = FALSE)
  phi.est.pow = x.sim.pow$coef
  sig2.est.pow = x.sim.pow$sigma2

  phi.est.v.pow[sim] = phi.est.pow
  sig2.est.v.pow[sim] = sig2.est.pow
  sp.ar1.est.pow = sp.ar1.f.est(w.unit = mu.j,phi.1 = phi.est.pow,sig2_f = 1)#sig2.est.pow

  x.sim = arima0(x.sim.ar1,order = c(1,0,0),method = "ML",include.mean = FALSE)
  phi.est = x.sim$coef
  sig2.est = x.sim$sigma2

  phi.est.v[sim] = phi.est
  sig2.est.v[sim] = sig2.est
  sp.ar1.est = sp.ar1.f.est(w.unit = mu.j,phi.1 = phi.est,sig2_f = 1)#sig2.est
  #######
  ####################
  
  #######Table 6######
  # # simulation of ARMA(1,1)
  # x.sim.arma.pow = arima.sim(list(order = c(1,0,1),ar = .8,ma = .2),n = n0)
  # prd.arma.pow = Mod(fft(x.sim.arma.pow- mean(x.sim.arma.pow)))^2/(2*pi*n0)
  # x.prd.arma.pow = prd.arma.pow[2:(n0/2)]
  # # simulation of ARFIMA(1,d,0)
  # x.sim.long = (fracdiff.sim(n0,ar = .1,d = .4))$series
  # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
  # x.prd.long = prd.long[2:(n0/2)]
  # 
  # ##MLE##
  # x.sim.pow = arfima(x.sim.arma.pow,order = c(1,0,0),quiet = TRUE)$modes
  # d.est.pow = x.sim.pow[[1]]$dfrac
  # phi.est.pow = x.sim.pow[[1]]$phi
  # sig2.est.pow = x.sim.pow[[1]]$sigma2
  # 
  # d.est.v.pow[sim] = d.est.pow
  # phi.est.v.pow[sim] = phi.est.pow
  # sig2.est.v.pow[sim] = sig2.est.pow
  # sp.arma.est = sp.long.f.est(w.loc = mu.j,phi.loc = phi.est.pow,d.loc = d.est.pow,sig2_f = 1)#sig2.est.pow
  # 
  # x.sim = arfima(x.sim.long,order = c(1,0,0),quiet = TRUE)$modes
  # d.est = x.sim[[1]]$dfrac
  # phi.est = x.sim[[1]]$phi
  # sig2.est = x.sim[[1]]$sigma2
  # 
  # d.est.v[sim] = d.est
  # phi.est.v[sim] = phi.est
  # sig2.est.v[sim] = sig2.est
  # sp.long.est = sp.long.f.est(w.loc = mu.j,phi.loc = phi.est,d.loc = d.est,sig2_f = 1)#sig2.est
  # #######
  ####################
  
  #######Table 7######
  # # simulation of ARFIMA(d)
  # x.sim.long = (fracdiff.sim(n0,d = 0.4))$series
  # prd = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
  # x.prd.long = prd[2:(n0/2)]
  # # simulation of ARMA(1,0,1)
  # x.sim.arma = arima.sim(list(order = c(1,0,1),ar = .8,ma = .2),n = n0)
  # prd.arma = Mod(fft(x.sim.arma- mean(x.sim.arma)))^2/(2*pi*n0)
  # x.prd.arma = prd.arma[2:(n0/2)]
  # 
  # ##MLE##
  # x.sim.pow = arima0(x.sim.long,order = c(1,0,1),method = "ML",include.mean = FALSE)
  # phi.est.pow = x.sim.pow$coef[1]
  # theta.est.pow = x.sim.pow$coef[2]
  # sig2.est.pow = x.sim.pow$sigma2
  # 
  # phi.est.v.pow[sim] = phi.est.pow
  # theta.est.v.pow[sim] = theta.est.pow
  # sig2.est.v.pow[sim] = sig2.est.pow
  # sp.ar2.est.pow = sp.ar2.f.est(w.loc = mu.j,phi.loc = phi.est.pow,theta.loc = theta.est.pow,sig2_f = sig2.est.pow)#
  # 
  # x.sim = arima0(x.sim.arma,order = c(1,0,1),method = "ML",include.mean = FALSE)
  # phi.est = x.sim$coef[1]
  # theta.est = x.sim$coef[2]
  # sig2.est = x.sim$sigma2
  # 
  # phi.est.v[sim] = phi.est
  # theta.est.v[sim] = theta.est
  # sig2.est.v[sim] = sig2.est
  # sp.ar2.est = sp.ar2.f.est(w.loc = mu.j,phi.loc = phi.est,theta.loc = theta.est,sig2_f = sig2.est)#
  # #######
  ####################
  
  #######Table 8######
  # # simulation of ARFIMA(1,d,0)
  # x.sim.long.pow = (fracdiff.sim(n0,ar = .1,d = .4))$series
  # prd.long.pow = Mod(fft(x.sim.long.pow- mean(x.sim.long.pow)))^2/(2*pi*n0)
  # x.prd.long.pow = prd.long.pow[2:(n0/2)]
  # # simulation of ARFIMA(d)
  # x.sim.long.bm = (fracdiff.sim(n0,d = .4))$series
  # prd.long.bm = Mod(fft(x.sim.long.bm- mean(x.sim.long.bm)))^2/(2*pi*n0)
  # x.prd.long.bm = prd.long.bm[2:(n0/2)]
  # 
  # ##MLE##
  # x.sim.pow = arfima(x.sim.long.pow,order = c(0,0,0),quiet = TRUE)$modes
  # d.est.pow = x.sim.pow[[1]]$dfrac
  # sig2.est.pow = x.sim.pow[[1]]$sigma2
  # d.est.v.pow[sim] = d.est.pow
  # sig2.est.v.pow[sim] = sig2.est.pow
  # sp.long.est.pow = sp.long.f.est(w.loc = mu.j,d.loc = d.est.pow,sig2_f = 1)#sig2.est.pow
  # 
  # x.sim.bm = arfima(x.sim.long.bm,order = c(0,0,0),quiet = TRUE)$modes
  # d.est.bm = x.sim.bm[[1]]$dfrac
  # sig2.est.bm = x.sim.bm[[1]]$sigma2
  # d.est.v[sim] = d.est.bm
  # sig2.est.v[sim] = sig2.est.bm
  # sp.long.est.bm = sp.long.f.est(w.loc = mu.j,d.loc = d.est.bm,sig2_f = 1)#sig2.est.pow
  # ########
  #####################
  
  ###Filter Periodogram###
  ###1/3###
  # x.prd.pow.filter.3 = x.prd.long.pow[seq(2,m,3)]
  # x.prd.bm.filter.3 = x.prd.long.bm[seq(2,m,3)]
  ###1/5###
  # x.prd.pow.filter.5 = x.prd.long.pow[seq(3,m,5)]
  # x.prd.bm.filter.5 = x.prd.long.bm[seq(3,m,5)]
  ###Filter Spectral Density###
  ###1/3###
  # sp.pow.filter.3 = sp.long.est.pow[seq(2,m,3)]
  # sp.bm.filter.3 = sp.long.est.bm[seq(2,m,3)]
  ###1/5###
  # sp.pow.filter.5 = sp.long.est.pow[seq(3,m,5)]
  # sp.bm.filter.5 = sp.long.est.bm[seq(3,m,5)]
  
  
  # Calculate beta.0.
  # This is for the case pars are unkown, please refer doc for details#######
  beta.0.pow[sim] = beta.0.f(I = x.prd.ar1.pow,sp = sp.ar1.est.pow)
  beta.0.emp[sim] = beta.0.f(I = x.prd.ar1,sp = sp.ar1.est)
  ###1/3###
  # beta.0.pow.3[sim] = beta.0.f(I = x.prd.pow.filter.3,sp = sp.pow.filter.3)
  # beta.0.emp.3[sim] = beta.0.f(I = x.prd.bm.filter.3,sp = sp.bm.filter.3)
  ###1/5###
  # beta.0.pow.5[sim] = beta.0.f(I = x.prd.pow.filter.5,sp = sp.pow.filter.5)
  # beta.0.emp.5[sim] = beta.0.f(I = x.prd.bm.filter.5,sp = sp.bm.filter.5)
  ###########################################################################
  
  #Calculate beta1, beta2, ...
  for (i in k0) {
    beta.k.pow[i,sim] = beta.k.f(k = i,I = x.prd.ar1.pow, sp = sp.ar1.est.pow, mu.j.in = mu.j)
    beta.k.emp[i,sim] = beta.k.f(k = i,I = x.prd.ar1, sp = sp.ar1.est, mu.j.in = mu.j)
    ###1/3###
    # beta.k.pow.3[i,sim] = beta.k.f(k = i,I = x.prd.pow.filter.3, sp = sp.pow.filter.3, mu.j.in = mu.j.filter.3)
    # beta.k.emp.3[i,sim] = beta.k.f(k = i,I = x.prd.bm.filter.3, sp = sp.bm.filter.3, mu.j.in = mu.j.filter.3)
    ###1/5###
    # beta.k.pow.5[i,sim] = beta.k.f(k = i,I = x.prd.pow.filter.5, sp = sp.pow.filter.5, mu.j.in = mu.j.filter.5)
    # beta.k.emp.5[i,sim] = beta.k.f(k = i,I = x.prd.bm.filter.5, sp = sp.bm.filter.5, mu.j.in = mu.j.filter.5)
  }
  
  progress(sim, max.value = N)
  if (sim == N) cat("Done!\n")
}
#exclude beta.0, starting at 1.
for (v0 in c(0:v)) {
  alpha.1 = .05
  alpha.2 = .1
  
  # when sig2 is unknown, beta.0 is excluded, df = v0
  chi.sq.true.1 = qchisq(alpha.1,df = v0 + 1,lower.tail = FALSE)
  chi.sq.true.2 = qchisq(alpha.2,df = v0 + 1,lower.tail = FALSE)
  
  # if (v0 == 1) {
  #   x.std.pow = n0/(2*pi)*(beta.k.pow[1,])^2
  #   x.std.emp = n0/(2*pi)*(beta.k.emp[1,])^2
  # } else {
  #   x.std.pow = n0/(2*pi)*colSums((beta.k.pow[c(1:v0),])^2)
  #   x.std.emp = n0/(2*pi)*colSums((beta.k.emp[c(1:v0),])^2)
  # }
  if (v0 == 0) {
    x.std.pow = n0/(2*pi)*(beta.0.pow - sqrt(pi))^2
    x.std.emp = n0/(2*pi)*(beta.0.emp - sqrt(pi))^2
    ###1/3###
    # x.std.pow.3 = n0.sub.3/(2*pi)*(beta.0.pow.3 - sqrt(pi))^2
    # x.std.emp.3 = n0.sub.3/(2*pi)*(beta.0.emp.3 - sqrt(pi))^2
    ###1/5###
    # x.std.pow.5 = n0.sub.5/(2*pi)*(beta.0.pow.5 - sqrt(pi))^2
    # x.std.emp.5 = n0.sub.5/(2*pi)*(beta.0.emp.5 - sqrt(pi))^2
  } else if (v0 == 1) {
    x.std.pow = n0/(2*pi)*(beta.0.pow - sqrt(pi))^2 + n0/(2*pi)*(beta.k.pow[1,])^2
    x.std.emp = n0/(2*pi)*(beta.0.emp - sqrt(pi))^2 + n0/(2*pi)*(beta.k.emp[1,])^2
    ###1/3###
    # x.std.pow.3 = n0.sub.3/(2*pi)*(beta.0.pow.3 - sqrt(pi))^2 + n0.sub.3/(2*pi)*(beta.k.pow.3[1,])^2
    # x.std.emp.3 = n0.sub.3/(2*pi)*(beta.0.emp.3 - sqrt(pi))^2 + n0.sub.3/(2*pi)*(beta.k.emp.3[1,])^2
    ###1/5###
    # x.std.pow.5 = n0.sub.5/(2*pi)*(beta.0.pow.5 - sqrt(pi))^2 + n0.sub.5/(2*pi)*(beta.k.pow.5[1,])^2
    # x.std.emp.5 = n0.sub.5/(2*pi)*(beta.0.emp.5 - sqrt(pi))^2 + n0.sub.5/(2*pi)*(beta.k.emp.5[1,])^2
  } else if (v0 %in% c(2:v)) {
    x.std.pow = n0/(2*pi)*(beta.0.pow - sqrt(pi))^2 + n0/(2*pi)*colSums((beta.k.pow[c(1:v0),])^2)
    x.std.emp = n0/(2*pi)*(beta.0.emp - sqrt(pi))^2 + n0/(2*pi)*colSums((beta.k.emp[c(1:v0),])^2)
    ###1/3###
    # x.std.pow.3 = n0.sub.3/(2*pi)*(beta.0.pow.3 - sqrt(pi))^2 + n0.sub.3/(2*pi)*colSums((beta.k.pow.3[c(1:v0),])^2)
    # x.std.emp.3 = n0.sub.3/(2*pi)*(beta.0.emp.3 - sqrt(pi))^2 + n0.sub.3/(2*pi)*colSums((beta.k.emp.3[c(1:v0),])^2)
    ###1/5###
    # x.std.pow.5 = n0.sub.5/(2*pi)*(beta.0.pow.5 - sqrt(pi))^2 + n0.sub.5/(2*pi)*colSums((beta.k.pow.5[c(1:v0),])^2)
    # x.std.emp.5 = n0.sub.5/(2*pi)*(beta.0.emp.5 - sqrt(pi))^2 + n0.sub.5/(2*pi)*colSums((beta.k.emp.5[c(1:v0),])^2)
  }

  #Comparing to the theoretical critical values.  
  r.rate.1[(v0 + 1)] = 100*length(x.std.pow[x.std.pow > chi.sq.true.1])/N
  r.rate.2[(v0 + 1)] = 100*length(x.std.pow[x.std.pow > chi.sq.true.2])/N
  # ###1/3###
  # r.rate.1.thr.third[(v0 + 1)] = 100*length(x.std.pow.3[x.std.pow.3 > chi.sq.true.1])/N
  # r.rate.2.thr.third[(v0 + 1)] = 100*length(x.std.pow.3[x.std.pow.3 > chi.sq.true.2])/N
  # ###1/5###
  # r.rate.1.thr.fifth[(v0 + 1)] = 100*length(x.std.pow.5[x.std.pow.5 > chi.sq.true.1])/N
  # r.rate.2.thr.fifth[(v0 + 1)] = 100*length(x.std.pow.5[x.std.pow.5 > chi.sq.true.2])/N
  #Comparing to the empirical critical values.
  chi.sq.emp.1 = quantile(x.std.emp, 1 - alpha.1, names = FALSE)
  chi.sq.emp.2 = quantile(x.std.emp, 1 - alpha.2, names = FALSE)
  r.rate.1.emp[(v0 + 1)] = 100*length(x.std.pow[x.std.pow > chi.sq.emp.1])/N
  r.rate.2.emp[(v0 + 1)] = 100*length(x.std.pow[x.std.pow > chi.sq.emp.2])/N
  # ###1/3###
  # chi.sq.emp.1.third = quantile(x.std.emp.3, 1 - alpha.1, names = FALSE)
  # chi.sq.emp.2.third = quantile(x.std.emp.3, 1 - alpha.2, names = FALSE)
  # r.rate.1.emp.third[(v0 + 1)] = 100*length(x.std.pow.3[x.std.pow.3 > chi.sq.emp.1.third])/N
  # r.rate.2.emp.third[(v0 + 1)] = 100*length(x.std.pow.3[x.std.pow.3 > chi.sq.emp.2.third])/N
  # ###1/5###
  # chi.sq.emp.1.fifth = quantile(x.std.emp.5, 1 - alpha.1, names = FALSE)
  # chi.sq.emp.2.fifth = quantile(x.std.emp.5, 1 - alpha.2, names = FALSE)
  # r.rate.1.emp.fifth[(v0 + 1)] = 100*length(x.std.pow.5[x.std.pow.5 > chi.sq.emp.1.fifth])/N
  # r.rate.2.emp.fifth[(v0 + 1)] = 100*length(x.std.pow.5[x.std.pow.5 > chi.sq.emp.2.fifth])/N
  
  # name = paste("Altr_DF_",v0 + 1,".png",sep="")
  # png(name)
  # if (v0 == 0) {
  #   hist(x.std.pow,prob = TRUE,nclass = 20,xlab = "X[m]",main = "Density estimate of Data")
  # } else {
  #   ylim_upper = max(dchisq(seq(0,80,.01),df = v0 + 1))
  #   hist(x.std.pow,prob = TRUE,nclass = 20,ylim = c(0,ylim_upper + .005),xlab = "X[m]",main = "Density estimate of Data")
  # }
  # lines(seq(0,80,.01),dchisq(seq(0,80,.01),df = v0 + 1),col = "red")
  # abline(v = c(chi.sq.true.1,chi.sq.true.2,chi.sq.emp.1,chi.sq.emp.2),col = c("blue","blue","purple","purple"))
  # dev.off()
  # 
  # name = paste("Null_DF_",v0 + 1,".png",sep="")
  # png(name)
  # hist(x.std.emp,prob = TRUE,nclass = 20,xlab = "X[m]",main = "Density estimate of Data")
  # lines(seq(0,80,.01),dchisq(seq(0,80,.01),df = v0 + 1),col = "red")
  # abline(v = c(chi.sq.true.1,chi.sq.true.2),col = c("blue","purple"))
  # dev.off()
  
}
r.rate.1.c = matrix(r.rate.1,nrow = (v + 1),ncol = 1)
r.rate.2.c = matrix(r.rate.2,nrow = (v + 1),ncol = 1)

r.rate.1.emp.c = matrix(r.rate.1.emp,nrow = (v + 1),ncol = 1)
r.rate.2.emp.c = matrix(r.rate.2.emp,nrow = (v + 1),ncol = 1)

# ###1/3###
# r.rate.1.thr.third.c = matrix(r.rate.1.thr.third,nrow = (v + 1),ncol = 1)
# r.rate.2.thr.third.c = matrix(r.rate.2.thr.third,nrow = (v + 1),ncol = 1)
# 
# r.rate.1.emp.third.c = matrix(r.rate.1.emp.third,nrow = (v + 1),ncol = 1)
# r.rate.2.emp.third.c = matrix(r.rate.2.emp.third,nrow = (v + 1),ncol = 1)
# ###1/5###
# r.rate.1.thr.fifth.c = matrix(r.rate.1.thr.fifth,nrow = (v + 1),ncol = 1)
# r.rate.2.thr.fifth.c = matrix(r.rate.2.thr.fifth,nrow = (v + 1),ncol = 1)
# 
# r.rate.1.emp.fifth.c = matrix(r.rate.1.emp.fifth,nrow = (v + 1),ncol = 1)
# r.rate.2.emp.fifth.c = matrix(r.rate.2.emp.fifth,nrow = (v + 1),ncol = 1)