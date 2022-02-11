#Simulation
library(astsa)
library(stats)
library(arfima)
library(fracdiff)
require(svMisc)

setwd("H:/TimeSeriesTest/LiAndBiao/Beta0_Additional(MLE)")

#Size of the sample
n0 = 1024

#The number of beta.k, we have beta.1, beta.2, ...
v = 30
k0 = c(1:v)

#The Fourier Frequencies
m = n0/2 - 1
j = c(1:m)
mu.j = 2*pi*j/n0

###1/3###
mu.j.filter = mu.j[seq(2,m,3)]
sub.length = m/3
# n0.sub = 2*(sub.length + 1)

###1/5###
# mu.j.filter = mu.j[seq(3,m,5)]
# sub.length = m/5
# n0.sub = 2*(sub.length + 1)

#The number of the simulation
N = 5000

#beta.0.ar1 is used to store all beta.0.
#beta.k.ar1 is used to store all beta.k, this is a 6-row, 5000-column matrix.
beta.0 = rep(0,N)
beta.k = matrix(0,nrow = v,ncol = N)
#The following two are for the cases when parameters are known.
beta.0.true = rep(0,N)
beta.k.true = matrix(0,nrow = v,ncol = N)
#1/3
beta.0.Avg = rep(0,N)
beta.k.Avg = matrix(0,nrow = v,ncol = N)
#The following two are for the cases when parameters are known.
beta.0.Avg.true = rep(0,N)
beta.k.Avg.true = matrix(0,nrow = v,ncol = N)

#Spectral Density Function
#@parameter w.unit: this is our fourier frequency
#@parameter phi.1: this is the coefficient of AR(1), the true value is 0.8
#@parameter sig2_f: the variance of noises, notated as sigma square.
####Short Memory####
sp.short.f.est = function(w.unit,phi.1,sig2_f){(sig2_f/((1 + phi.1^2) - 2*phi.1*cos(w.unit)))/(2*pi)}
# Table 1
sp.short.true = sp.short.f.est(w.unit = mu.j, phi.1 = .8, sig2_f = 1)
# Table 3
# sp.short.true = sp.short.f.est(w.unit = mu.j, phi.1 = .8, sig2_f = 9/7)
####Additional####
# ARMA(1,1)
# sp.short.f.est = function(w.loc,phi.loc,theta.loc,sig2_f){sig2_f/(2*pi)*(1 + theta.loc^2 + 2*theta.loc*cos(w.loc))/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}
# sp.short.true = sp.short.f.est(w.loc = mu.j, phi.loc = .8, theta.loc = -.15, sig2_f = 1)
# AR(2) P179
# sp.short.f.est = function(w.unit,phi.1,phi.2,sig2_f){(sig2_f/((1 + phi.1^2 + phi.2^2) + 2*(phi.1*phi.2 - phi.1)*cos(w.unit) - 2*phi.2*cos(2*w.unit)))/(2*pi)}
# sp.short.true = sp.short.f.est(w.unit = mu.j, phi.1 = 1, phi.2 = -.9, sig2_f = 1)
# MA(2)
# sp.short.f.est = function(w.loc,theta1.loc,theta2.loc,sig2_f){sig2_f/(2*pi)*(1 + theta1.loc^2 + theta2.loc^2 + 2*(theta1.loc*theta2.loc + theta1.loc)*cos(w.loc) + 2*theta2.loc*cos(2*w.loc))}
# sp.short.true = sp.short.f.est(w.loc = mu.j, theta1.loc = .9, theta2.loc = .18, sig2_f = 1)
# ARMA(2,2)
# sp.short.f.est = function(w.loc,phi1.loc,phi2.loc,theta1.loc,theta2.loc,sig2_f){sig2_f/(2*pi)*(1 + theta1.loc^2 + theta2.loc^2 + 2*(theta1.loc*theta2.loc + theta1.loc)*cos(w.loc) + 2*theta2.loc*cos(2*w.loc))/(1 + phi1.loc^2 + phi2.loc^2 + 2*(phi1.loc*phi2.loc - phi1.loc)*cos(w.loc) - 2*phi2.loc*cos(2*w.loc))}
# sp.short.true = sp.short.f.est(w.loc = mu.j, phi1.loc = -1, phi2.loc = -.25, theta1.loc = 1, theta2.loc = .24, sig2_f = 1)


####Long Memory####
# sp.long.f.est = function(w.loc,d.loc,sig2_f){sig2_f*(4*(sin(w.loc/2))^2)^(-d.loc)/(2*pi)}
# Table 2
# sp.long.true = sp.long.f.est(w.loc = mu.j,d.loc = .4,sig2_f = 1)
# Table 4
# sp.long.true = sp.long.f.est(w.loc = mu.j,d.loc = .4,sig2_f = 9/7)
####Additional####
# ARFIMA(1,1)
# sp.long.f.est = function(w.loc,phi.loc,theta.loc,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)*(1 + theta.loc^2 + 2*theta.loc*cos(w.loc))/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}
# sp.long.true = sp.long.f.est(w.loc = mu.j, phi.loc = .6, theta.loc = -.25, d.loc = .3, sig2_f = 1)
# ARFIMA(1,d,0)
# sp.long.f.est = function(w.loc,phi.loc,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}
# sp.long.true = sp.long.f.est(w.loc = mu.j, phi.loc = -.5, d.loc = .4, sig2_f = 1)
# ARFIMA(2,d,0) 
# sp.long.f.est = function(w.loc,phi.1,phi.2,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)/((1 + phi.1^2 + phi.2^2) + 2*(phi.1*phi.2 - phi.1)*cos(w.loc) - 2*phi.2*cos(2*w.loc))}
# sp.long.true = sp.long.f.est(w.loc = mu.j, phi.1 = -1.2, phi.2 = -.32, d.loc = .3, sig2_f = 1)
# ARFIMA(0,d,2)
# sp.long.f.est = function(w.loc,theta1.loc,theta2.loc,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)*(1 + theta1.loc^2 + theta2.loc^2 + 2*(theta1.loc*theta2.loc + theta1.loc)*cos(w.loc) + 2*theta2.loc*cos(2*w.loc))}
# sp.long.true = sp.long.f.est(w.loc = mu.j, theta1.loc = -1.2, theta2.loc = .27, d.loc = .3, sig2_f = 1)
# ARFIMA(2,d,2)
# sp.long.f.est = function(w.loc,phi1.loc,phi2.loc,theta1.loc,theta2.loc,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)*(1 + theta1.loc^2 + theta2.loc^2 + 2*(theta1.loc*theta2.loc + theta1.loc)*cos(w.loc) + 2*theta2.loc*cos(2*w.loc))/(1 + phi1.loc^2 + phi2.loc^2 + 2*(phi1.loc*phi2.loc - phi1.loc)*cos(w.loc) - 2*phi2.loc*cos(2*w.loc))}
# sp.long.true = sp.long.f.est(w.loc = mu.j, phi1.loc = -.7, phi2.loc = -.1, theta1.loc = -.65, theta2.loc = .105, d.loc = .3, sig2_f = 1)

#Beta.0 Function
#It is used to calculate beta.0
#@parameter I: The Periodogram
#@parameter sp: The Spectral Density
beta.0.f = function(I,sp){
  sum(I/sp)*2*sqrt(pi)/(2*(length(I) + 1))
}

#Beta.k Function
#It is used to calculate beta.k
#@parameter I: The Periodogram
#@parameter sp: The Spectral Density
beta.k.f = function(k,I,sp,mu.j.in){
  sum(cos(k*mu.j.in)*I/sp)*2*sqrt(2*pi)/(2*(length(I) + 1))
}

#r.rate.1 is to store the rejection rate at 5%
#r,rate.2 is to store the rejection rate at 10%
r.rate.1 = rep(0,(v + 1))
r.rate.2 = rep(0,(v + 1))
#The following two are for the cases when parameters are known.
r.rate.1.true = rep(0,(v + 1))
r.rate.2.true = rep(0,(v + 1))
#1/3
r.rate.1.Avg = rep(0,(v + 1))
r.rate.2.Avg = rep(0,(v + 1))
#The following two are for the cases when parameters are known.
r.rate.1.Avg.true = rep(0,(v + 1))
r.rate.2.Avg.true = rep(0,(v + 1))


#phi.est.v is to store the coeff est in AR(1)
phi.est.v = rep(0,N)
# phi1.est.v = rep(0,N)
# phi2.est.v = rep(0,N)
# d.est.v = rep(0,N)
# theta.est.v = rep(0,N)
# theta1.est.v = rep(0,N)
# theta2.est.v = rep(0,N)
sig2.est.v = rep(0,N)


###############
##Simulations##
###############


    
    for(sim in c(1:N)){
      ####### Table 1 & 3###### For t, use innov = rt(n0,9)
      x.sim.short = arima.sim(list(order = c(1,0,0),ar = .8), n = n0)#, innov = rt(n0,9))

      # Periodogram
      # fft function returns a vector, contains 512 complex numbers,
      # they are corresponding to the frequencies 2*pi*0/n, 2*pi*1/n, ... , 2*pi*(n - 1)/n.
      # Here we are looking for frequencies 2*pi*1/n, ... , 2*pi*(n/2 - 1)/n,
      # therefore we pick up sub vector whose index is from 2 to n/2.
      # Details could be checked in fft function Documents.
      prd = Mod(fft(x.sim.short- mean(x.sim.short)))^2/(2*pi*n0)
      #The periodogram starting at 1
      x.prd.short = prd[2:(n0/2)]
      
      #######Additional######
      # simulation of ARMA(1,0,1)
      # x.sim.short = arima.sim(list(order = c(1,0,1),ar = .8,ma = -.15),n = n0)
      # prd.arma = Mod(fft(x.sim.short- mean(x.sim.short)))^2/(2*pi*n0)
      # x.prd.short = prd.arma[2:(n0/2)]
      # simulation of AR(2)
      # x.sim.short = arima.sim(list(order = c(2,0,0),ar = c(1, -.9)),n = n0)
      # prd.arma = Mod(fft(x.sim.short- mean(x.sim.short)))^2/(2*pi*n0)
      # x.prd.short = prd.arma[2:(n0/2)]
      # simulation of MA(2)
      # x.sim.short = arima.sim(list(order = c(0,0,2),ma = c(.9, .18)),n = n0)
      # prd.arma = Mod(fft(x.sim.short- mean(x.sim.short)))^2/(2*pi*n0)
      # x.prd.short = prd.arma[2:(n0/2)]
      # simulation of ARMA(2,0,2)
      # x.sim.short = arima.sim(list(order = c(2,0,2), ar = c(-1,-.25),ma = c(1,.24)),n = n0)
      # prd.arma = Mod(fft(x.sim.short- mean(x.sim.short)))^2/(2*pi*n0)
      # x.prd.short = prd.arma[2:(n0/2)]
      
      ####### Table 2 & 4###### For t, use innov = rt(n0,9)
      # x.sim.long = (fracdiff.sim(n0,d = 0.4))$series
      # prd = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
      # x.prd.long = prd[2:(n0/2)]
      
      #######Additional######
      # simulation of ARFIMA(1,d,1)
      # x.sim.long = (fracdiff.sim(n0,ar = .6, ma = .25,d = .3))$series
      # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
      # x.prd.long = prd.long[2:(n0/2)]
      # simulation of ARFIMA(1,d,0)
      # x.sim.long = (fracdiff.sim(n0,ar = -.5,d = .4))$series
      # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
      # x.prd.long = prd.long[2:(n0/2)]
      # simulation of ARFIMA(2,d,0)
      # x.sim.long = (fracdiff.sim(n0,ar = c(-1.2,-.32),d = .3))$series
      # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
      # x.prd.long = prd.long[2:(n0/2)]
      # simulation of ARFIMA(0,d,2), ma is in the opposite sign ! ! !
      # x.sim.long = (fracdiff.sim(n0,ma = c(1.2,-.27),d = .3))$series
      # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
      # x.prd.long = prd.long[2:(n0/2)]
      # simulation of ARFIMA(2,0,2), ma is in the opposite sign ! ! 
      # x.sim.long = (fracdiff.sim(n0,ar = c(-.7,-.1), ma = c(.65,-.105),d = .3))$series
      # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
      # x.prd.long = prd.long[2:(n0/2)]
      
      #######MLE Table 1 & 3##########
      x.sim.size = arima0(x.sim.short,order = c(1,0,0),method = "ML",include.mean = FALSE)
      phi.size = x.sim.size$coef

      attributes(phi.size) = NULL

      sig2.size = x.sim.size$sigma2
      phi.est.v[sim] = phi.size
      sig2.est.v[sim] = sig2.size
      sp.short.est = sp.short.f.est(w.unit = mu.j,phi.1 = phi.size,sig2_f = 1)#sig2.size
      ################################
      
      #######MLE Table 2 & 4##########
      # x.MLE = arfima(x.sim.long,order = c(0,0,0),quiet = TRUE)$modes
      # d_l = x.MLE[[1]]$dfrac
      # sig2.est = x.MLE[[1]]$sigma2
      # d.est.v[sim] = d_l
      # sig2.est.v[sim] = sig2.est
      # sp.long.est = sp.long.f.est(w.loc = mu.j,d.loc = d_l,sig2_f = 1)#sig2.est
      ################################
      
      #######MLE Additional ARMA(1,1)##########
      # x.sim = arima0(x.sim.short,order = c(1,0,1),method = "ML",include.mean = FALSE)
      # # x.sim = sarima(x.sim.short,1,0,1,details = FALSE,no.constant = TRUE)$fit
      # phi.est = x.sim$coef[1]
      # theta.est = x.sim$coef[2]
      # sig2.est = x.sim$sigma2
      # 
      # #For sarima function
      # # attributes(phi.est) = NULL
      # # attributes(theta.est) = NULL
      # 
      # phi.est.v[sim] = phi.est
      # theta.est.v[sim] = theta.est
      # sig2.est.v[sim] = sig2.est
      # sp.short.est = sp.short.f.est(w.loc = mu.j,phi.loc = phi.est,theta.loc = theta.est,sig2_f = sig2.est)
      ##########################################
      
      #######MLE Additional ARFIMA(1,d,0)##########
      # x.sim = arfima(x.sim.long,order = c(1,0,0),quiet = TRUE)$modes
      # d.est = x.sim[[1]]$dfrac
      # phi.est = x.sim[[1]]$phi
      # sig2.est = x.sim[[1]]$sigma2
      # 
      # d.est.v[sim] = d.est
      # phi.est.v[sim] = phi.est
      # sig2.est.v[sim] = sig2.est
      # sp.long.est = sp.long.f.est(w.loc = mu.j,phi.loc = phi.est,d.loc = d.est,sig2_f = sig2.est)
      ##########################################
      
      #######MLE Additional AR(2)##########
      # x.sim.size = arima0(x.sim.short,order = c(2,0,0),method = "ML",include.mean = FALSE)
      # x.sim.size = sarima(x.sim.short,2,0,0,details = FALSE,no.constant = TRUE)$fit
      # phi1.size = x.sim.size$coef[1]
      # phi2.size = x.sim.size$coef[2]
      # 
      # attributes(phi1.size) = NULL
      # attributes(phi2.size) = NULL
      # 
      # sig2.size = x.sim.size$sigma2
      # phi1.est.v[sim] = phi1.size
      # phi2.est.v[sim] = phi2.size
      # sig2.est.v[sim] = sig2.size
      # sp.short.est = sp.short.f.est(w.unit = mu.j,phi.1 = phi1.size,phi.2 = phi2.size,sig2_f = sig2.size)
      ######################################
      
      #######MLE Additional ARMA(2,2)##########
      # x.sim.size = arima0(x.sim.short,order = c(2,0,2),method = "ML",include.mean = FALSE)
      # # x.sim.size = sarima(x.sim.short,2,0,2,details = FALSE,no.constant = TRUE)$fit
      # phi1.size = x.sim.size$coef[1]
      # phi2.size = x.sim.size$coef[2]
      # theta1.size = x.sim.size$coef[3]
      # theta2.size = x.sim.size$coef[4]
      # 
      # attributes(phi1.size) = NULL
      # attributes(phi2.size) = NULL
      # attributes(theta1.size) = NULL
      # attributes(theta2.size) = NULL
      # 
      # sig2.size = x.sim.size$sigma2
      # phi1.est.v[sim] = phi1.size
      # phi2.est.v[sim] = phi2.size
      # theta1.est.v[sim] = theta1.size
      # theta2.est.v[sim] = theta2.size
      # sig2.est.v[sim] = sig2.size
      # sp.short.est = sp.short.f.est(w.loc = mu.j, phi1.loc = phi1.size, phi2.loc = phi2.size, theta1.loc = theta1.size, theta2.loc = theta2.size, sig2_f = sig2.size)
      ######################################
      
      ###Avg Periodogram###
      x.prd.avg = rep(0,sub.length)
      k_s = 0
      for (j_s in seq(1,m,3)) {
        k_s = k_s + 1
        x.prd.avg[k_s] = (x.prd.short[j_s] + x.prd.short[j_s + 1] + x.prd.short[j_s + 2])/3
      }
      
      ###Filter Periodogram###
      # x.prd.filter = x.prd.long[seq(2,m,3)]
      
      ###Filter Spectral Density###
      sp.filter.est = sp.short.est[seq(2,m,3)]
      sp.filter.true = sp.short.true[seq(2,m,3)]
      
      # Calculate beta.0.
      # This is for the case pars are unkown, please refer doc for details#######
      beta.0[sim] = beta.0.f(I = x.prd.short,sp = sp.short.est)
      beta.0.true[sim] = beta.0.f(I = x.prd.short,sp = sp.short.true)
      #1/3
      beta.0.Avg[sim] = beta.0.f(I = x.prd.avg,sp = sp.filter.est)
      beta.0.Avg.true[sim] = beta.0.f(I = x.prd.avg,sp = sp.filter.true)
      ###########################################################################
      
      #Calculate beta1, beta2, ...
      for (i in k0) {
        beta.k[i,sim] = beta.k.f(k = i,I = x.prd.short, sp = sp.short.est, mu.j.in = mu.j)
        beta.k.true[i,sim] = beta.k.f(k = i,I = x.prd.short, sp = sp.short.true, mu.j.in = mu.j)
        #1/3
        beta.k.Avg[i,sim] = beta.k.f(k = i,I = x.prd.avg, sp = sp.filter.est, mu.j.in = mu.j.filter)
        beta.k.Avg.true[i,sim] = beta.k.f(k = i,I = x.prd.avg, sp = sp.filter.true, mu.j.in = mu.j.filter)
      }
      progress(sim, max.value = N)
      if (sim == N) cat("Done!\n")
    }

for (v0 in c(0:v)) {
  alpha.1 = .05
  alpha.2 = .1

  # when sig2 is unknown, beta.0 is excluded, df = v0
  chi.sq.true.1 = qchisq(alpha.1,df = v0 + 1,lower.tail = FALSE)
  chi.sq.true.2 = qchisq(alpha.2,df = v0 + 1,lower.tail = FALSE)

  if (v0 == 0) {
    # If it's t, here should be changed ! ! ! q0 = 1 + 243/245 = 488/245
    x.std = n0/(2*pi)*(beta.0 - sqrt(pi))^2
    x.std.true = n0/(2*pi)*(beta.0.true - sqrt(pi))^2
    #1/3
    x.std.Avg = n0/(2*pi)*(beta.0.Avg - sqrt(pi))^2
    x.std.Avg.true = n0/(2*pi)*(beta.0.Avg.true - sqrt(pi))^2
  } else if (v0 == 1) {
    # If it's t, here should be changed ! ! ! q0 = 1 + 243/245 = 488/245
    x.std = n0/(2*pi)*(beta.0 - sqrt(pi))^2 + n0/(2*pi)*(beta.k[1,])^2
    x.std.true = n0/(2*pi)*(beta.0.true - sqrt(pi))^2 + n0/(2*pi)*(beta.k.true[1,])^2
    #1/3
    x.std.Avg = n0/(2*pi)*(beta.0.Avg - sqrt(pi))^2 + n0/(2*pi)*(beta.k.Avg[1,])^2
    x.std.Avg.true = n0/(2*pi)*(beta.0.Avg.true - sqrt(pi))^2 + n0/(2*pi)*(beta.k.Avg.true[1,])^2
  } else if (v0 %in% c(2:v)) {
    # If it's t, here should be changed ! ! ! q0 = 1 + 243/245 = 488/245
    x.std = n0/(2*pi)*(beta.0 - sqrt(pi))^2 + n0/(2*pi)*colSums((beta.k[c(1:v0),])^2)
    x.std.true = n0/(2*pi)*(beta.0.true - sqrt(pi))^2 + n0/(2*pi)*colSums((beta.k.true[c(1:v0),])^2)
    #1/3
    x.std.Avg = n0/(2*pi)*(beta.0.Avg - sqrt(pi))^2 + n0/(2*pi)*colSums((beta.k.Avg[c(1:v0),])^2)
    x.std.Avg.true = n0/(2*pi)*(beta.0.Avg.true - sqrt(pi))^2 + n0/(2*pi)*colSums((beta.k.Avg.true[c(1:v0),])^2)
  }

  # name = paste("Est_DF_",v0 + 1,".png",sep="")
  # png(name)
  # hist(x.std,prob = TRUE,nclass = 20,xlab = "X[m]",main = "Density estimate of Data")
  # lines(seq(0,80,.01),dchisq(seq(0,80,.01),df = v0 + 1),col = "red")
  # abline(v = c(chi.sq.true.1,chi.sq.true.2),col = c("blue","red"))
  # dev.off()
  #
  # name = paste("TRUE_DF_",v0 + 1,".png",sep="")
  # png(name)
  # hist(x.std.true,prob = TRUE,nclass = 20,xlab = "X[m]",main = "Density estimate of Data")
  # lines(seq(0,80,.01),dchisq(seq(0,80,.01),df = v0 + 1),col = "red")
  # abline(v = c(chi.sq.true.1,chi.sq.true.2),col = c("blue","red"))
  # dev.off()

  #Pars are unkown
  r.rate.1[(v0 + 1)] = 100*length(x.std[x.std > chi.sq.true.1])/N
  r.rate.2[(v0 + 1)] = 100*length(x.std[x.std > chi.sq.true.2])/N
  #Pars are given
  r.rate.1.true[(v0 + 1)] = 100*length(x.std.true[x.std.true > chi.sq.true.1])/N
  r.rate.2.true[(v0 + 1)] = 100*length(x.std.true[x.std.true > chi.sq.true.2])/N

  #Pars are unkown
  r.rate.1.Avg[(v0 + 1)] = 100*length(x.std.Avg[x.std.Avg > chi.sq.true.1])/N
  r.rate.2.Avg[(v0 + 1)] = 100*length(x.std.Avg[x.std.Avg > chi.sq.true.2])/N
  #Pars are given
  r.rate.1.Avg.true[(v0 + 1)] = 100*length(x.std.Avg.true[x.std.Avg.true > chi.sq.true.1])/N
  r.rate.2.Avg.true[(v0 + 1)] = 100*length(x.std.Avg.true[x.std.Avg.true > chi.sq.true.2])/N
}

r.rate.1.c = matrix(r.rate.1,nrow = (v + 1),ncol = 1)
r.rate.2.c = matrix(r.rate.2,nrow = (v + 1),ncol = 1)

r.rate.1.true.c = matrix(r.rate.1.true,nrow = (v + 1),ncol = 1)
r.rate.2.true.c = matrix(r.rate.2.true,nrow = (v + 1),ncol = 1)
#1/3
r.rate.1.Avg.c = matrix(r.rate.1.Avg,nrow = (v + 1),ncol = 1)
r.rate.2.Avg.c = matrix(r.rate.2.Avg,nrow = (v + 1),ncol = 1)

r.rate.1.Avg.true.c = matrix(r.rate.1.Avg.true,nrow = (v + 1),ncol = 1)
r.rate.2.Avg.true.c = matrix(r.rate.2.Avg.true,nrow = (v + 1),ncol = 1)

#####################################################
##The following part is to verify Theorem 2.4      ##
#####################################################
##Note: The following stats are for analysis, they
##are obtained in two cases, sig2 is unknown, and 
##sig2 is given as 1.

# Combine the standardized beta.0 and beta.k together.
# This gives us a 7-row, 5000-column matrix.
beta.st = rbind(sqrt(n0)*(beta.0 - sqrt(pi)),beta.k*sqrt(n0))
beta.st.true = rbind(sqrt(n0)*(beta.0.true - sqrt(pi)),beta.k.true*sqrt(n0))
#1/3
beta.st.Avg = rbind(sqrt(n0)*(beta.0.Avg - sqrt(pi)),beta.k.Avg*sqrt(n0))
beta.st.Avg.true = rbind(sqrt(n0)*(beta.0.Avg.true - sqrt(pi)),beta.k.Avg.true*sqrt(n0))
# Find the covaraince of all beta.
cov.beta = cov(t(beta.st))
cov.beta.true = cov(t(beta.st.true))
#1/3
cov.beta.Avg = cov(t(beta.st.Avg))
cov.beta.Avg.true = cov(t(beta.st.Avg.true))
# The correlation matrix to store the correlations among beta.
cor.beta = matrix(rep(0,(1 + v)^2),nrow = 1 + v,ncol = 1 + v)
cor.beta.true = matrix(rep(0,(1 + v)^2),nrow = 1 + v,ncol = 1 + v)

cor.beta.Avg = matrix(rep(0,(1 + v)^2),nrow = 1 + v,ncol = 1 + v)
cor.beta.Avg.true = matrix(rep(0,(1 + v)^2),nrow = 1 + v,ncol = 1 + v)
# Find the corrleation of all beta,
# this gives you an upper right triangular matrix, 
# the dignal stores variances.
for(i in c(1:(1 + v))){
  for(j in c(i:(1 + v))){
    
    if(i == j){
      cor.beta[i,j] = cov.beta[i,j]
      cor.beta.true[i,j] = cov.beta.true[i,j]
      #1/3
      cor.beta.Avg[i,j] = cov.beta.Avg[i,j]
      cor.beta.Avg.true[i,j] = cov.beta.Avg.true[i,j]
    } else {
      # if(i != 1){
      cor.beta[i,j] = cov.beta[i,j]/sqrt(cov.beta[i,i]*cov.beta[j,j])
      # }
      cor.beta.true[i,j] = cov.beta.true[i,j]/sqrt(cov.beta.true[i,i]*cov.beta.true[j,j])
      #1/3
      cor.beta.Avg[i,j] = cov.beta.Avg[i,j]/sqrt(cov.beta.Avg[i,i]*cov.beta.Avg[j,j])
      cor.beta.Avg.true[i,j] = cov.beta.Avg.true[i,j]/sqrt(cov.beta.Avg.true[i,i]*cov.beta.Avg.true[j,j])
    } 
    
  }
}

# Fine the mean of the standardized beta.0 and beta.k.
# They should be close to 0.
u = sapply(c(1:v),function(i){mean(sqrt(n0)*beta.k[i,])})
u = matrix(c(mean(sqrt(n0)*(beta.0 - sqrt(pi))),u),nrow = 1)

u.true = sapply(c(1:v),function(i){mean(sqrt(n0)*beta.k.true[i,])})
u.true = matrix(c(mean(sqrt(n0)*(beta.0.true - sqrt(pi))),u.true),nrow = 1)

hist(beta.0, nclass = 20, prob = TRUE, main = "Sig2 is unknown")
hist(sqrt(n0)*(beta.0 - sqrt(pi)),xlim = c(-6,6), nclass = 20, prob = TRUE, main = "Sig2 is unknown")
lines(seq(-10,10,.01),dnorm(seq(-10,10,.01)),col = "red")

hist(sqrt(n0)*beta.k[1,], nclass = 20, xlim = c(-10,10), prob = TRUE, main = "Sig2 is unknown")
lines(seq(-10,10,.01),dnorm(seq(-10,10,.01),sd = sqrt(2*pi)),col = "red")

hist(sqrt(n0)*(beta.0.true - sqrt(pi)), nclass = 20, prob = TRUE, main = "Sig2 is unknown")
lines(seq(-10,10,.01),dnorm(seq(-10,10,.01),sd = sqrt(2*pi)),col = "red")

#Check the dist of the est phi
hist(phi.est.v, nclass = 20, prob = TRUE, main = "Sig2 is unknown")
hist(phi.est.sig2.v, nclass = 20, prob = TRUE, main = "Sig2 is fixed")