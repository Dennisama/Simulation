library(fracdiff)
library(astsa)
library(arfima)
require(svMisc)
#512, 1046
n0 = 512
N = 5000

#####################Chen And Deo#####################
j = c(1:(n0-1))
mu.j = 2*pi*j/n0

#####################Li And Biao######################
v = 30
k0 = c(1:v)

m = n0/2 - 1
j.LB = c(1:m)
mu.j.LB = 2*pi*j.LB/n0

mu.j.LB.filter = mu.j.LB[seq(2,m,3)]
sub.length = m/3

#####################Chen And Deo#####################
Bartlett = function(z){
  ifelse(abs(z) <= 1,1 - abs(z),0)
}

Tukey = function(z){
  ifelse(abs(z) <= 1,.5*(cos(z*pi) + 1),0)
}

QS = function(z){
  ifelse(z == 0,1,(sin(6*pi*z/5)/(6*pi*z/5) - cos(6*pi*z/5))*25/(12*pi^2*z^2)) 
}

#####Bartlett#####
W_BAR_37.fft = function(j_f){
  1/(2*pi)*(fft(Bartlett(c(0:(n0-1))/37)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Bartlett(c(0:(n0-1))/37)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Bartlett(0))
}
w_BAR_37.fft.test = sapply(j, W_BAR_37.fft)
w_BAR_37.fft.test = t(Re(w_BAR_37.fft.test))

W_BAR_20.fft = function(j_f){
  1/(2*pi)*(fft(Bartlett(c(0:(n0-1))/20)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Bartlett(c(0:(n0-1))/20)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Bartlett(0))
}
w_BAR_20.fft.test = sapply(j, W_BAR_20.fft)
w_BAR_20.fft.test = t(Re(w_BAR_20.fft.test))

W_BAR_11.fft = function(j_f){
  1/(2*pi)*(fft(Bartlett(c(0:(n0-1))/11)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Bartlett(c(0:(n0-1))/11)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Bartlett(0))
}
w_BAR_11.fft.test = sapply(j, W_BAR_11.fft)
w_BAR_11.fft.test = t(Re(w_BAR_11.fft.test))
##################


######Tukey#######
W_TUK_37.fft = function(j_f){
  1/(2*pi)*(fft(Tukey(c(0:(n0-1))/37)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Tukey(c(0:(n0-1))/37)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Tukey(0))
}
w_TUK_37.fft.test = sapply(j, W_TUK_37.fft)
w_TUK_37.fft.test = t(Re(w_TUK_37.fft.test))

W_TUK_20.fft = function(j_f){
  1/(2*pi)*(fft(Tukey(c(0:(n0-1))/20)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Tukey(c(0:(n0-1))/20)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Tukey(0))
}
w_TUK_20.fft.test = sapply(j, W_TUK_20.fft)
w_TUK_20.fft.test = t(Re(w_TUK_20.fft.test))

W_TUK_11.fft = function(j_f){
  1/(2*pi)*(fft(Tukey(c(0:(n0-1))/11)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Tukey(c(0:(n0-1))/11)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Tukey(0))
}
w_TUK_11.fft.test = sapply(j, W_TUK_11.fft)
w_TUK_11.fft.test = t(Re(w_TUK_11.fft.test))
##################

########QS########
W_QS_37.fft = function(j_f){
  1/(2*pi)*(fft(QS(c(0:(n0-1))/37)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(QS(c(0:(n0-1))/37)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -QS(0))
}
w_QS_37.fft.test = sapply(j, W_QS_37.fft)
w_QS_37.fft.test = t(Re(w_QS_37.fft.test))

W_QS_20.fft = function(j_f){
  1/(2*pi)*(fft(QS(c(0:(n0-1))/20)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(QS(c(0:(n0-1))/20)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -QS(0))
}
w_QS_20.fft.test = sapply(j, W_QS_20.fft)
w_QS_20.fft.test = t(Re(w_QS_20.fft.test))

W_QS_11.fft = function(j_f){
  1/(2*pi)*(fft(QS(c(0:(n0-1))/11)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(QS(c(0:(n0-1))/11)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -QS(0))
}
w_QS_11.fft.test = sapply(j, W_QS_11.fft)
w_QS_11.fft.test = t(Re(w_QS_11.fft.test))
##################

#####################Li And Biao######################
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


####Short Memory####
sp.short.f.est = function(w.unit,phi.1,sig2_f){(sig2_f/((1 + phi.1^2) - 2*phi.1*cos(w.unit)))/(2*pi)}
# Table 1
sp.short.true = sp.short.f.est(w.unit = mu.j, phi.1 = .8, sig2_f = 1)
sp.short.true.LB = sp.short.f.est(w.unit = mu.j.LB, phi.1 = .8, sig2_f = 1)
# Table 3
# sp.short.true = sp.short.f.est(w.unit = mu.j, phi.1 = .8, sig2_f = 9/7)
##Long Memory####
# sp.long.f.est = function(w.loc,d.loc,sig2_f){sig2_f*(4*(sin(w.loc/2))^2)^(-d.loc)/(2*pi)}
# Table 2
# sp.long.true = sp.long.f.est(w.loc = mu.j,d.loc = .4,sig2_f = 1)
# Table 4
# sp.long.true = sp.long.f.est(w.loc = mu.j,d.loc = .4,sig2_f = 9/7)

d.est.v = rep(0,N)
phi.est.v = rep(0,N)
sig2.est.v = rep(0,N)

alpha.1 = .05
alpha.2 = .1

#####################Chen And Deo#####################
###BAR###
Tn.BAR_37 = rep(0,N)
Tn_q.BAR_37 = rep(0,N)

Tn.BAR_20 = rep(0,N)
Tn_q.BAR_20 = rep(0,N)

Tn.BAR_11 = rep(0,N)
Tn_q.BAR_11 = rep(0,N)
#########

###TUK###
Tn.TUK_37 = rep(0,N)
Tn_q.TUK_37 = rep(0,N)

Tn.TUK_20 = rep(0,N)
Tn_q.TUK_20 = rep(0,N)

Tn.TUK_11 = rep(0,N)
Tn_q.TUK_11 = rep(0,N)
#########

###QS###
Tn.QS_37 = rep(0,N)
Tn_q.QS_37 = rep(0,N)

Tn.QS_20 = rep(0,N)
Tn_q.QS_20 = rep(0,N)

Tn.QS_11 = rep(0,N)
Tn_q.QS_11 = rep(0,N)
#########

####BAR####
Cn_BAR_37 = sum((1 - c(1:(n0-1))/n0)*Bartlett(c(1:(n0-1))/37)^2)/(n0*pi) + 1/(2*pi)
Dn_BAR_37 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Bartlett(c(1:(n0-2))/37)^4)*2/pi^2

Cn_BAR_20 = sum((1 - c(1:(n0-1))/n0)*Bartlett(c(1:(n0-1))/20)^2)/(n0*pi) + 1/(2*pi)
Dn_BAR_20 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Bartlett(c(1:(n0-2))/20)^4)*2/pi^2

Cn_BAR_11 = sum((1 - c(1:(n0-1))/n0)*Bartlett(c(1:(n0-1))/11)^2)/(n0*pi) + 1/(2*pi)
Dn_BAR_11 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Bartlett(c(1:(n0-2))/11)^4)*2/pi^2
###########

####TUK####
Cn_TUK_37 = sum((1 - c(1:(n0-1))/n0)*Tukey(c(1:(n0-1))/37)^2)/(n0*pi) + 1/(2*pi)
Dn_TUK_37 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Tukey(c(1:(n0-2))/37)^4)*2/pi^2

Cn_TUK_20 = sum((1 - c(1:(n0-1))/n0)*Tukey(c(1:(n0-1))/20)^2)/(n0*pi) + 1/(2*pi)
Dn_TUK_20 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Tukey(c(1:(n0-2))/20)^4)*2/pi^2

Cn_TUK_11 = sum((1 - c(1:(n0-1))/n0)*Tukey(c(1:(n0-1))/11)^2)/(n0*pi) + 1/(2*pi)
Dn_TUK_11 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Tukey(c(1:(n0-2))/11)^4)*2/pi^2
###########

#####QS####
Cn_QS_37 = sum((1 - c(1:(n0-1))/n0)*QS(c(1:(n0-1))/37)^2)/(n0*pi) + 1/(2*pi)
Dn_QS_37 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*QS(c(1:(n0-2))/37)^4)*2/pi^2

Cn_QS_20 = sum((1 - c(1:(n0-1))/n0)*QS(c(1:(n0-1))/20)^2)/(n0*pi) + 1/(2*pi)
Dn_QS_20 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*QS(c(1:(n0-2))/20)^4)*2/pi^2

Cn_QS_11 = sum((1 - c(1:(n0-1))/n0)*QS(c(1:(n0-1))/11)^2)/(n0*pi) + 1/(2*pi)
Dn_QS_11 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*QS(c(1:(n0-2))/11)^4)*2/pi^2
###########

#####################Li And Biao######################
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


for(sim in c(1:N)){
  
  ####### Table 1 & 3###### For t, use innov = rt(n0,9)
  x.sim.short = arima.sim(list(order = c(1,0,0),ar = .8),n = n0)#,innov = rt(n0,9))
  prd = Mod(fft(x.sim.short - mean(x.sim.short)))^2/(2*pi*n0)
  x.prd.short = prd[2:n0]
  x.prd.short.LB = prd[2:(n0/2)]

  x.sim.size = arima0(x.sim.short,order = c(1,0,0),method = "ML",include.mean = FALSE)
  phi.size = x.sim.size$coef
  sig2.size = x.sim.size$sigma2

  phi.est.v[sim] = phi.size
  sig2.est.v[sim] = sig2.size
  sp.short.est = sp.short.f.est(w.unit = mu.j,phi.1 = phi.size,sig2_f = 1)#sig2.size
  sp.short.est.LB = sp.short.f.est(w.unit = mu.j.LB,phi.1 = phi.size,sig2_f = 1)
  #######################
  
  ####### Table 2 & 4###### For t, use innov = rt(n0,9)
  # x.sim.long = (fracdiff.sim(n0,d = 0.4))$series
  # prd = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
  # x.prd.long = prd[2:n0]
  # 
  # x.MLE = arfima(x.sim.long,order = c(0,0,0),quiet = TRUE)$modes
  # d_l = x.MLE[[1]]$dfrac
  # sig2.est = x.MLE[[1]]$sigma2
  # d.est.v[sim] = d_l
  # sig2.est.v[sim] = sig2.est
  # sp.long.est = sp.long.f.est(w.loc = mu.j,d.loc = d_l,sig2_f = 1)#sig2.est
  #######################
  
  #####################Li And Biao######################
  ###Avg Periodogram###
  x.prd.avg = rep(0,sub.length)
  k_s = 0
  for (j_s in seq(1,m,3)) {
    k_s = k_s + 1
    x.prd.avg[k_s] = (x.prd.short.LB[j_s] + x.prd.short.LB[j_s + 1] + x.prd.short.LB[j_s + 2])/3
  }
  

  ###Filter Spectral Density###
  sp.filter.est = sp.short.est.LB[seq(2,m,3)]
  sp.filter.true = sp.short.true.LB[seq(2,m,3)]
  
  # Calculate beta.0.
  # This is for the case pars are unkown, please refer doc for details#######
  beta.0[sim] = beta.0.f(I = x.prd.short.LB,sp = sp.short.est.LB)
  beta.0.true[sim] = beta.0.f(I = x.prd.short.LB,sp = sp.short.true.LB)
  #1/3
  beta.0.Avg[sim] = beta.0.f(I = x.prd.avg,sp = sp.filter.est)
  beta.0.Avg.true[sim] = beta.0.f(I = x.prd.avg,sp = sp.filter.true)
  ###########################################################################
  
  #Calculate beta1, beta2, ...
  for (i in k0) {
    beta.k[i,sim] = beta.k.f(k = i,I = x.prd.short.LB, sp = sp.short.est.LB, mu.j.in = mu.j.LB)
    beta.k.true[i,sim] = beta.k.f(k = i,I = x.prd.short.LB, sp = sp.short.true.LB, mu.j.in = mu.j.LB)
    #1/3
    beta.k.Avg[i,sim] = beta.k.f(k = i,I = x.prd.avg, sp = sp.filter.est, mu.j.in = mu.j.LB.filter)
    beta.k.Avg.true[i,sim] = beta.k.f(k = i,I = x.prd.avg, sp = sp.filter.true, mu.j.in = mu.j.LB.filter)
  }
  
  #####################Chen And Deo#####################
  ###BAR###
  sp.disc.size.BAR_11 = (x.prd.short/sp.short.est)%*%w_BAR_11.fft.test*2*pi/n0
  sp.disc.BAR_11 = (x.prd.short/sp.short.true)%*%w_BAR_11.fft.test*2*pi/n0
  
  Tn.BAR_11[sim] = 2*pi/n0*sum(sp.disc.size.BAR_11^2)/(2*pi/n0*sum(sp.disc.size.BAR_11))^2
  Tn_q.BAR_11[sim] = 2*pi/n0*sum(sp.disc.BAR_11^2)/(2*pi/n0*sum(sp.disc.BAR_11))^2
  
  sp.disc.size.BAR_20 = (x.prd.short/sp.short.est)%*%w_BAR_20.fft.test*2*pi/n0
  sp.disc.BAR_20 = (x.prd.short/sp.short.true)%*%w_BAR_20.fft.test*2*pi/n0
  
  Tn.BAR_20[sim] = 2*pi/n0*sum(sp.disc.size.BAR_20^2)/(2*pi/n0*sum(sp.disc.size.BAR_20))^2
  Tn_q.BAR_20[sim] = 2*pi/n0*sum(sp.disc.BAR_20^2)/(2*pi/n0*sum(sp.disc.BAR_20))^2
  
  sp.disc.size.BAR_37 = (x.prd.short/sp.short.est)%*%w_BAR_37.fft.test*2*pi/n0
  sp.disc.BAR_37 = (x.prd.short/sp.short.true)%*%w_BAR_37.fft.test*2*pi/n0
  
  Tn.BAR_37[sim] = 2*pi/n0*sum(sp.disc.size.BAR_37^2)/(2*pi/n0*sum(sp.disc.size.BAR_37))^2
  Tn_q.BAR_37[sim] = 2*pi/n0*sum(sp.disc.BAR_37^2)/(2*pi/n0*sum(sp.disc.BAR_37))^2
  #########
  
  ###TUK###
  sp.disc.size.TUK_11 = (x.prd.short/sp.short.est)%*%w_TUK_11.fft.test*2*pi/n0
  sp.disc.TUK_11 = (x.prd.short/sp.short.true)%*%w_TUK_11.fft.test*2*pi/n0
  
  Tn.TUK_11[sim] = 2*pi/n0*sum(sp.disc.size.TUK_11^2)/(2*pi/n0*sum(sp.disc.size.TUK_11))^2
  Tn_q.TUK_11[sim] = 2*pi/n0*sum(sp.disc.TUK_11^2)/(2*pi/n0*sum(sp.disc.TUK_11))^2
  
  sp.disc.size.TUK_20 = (x.prd.short/sp.short.est)%*%w_TUK_20.fft.test*2*pi/n0
  sp.disc.TUK_20 = (x.prd.short/sp.short.true)%*%w_TUK_20.fft.test*2*pi/n0
  
  Tn.TUK_20[sim] = 2*pi/n0*sum(sp.disc.size.TUK_20^2)/(2*pi/n0*sum(sp.disc.size.TUK_20))^2
  Tn_q.TUK_20[sim] = 2*pi/n0*sum(sp.disc.TUK_20^2)/(2*pi/n0*sum(sp.disc.TUK_20))^2
  
  sp.disc.size.TUK_37 = (x.prd.short/sp.short.est)%*%w_TUK_37.fft.test*2*pi/n0
  sp.disc.TUK_37 = (x.prd.short/sp.short.true)%*%w_TUK_37.fft.test*2*pi/n0
  
  Tn.TUK_37[sim] = 2*pi/n0*sum(sp.disc.size.TUK_37^2)/(2*pi/n0*sum(sp.disc.size.TUK_37))^2
  Tn_q.TUK_37[sim] = 2*pi/n0*sum(sp.disc.TUK_37^2)/(2*pi/n0*sum(sp.disc.TUK_37))^2
  #########
  
  ###QS###
  sp.disc.size.QS_11 = (x.prd.short/sp.short.est)%*%w_QS_11.fft.test*2*pi/n0
  sp.disc.QS_11 = (x.prd.short/sp.short.true)%*%w_QS_11.fft.test*2*pi/n0
  
  Tn.QS_11[sim] = 2*pi/n0*sum(sp.disc.size.QS_11^2)/(2*pi/n0*sum(sp.disc.size.QS_11))^2
  Tn_q.QS_11[sim] = 2*pi/n0*sum(sp.disc.QS_11^2)/(2*pi/n0*sum(sp.disc.QS_11))^2
  
  sp.disc.size.QS_20 = (x.prd.short/sp.short.est)%*%w_QS_20.fft.test*2*pi/n0
  sp.disc.QS_20 = (x.prd.short/sp.short.true)%*%w_QS_20.fft.test*2*pi/n0
  
  Tn.QS_20[sim] = 2*pi/n0*sum(sp.disc.size.QS_20^2)/(2*pi/n0*sum(sp.disc.size.QS_20))^2
  Tn_q.QS_20[sim] = 2*pi/n0*sum(sp.disc.QS_20^2)/(2*pi/n0*sum(sp.disc.QS_20))^2
  
  sp.disc.size.QS_37 = (x.prd.short/sp.short.est)%*%w_QS_37.fft.test*2*pi/n0
  sp.disc.QS_37 = (x.prd.short/sp.short.true)%*%w_QS_37.fft.test*2*pi/n0
  
  Tn.QS_37[sim] = 2*pi/n0*sum(sp.disc.size.QS_37^2)/(2*pi/n0*sum(sp.disc.size.QS_37))^2
  Tn_q.QS_37[sim] = 2*pi/n0*sum(sp.disc.QS_37^2)/(2*pi/n0*sum(sp.disc.QS_37))^2
  #########
  
  progress(sim, max.value = N)
  if (sim == N) cat("Done!\n")
}


for (v0 in c(0:v)) {
  # alpha.1 = .05
  # alpha.2 = .1
  
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

#####################Chen And Deo#####################
###BAR###
Tn.standard.BAR_37 = n0*(Tn.BAR_37 - Cn_BAR_37)/sqrt(Dn_BAR_37)
Tn_q.standard.BAR_37 = n0*(Tn_q.BAR_37 - Cn_BAR_37)/sqrt(Dn_BAR_37)

Tn.standard.BAR_20 = n0*(Tn.BAR_20 - Cn_BAR_20)/sqrt(Dn_BAR_20)
Tn_q.standard.BAR_20 = n0*(Tn_q.BAR_20 - Cn_BAR_20)/sqrt(Dn_BAR_20)

Tn.standard.BAR_11 = n0*(Tn.BAR_11 - Cn_BAR_11)/sqrt(Dn_BAR_11)
Tn_q.standard.BAR_11 = n0*(Tn_q.BAR_11 - Cn_BAR_11)/sqrt(Dn_BAR_11)
#########

###TUK###
Tn.standard.TUK_37 = n0*(Tn.TUK_37 - Cn_TUK_37)/sqrt(Dn_TUK_37)
Tn_q.standard.TUK_37 = n0*(Tn_q.TUK_37 - Cn_TUK_37)/sqrt(Dn_TUK_37)

Tn.standard.TUK_20 = n0*(Tn.TUK_20 - Cn_TUK_20)/sqrt(Dn_TUK_20)
Tn_q.standard.TUK_20 = n0*(Tn_q.TUK_20 - Cn_TUK_20)/sqrt(Dn_TUK_20)

Tn.standard.TUK_11 = n0*(Tn.TUK_11 - Cn_TUK_11)/sqrt(Dn_TUK_11)
Tn_q.standard.TUK_11 = n0*(Tn_q.TUK_11 - Cn_TUK_11)/sqrt(Dn_TUK_11)
#########

###QS###
Tn.standard.QS_37 = n0*(Tn.QS_37 - Cn_QS_37)/sqrt(Dn_QS_37)
Tn_q.standard.QS_37 = n0*(Tn_q.QS_37 - Cn_QS_37)/sqrt(Dn_QS_37)

Tn.standard.QS_20 = n0*(Tn.QS_20 - Cn_QS_20)/sqrt(Dn_QS_20)
Tn_q.standard.QS_20 = n0*(Tn_q.QS_20 - Cn_QS_20)/sqrt(Dn_QS_20)

Tn.standard.QS_11 = n0*(Tn.QS_11 - Cn_QS_11)/sqrt(Dn_QS_11)
Tn_q.standard.QS_11 = n0*(Tn_q.QS_11 - Cn_QS_11)/sqrt(Dn_QS_11)
#########

#95th and 90th Quantile
###BAR###
sig_05.BAR_11 = quantile(Tn_q.standard.BAR_11, 1 - alpha.1, names = FALSE)
sig_10.BAR_11 = quantile(Tn_q.standard.BAR_11, 1 - alpha.2, names = FALSE)

r.rate.1.emp.BAR_11 = 100*length(Tn.standard.BAR_11[Tn.standard.BAR_11 > sig_05.BAR_11])/N
r.rate.2.emp.BAR_11 = 100*length(Tn.standard.BAR_11[Tn.standard.BAR_11 > sig_10.BAR_11])/N

sig_05.BAR_20 = quantile(Tn_q.standard.BAR_20, 1 - alpha.1, names = FALSE)
sig_10.BAR_20 = quantile(Tn_q.standard.BAR_20, 1 - alpha.2, names = FALSE)

r.rate.1.emp.BAR_20 = 100*length(Tn.standard.BAR_20[Tn.standard.BAR_20 > sig_05.BAR_20])/N
r.rate.2.emp.BAR_20 = 100*length(Tn.standard.BAR_20[Tn.standard.BAR_20 > sig_10.BAR_20])/N

sig_05.BAR_37 = quantile(Tn_q.standard.BAR_37, 1 - alpha.1, names = FALSE)
sig_10.BAR_37 = quantile(Tn_q.standard.BAR_37, 1 - alpha.2, names = FALSE)

r.rate.1.emp.BAR_37 = 100*length(Tn.standard.BAR_37[Tn.standard.BAR_37 > sig_05.BAR_37])/N
r.rate.2.emp.BAR_37 = 100*length(Tn.standard.BAR_37[Tn.standard.BAR_37 > sig_10.BAR_37])/N
#########

###TUK###
sig_05.TUK_11 = quantile(Tn_q.standard.TUK_11, 1 - alpha.1, names = FALSE)
sig_10.TUK_11 = quantile(Tn_q.standard.TUK_11, 1 - alpha.2, names = FALSE)

r.rate.1.emp.TUK_11 = 100*length(Tn.standard.TUK_11[Tn.standard.TUK_11 > sig_05.TUK_11])/N
r.rate.2.emp.TUK_11 = 100*length(Tn.standard.TUK_11[Tn.standard.TUK_11 > sig_10.TUK_11])/N

sig_05.TUK_20 = quantile(Tn_q.standard.TUK_20, 1 - alpha.1, names = FALSE)
sig_10.TUK_20 = quantile(Tn_q.standard.TUK_20, 1 - alpha.2, names = FALSE)

r.rate.1.emp.TUK_20 = 100*length(Tn.standard.TUK_20[Tn.standard.TUK_20 > sig_05.TUK_20])/N
r.rate.2.emp.TUK_20 = 100*length(Tn.standard.TUK_20[Tn.standard.TUK_20 > sig_10.TUK_20])/N

sig_05.TUK_37 = quantile(Tn_q.standard.TUK_37, 1 - alpha.1, names = FALSE)
sig_10.TUK_37 = quantile(Tn_q.standard.TUK_37, 1 - alpha.2, names = FALSE)

r.rate.1.emp.TUK_37 = 100*length(Tn.standard.TUK_37[Tn.standard.TUK_37 > sig_05.TUK_37])/N
r.rate.2.emp.TUK_37 = 100*length(Tn.standard.TUK_37[Tn.standard.TUK_37 > sig_10.TUK_37])/N
#########

###QS###
sig_05.QS_11 = quantile(Tn_q.standard.QS_11, 1 - alpha.1, names = FALSE)
sig_10.QS_11 = quantile(Tn_q.standard.QS_11, 1 - alpha.2, names = FALSE)

r.rate.1.emp.QS_11 = 100*length(Tn.standard.QS_11[Tn.standard.QS_11 > sig_05.QS_11])/N
r.rate.2.emp.QS_11 = 100*length(Tn.standard.QS_11[Tn.standard.QS_11 > sig_10.QS_11])/N

sig_05.QS_20 = quantile(Tn_q.standard.QS_20, 1 - alpha.1, names = FALSE)
sig_10.QS_20 = quantile(Tn_q.standard.QS_20, 1 - alpha.2, names = FALSE)

r.rate.1.emp.QS_20 = 100*length(Tn.standard.QS_20[Tn.standard.QS_20 > sig_05.QS_20])/N
r.rate.2.emp.QS_20 = 100*length(Tn.standard.QS_20[Tn.standard.QS_20 > sig_10.QS_20])/N

sig_05.QS_37 = quantile(Tn_q.standard.QS_37, 1 - alpha.1, names = FALSE)
sig_10.QS_37 = quantile(Tn_q.standard.QS_37, 1 - alpha.2, names = FALSE)

r.rate.1.emp.QS_37 = 100*length(Tn.standard.QS_37[Tn.standard.QS_37 > sig_05.QS_37])/N
r.rate.2.emp.QS_37 = 100*length(Tn.standard.QS_37[Tn.standard.QS_37 > sig_10.QS_37])/N
#########

#95th and 90th Nromal Quantile
sig_05.true = qnorm(alpha.1, lower.tail = FALSE)
sig_10.true = qnorm(alpha.2, lower.tail = FALSE)

###BAR###
r.rate.1.BAR_11 = 100*length(Tn.standard.BAR_11[Tn.standard.BAR_11 > sig_05.true])/N
r.rate.2.BAR_11 = 100*length(Tn.standard.BAR_11[Tn.standard.BAR_11 > sig_10.true])/N

r.rate.1.BAR_20 = 100*length(Tn.standard.BAR_20[Tn.standard.BAR_20 > sig_05.true])/N
r.rate.2.BAR_20 = 100*length(Tn.standard.BAR_20[Tn.standard.BAR_20 > sig_10.true])/N

r.rate.1.BAR_37 = 100*length(Tn.standard.BAR_37[Tn.standard.BAR_37 > sig_05.true])/N
r.rate.2.BAR_37 = 100*length(Tn.standard.BAR_37[Tn.standard.BAR_37 > sig_10.true])/N
#########

###TUK###
r.rate.1.TUK_11 = 100*length(Tn.standard.TUK_11[Tn.standard.TUK_11 > sig_05.true])/N
r.rate.2.TUK_11 = 100*length(Tn.standard.TUK_11[Tn.standard.TUK_11 > sig_10.true])/N

r.rate.1.TUK_20 = 100*length(Tn.standard.TUK_20[Tn.standard.TUK_20 > sig_05.true])/N
r.rate.2.TUK_20 = 100*length(Tn.standard.TUK_20[Tn.standard.TUK_20 > sig_10.true])/N

r.rate.1.TUK_37 = 100*length(Tn.standard.TUK_37[Tn.standard.TUK_37 > sig_05.true])/N
r.rate.2.TUK_37 = 100*length(Tn.standard.TUK_37[Tn.standard.TUK_37 > sig_10.true])/N
#########

###QS###
r.rate.1.QS_11 = 100*length(Tn.standard.QS_11[Tn.standard.QS_11 > sig_05.true])/N
r.rate.2.QS_11 = 100*length(Tn.standard.QS_11[Tn.standard.QS_11 > sig_10.true])/N

r.rate.1.QS_20 = 100*length(Tn.standard.QS_20[Tn.standard.QS_20 > sig_05.true])/N
r.rate.2.QS_20 = 100*length(Tn.standard.QS_20[Tn.standard.QS_20 > sig_10.true])/N

r.rate.1.QS_37 = 100*length(Tn.standard.QS_37[Tn.standard.QS_37 > sig_05.true])/N
r.rate.2.QS_37 = 100*length(Tn.standard.QS_37[Tn.standard.QS_37 > sig_10.true])/N
#########

size.emp.table = data.frame(
  B11_5perct = c(r.rate.1.emp.BAR_11, r.rate.1.emp.TUK_11, r.rate.1.emp.QS_11),
  B11_10perct = c(r.rate.2.emp.BAR_11, r.rate.2.emp.TUK_11, r.rate.2.emp.QS_11),
  B20_5perct = c(r.rate.1.emp.BAR_20, r.rate.1.emp.TUK_20, r.rate.1.emp.QS_20),
  B20_10perct = c(r.rate.2.emp.BAR_20, r.rate.2.emp.TUK_20, r.rate.2.emp.QS_20),
  B37_5perct = c(r.rate.1.emp.BAR_37, r.rate.1.emp.TUK_37, r.rate.1.emp.QS_37),
  B37_10perct = c(r.rate.2.emp.BAR_37, r.rate.2.emp.TUK_37, r.rate.2.emp.QS_37)
)

size.table = data.frame(
  B11_5perct = c(r.rate.1.BAR_11, r.rate.1.TUK_11, r.rate.1.QS_11),
  B11_10perct = c(r.rate.2.BAR_11, r.rate.2.TUK_11, r.rate.2.QS_11),
  B20_5perct = c(r.rate.1.BAR_20, r.rate.1.TUK_20, r.rate.1.QS_20),
  B20_10perct = c(r.rate.2.BAR_20, r.rate.2.TUK_20, r.rate.2.QS_20),
  B37_5perct = c(r.rate.1.BAR_37, r.rate.1.TUK_37, r.rate.1.QS_37),
  B37_10perct = c(r.rate.2.BAR_37, r.rate.2.TUK_37, r.rate.2.QS_37)
)

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