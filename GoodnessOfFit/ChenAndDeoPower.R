library(fracdiff)
library(astsa)
library(arfima)
require(svMisc)

n0 = 1024
j = c(1:(n0-1))
mu.j = 2*pi*j/n0


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
W_BAR_48.fft = function(j_f){
  1/(2*pi)*(fft(Bartlett(c(0:(n0-1))/48)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Bartlett(c(0:(n0-1))/48)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Bartlett(0))
}
w_BAR_48.fft.test = sapply(j, W_BAR_48.fft)
w_BAR_48.fft.test = t(Re(w_BAR_48.fft.test))

W_BAR_24.fft = function(j_f){
  1/(2*pi)*(fft(Bartlett(c(0:(n0-1))/24)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Bartlett(c(0:(n0-1))/24)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Bartlett(0))
}
w_BAR_24.fft.test = sapply(j, W_BAR_24.fft)
w_BAR_24.fft.test = t(Re(w_BAR_24.fft.test))

W_BAR_12.fft = function(j_f){
  1/(2*pi)*(fft(Bartlett(c(0:(n0-1))/12)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Bartlett(c(0:(n0-1))/12)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Bartlett(0))
}
w_BAR_12.fft.test = sapply(j, W_BAR_12.fft)
w_BAR_12.fft.test = t(Re(w_BAR_12.fft.test))
##################


######Tukey#######
W_TUK_48.fft = function(j_f){
  1/(2*pi)*(fft(Tukey(c(0:(n0-1))/48)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Tukey(c(0:(n0-1))/48)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Tukey(0))
}
w_TUK_48.fft.test = sapply(j, W_TUK_48.fft)
w_TUK_48.fft.test = t(Re(w_TUK_48.fft.test))

W_TUK_24.fft = function(j_f){
  1/(2*pi)*(fft(Tukey(c(0:(n0-1))/24)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Tukey(c(0:(n0-1))/24)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Tukey(0))
}
w_TUK_24.fft.test = sapply(j, W_TUK_24.fft)
w_TUK_24.fft.test = t(Re(w_TUK_24.fft.test))

W_TUK_12.fft = function(j_f){
  1/(2*pi)*(fft(Tukey(c(0:(n0-1))/12)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(Tukey(c(0:(n0-1))/12)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -Tukey(0))
}
w_TUK_12.fft.test = sapply(j, W_TUK_12.fft)
w_TUK_12.fft.test = t(Re(w_TUK_12.fft.test))
##################

########QS########
W_QS_48.fft = function(j_f){
  1/(2*pi)*(fft(QS(c(0:(n0-1))/48)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(QS(c(0:(n0-1))/48)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -QS(0))
}
w_QS_48.fft.test = sapply(j, W_QS_48.fft)
w_QS_48.fft.test = t(Re(w_QS_48.fft.test))

W_QS_24.fft = function(j_f){
  1/(2*pi)*(fft(QS(c(0:(n0-1))/24)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(QS(c(0:(n0-1))/24)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -QS(0))
}
w_QS_24.fft.test = sapply(j, W_QS_24.fft)
w_QS_24.fft.test = t(Re(w_QS_24.fft.test))

W_QS_12.fft = function(j_f){
  1/(2*pi)*(fft(QS(c(0:(n0-1))/12)*exp(complex(imaginary = 2*pi*c(0:(n0-1))*j_f/n0)))
            +fft(QS(c(0:(n0-1))/12)*exp(complex(imaginary = -2*pi*c(0:(n0-1))*j_f/n0)),inverse = TRUE)
            -QS(0))
}
w_QS_12.fft.test = sapply(j, W_QS_12.fft)
w_QS_12.fft.test = t(Re(w_QS_12.fft.test))
##################

####Table 5####
sp.ar1.f.est = function(w.unit,phi.1,sig2_f){(sig2_f/((1 + phi.1^2) - 2*phi.1*cos(w.unit)))/(2*pi)}

####Table 6####
# sp.long.f.est = function(w.loc,phi.loc,d.loc,sig2_f){sig2_f/(2*pi)*(4*(sin(w.loc/2))^2)^(-d.loc)/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}

####Table 7####
# sp.ar2.f.est = function(w.loc,phi.loc,theta.loc,sig2_f){sig2_f/(2*pi)*(1 + theta.loc^2 + 2*theta.loc*cos(w.loc))/(1 + phi.loc^2 - 2*phi.loc*cos(w.loc))}

####Table 8####
# sp.long.f.est = function(w.loc,d.loc,sig2_f){sig2_f*(4*(sin(w.loc/2))^2)^(-d.loc)/(2*pi)}

N = 5000
###BAR###
Tn.BAR_48 = rep(0,N)
Tn_q.BAR_48 = rep(0,N)

Tn.BAR_24 = rep(0,N)
Tn_q.BAR_24 = rep(0,N)

Tn.BAR_12 = rep(0,N)
Tn_q.BAR_12 = rep(0,N)
#########

###TUK###
Tn.TUK_48 = rep(0,N)
Tn_q.TUK_48 = rep(0,N)

Tn.TUK_24 = rep(0,N)
Tn_q.TUK_24 = rep(0,N)

Tn.TUK_12 = rep(0,N)
Tn_q.TUK_12 = rep(0,N)
#########

###QS###
Tn.QS_48 = rep(0,N)
Tn_q.QS_48 = rep(0,N)

Tn.QS_24 = rep(0,N)
Tn_q.QS_24 = rep(0,N)

Tn.QS_12 = rep(0,N)
Tn_q.QS_12 = rep(0,N)
#########
d.est.v.pow = rep(0,N)
d.est.v = rep(0,N)

phi.est.v.pow = rep(0,N)
phi.est.v = rep(0,N)

theta.est.v.pow = rep(0,N)
theta.est.v = rep(0,N)

sig2.est.v.pow = rep(0,N)
sig2.est.v = rep(0,N)

alpha.1 = .05
alpha.2 = .1

####BAR####
Cn_BAR_48 = sum((1 - c(1:(n0-1))/n0)*Bartlett(c(1:(n0-1))/48)^2)/(n0*pi) + 1/(2*pi)
Dn_BAR_48 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Bartlett(c(1:(n0-2))/48)^4)*2/pi^2

Cn_BAR_24 = sum((1 - c(1:(n0-1))/n0)*Bartlett(c(1:(n0-1))/24)^2)/(n0*pi) + 1/(2*pi)
Dn_BAR_24 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Bartlett(c(1:(n0-2))/24)^4)*2/pi^2

Cn_BAR_12 = sum((1 - c(1:(n0-1))/n0)*Bartlett(c(1:(n0-1))/12)^2)/(n0*pi) + 1/(2*pi)
Dn_BAR_12 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Bartlett(c(1:(n0-2))/12)^4)*2/pi^2
###########

####TUK####
Cn_TUK_48 = sum((1 - c(1:(n0-1))/n0)*Tukey(c(1:(n0-1))/48)^2)/(n0*pi) + 1/(2*pi)
Dn_TUK_48 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Tukey(c(1:(n0-2))/48)^4)*2/pi^2

Cn_TUK_24 = sum((1 - c(1:(n0-1))/n0)*Tukey(c(1:(n0-1))/24)^2)/(n0*pi) + 1/(2*pi)
Dn_TUK_24 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Tukey(c(1:(n0-2))/24)^4)*2/pi^2

Cn_TUK_12 = sum((1 - c(1:(n0-1))/n0)*Tukey(c(1:(n0-1))/12)^2)/(n0*pi) + 1/(2*pi)
Dn_TUK_12 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*Tukey(c(1:(n0-2))/12)^4)*2/pi^2
###########

#####QS####
Cn_QS_48 = sum((1 - c(1:(n0-1))/n0)*QS(c(1:(n0-1))/48)^2)/(n0*pi) + 1/(2*pi)
Dn_QS_48 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*QS(c(1:(n0-2))/48)^4)*2/pi^2

Cn_QS_24 = sum((1 - c(1:(n0-1))/n0)*QS(c(1:(n0-1))/24)^2)/(n0*pi) + 1/(2*pi)
Dn_QS_24 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*QS(c(1:(n0-2))/24)^4)*2/pi^2

Cn_QS_12 = sum((1 - c(1:(n0-1))/n0)*QS(c(1:(n0-1))/12)^2)/(n0*pi) + 1/(2*pi)
Dn_QS_12 = sum((1 - c(1:(n0-2))/n0)*(1 - (c(1:(n0-2))+1)/n0)*QS(c(1:(n0-2))/12)^4)*2/pi^2
###########
# for(re in c(1:3)){

for(i in c(1:N)){
  #######Table 5######
  # simulation of AR(2)
  x.sim.ar1.pow = arima.sim(list(order = c(2,0,0),ar = c(.8,-.15)),n = n0)
  prd.pow = Mod(fft(x.sim.ar1.pow- mean(x.sim.ar1.pow)))^2/(2*pi*n0)
  x.prd.ar1.pow = prd.pow[2:n0]
  # simulation of AR(1)
  x.sim.ar1 = arima.sim(list(order = c(1,0,0),ar = .8),n = n0)
  prd = Mod(fft(x.sim.ar1 - mean(x.sim.ar1)))^2/(2*pi*n0)
  x.prd.ar1 = prd[2:n0]

  ##MLE##
  x.sim.pow = arima0(x.sim.ar1.pow,order = c(1,0,0),method = "ML",include.mean = FALSE)
  phi.whit.pow = x.sim.pow$coef
  sig2.whit.pow = x.sim.pow$sigma2

  phi.est.v.pow[i] = phi.whit.pow
  sig2.est.v.pow[i] = sig2.whit.pow
  sp.ar1.est.pow = sp.ar1.f.est(w.unit = mu.j,phi.1 = phi.whit.pow,sig2_f = 1)#sig2.whit.pow

  x.sim = arima0(x.sim.ar1,order = c(1,0,0),method = "ML",include.mean = FALSE)
  phi.whit = x.sim$coef
  sig2.whit = x.sim$sigma2

  phi.est.v[i] = phi.whit
  sig2.est.v[i] = sig2.whit
  sp.ar1.est = sp.ar1.f.est(w.unit = mu.j,phi.1 = phi.whit,sig2_f = 1)#sig2.whit
  #######
  ####################

  #######Table 6######
  # # simulation of ARMA(1,1)
  # x.sim.arma.pow = arima.sim(list(order = c(1,0,1),ar = .8,ma = .2),n = n0)
  # prd.arma.pow = Mod(fft(x.sim.arma.pow- mean(x.sim.arma.pow)))^2/(2*pi*n0)
  # x.prd.arma.pow = prd.arma.pow[2:n0]
  # # simulation of ARFIMA(1,d,0)
  # x.sim.long = (fracdiff.sim(n0,ar = .1,d = .4))$series
  # prd.long = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
  # x.prd.long = prd.long[2:n0]
  # 
  # ##MLE##
  # x.sim.pow = arfima(x.sim.arma.pow,order = c(1,0,0),quiet = TRUE)$modes
  # d.whit.pow = x.sim.pow[[1]]$dfrac
  # phi.whit.pow = x.sim.pow[[1]]$phi
  # sig2.whit.pow = x.sim.pow[[1]]$sigma2
  # 
  # d.est.v.pow[i] = d.whit.pow
  # phi.est.v.pow[i] = phi.whit.pow
  # sig2.est.v.pow[i] = sig2.whit.pow
  # sp.arma.est = sp.long.f.est(w.loc = mu.j,phi.loc = phi.whit.pow,d.loc = d.whit.pow,sig2_f = 1)#sig2.whit.pow
  # 
  # x.sim = arfima(x.sim.long,order = c(1,0,0),quiet = TRUE)$modes
  # d.whit = x.sim[[1]]$dfrac
  # phi.whit = x.sim[[1]]$phi
  # sig2.whit = x.sim[[1]]$sigma2
  # 
  # d.est.v[i] = d.whit
  # phi.est.v[i] = phi.whit
  # sig2.est.v[i] = sig2.whit
  # sp.long.est = sp.long.f.est(w.loc = mu.j,phi.loc = phi.whit,d.loc = d.whit,sig2_f = 1)#sig2.whit
  # #######
  ####################
  
  #######Table 7######
  # # simulation of ARFIMA(d)
  # x.sim.long = (fracdiff.sim(n0,d = 0.4))$series
  # prd = Mod(fft(x.sim.long- mean(x.sim.long)))^2/(2*pi*n0)
  # x.prd.long = prd[2:n0]
  # # simulation of ARMA(1,0,1)
  # x.sim.arma = arima.sim(list(order = c(1,0,1),ar = .8,ma = .2),n = n0)
  # prd.arma = Mod(fft(x.sim.arma- mean(x.sim.arma)))^2/(2*pi*n0)
  # x.prd.arma = prd.arma[2:n0]
  # 
  # ##MLE##
  # x.sim.pow = arima0(x.sim.long,order = c(1,0,1),method = "ML",include.mean = FALSE)
  # phi.whit.pow = x.sim.pow$coef[1]
  # theta.whit.pow = x.sim.pow$coef[2]
  # sig2.whit.pow = x.sim.pow$sigma2
  # 
  # phi.est.v.pow[i] = phi.whit.pow
  # theta.est.v.pow[i] = theta.whit.pow
  # sig2.est.v.pow[i] = sig2.whit.pow
  # sp.ar2.est.pow = sp.ar2.f.est(w.loc = mu.j,phi.loc = phi.whit.pow,theta.loc = theta.whit.pow,sig2_f = 1)#sig2.whit.pow
  # 
  # x.sim = arima0(x.sim.arma,order = c(1,0,1),method = "ML",include.mean = FALSE)
  # phi.whit = x.sim$coef[1]
  # theta.whit = x.sim$coef[2]
  # sig2.whit = x.sim$sigma2
  # 
  # phi.est.v[i] = phi.whit
  # theta.est.v[i] = theta.whit
  # sig2.est.v[i] = sig2.whit
  # sp.ar2.est = sp.ar2.f.est(w.loc = mu.j,phi.loc = phi.whit,theta.loc = theta.whit,sig2_f = 1)#sig2.whit
  # #######
  ####################
  
    
  #######Table 8######
  # # simulation of ARFIMA(1,d,0)
  # x.sim.long.pow = (fracdiff.sim(n0,ar = .1,d = .4))$series
  # prd.long.pow = Mod(fft(x.sim.long.pow- mean(x.sim.long.pow)))^2/(2*pi*n0)
  # x.prd.long.pow = prd.long.pow[2:n0]
  # # simulation of ARFIMA(d)
  # x.sim.long.bm = (fracdiff.sim(n0,d = .4))$series
  # prd.long.bm = Mod(fft(x.sim.long.bm- mean(x.sim.long.bm)))^2/(2*pi*n0)
  # x.prd.long.bm = prd.long.bm[2:n0]
  # 
  # ##MLE##
  # x.sim.pow = arfima(x.sim.long.pow,order = c(0,0,0),quiet = TRUE)$modes
  # d.whit.pow = x.sim.pow[[1]]$dfrac
  # sig2.whit.pow = x.sim.pow[[1]]$sigma2
  # d.est.v.pow[i] = d.whit.pow
  # sig2.est.v.pow[i] = sig2.whit.pow
  # sp.long.est.pow = sp.long.f.est(w.loc = mu.j,d.loc = d.whit.pow,sig2_f = 1)#sig2.whit.pow
  # 
  # x.sim.bm = arfima(x.sim.long.bm,order = c(0,0,0),quiet = TRUE)$modes
  # d.whit.bm = x.sim.bm[[1]]$dfrac
  # sig2.whit.bm = x.sim.bm[[1]]$sigma2
  # d.est.v[i] = d.whit.bm
  # sig2.est.v[i] = sig2.whit.bm
  # sp.long.est.bm = sp.long.f.est(w.loc = mu.j,d.loc = d.whit.bm,sig2_f = 1)#sig2.whit.bm
  # ########
  #####################
  
  ###BAR###
  sp.disc.pow.BAR_12 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_BAR_12.fft.test*2*pi/n0
  sp.disc.BAR_12 = (x.prd.ar1/sp.ar1.est)%*%w_BAR_12.fft.test*2*pi/n0
  
  Tn.BAR_12[i] = 2*pi/n0*sum(sp.disc.pow.BAR_12^2)/(2*pi/n0*sum(sp.disc.pow.BAR_12))^2
  Tn_q.BAR_12[i] = 2*pi/n0*sum(sp.disc.BAR_12^2)/(2*pi/n0*sum(sp.disc.BAR_12))^2
  
  sp.disc.pow.BAR_24 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_BAR_24.fft.test*2*pi/n0
  sp.disc.BAR_24 = (x.prd.ar1/sp.ar1.est)%*%w_BAR_24.fft.test*2*pi/n0
  
  Tn.BAR_24[i] = 2*pi/n0*sum(sp.disc.pow.BAR_24^2)/(2*pi/n0*sum(sp.disc.pow.BAR_24))^2
  Tn_q.BAR_24[i] = 2*pi/n0*sum(sp.disc.BAR_24^2)/(2*pi/n0*sum(sp.disc.BAR_24))^2
  
  sp.disc.pow.BAR_48 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_BAR_48.fft.test*2*pi/n0
  sp.disc.BAR_48 = (x.prd.ar1/sp.ar1.est)%*%w_BAR_48.fft.test*2*pi/n0
  
  Tn.BAR_48[i] = 2*pi/n0*sum(sp.disc.pow.BAR_48^2)/(2*pi/n0*sum(sp.disc.pow.BAR_48))^2
  Tn_q.BAR_48[i] = 2*pi/n0*sum(sp.disc.BAR_48^2)/(2*pi/n0*sum(sp.disc.BAR_48))^2
  #########
  
  ###TUK###
  sp.disc.pow.TUK_12 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_TUK_12.fft.test*2*pi/n0
  sp.disc.TUK_12 = (x.prd.ar1/sp.ar1.est)%*%w_TUK_12.fft.test*2*pi/n0
  
  Tn.TUK_12[i] = 2*pi/n0*sum(sp.disc.pow.TUK_12^2)/(2*pi/n0*sum(sp.disc.pow.TUK_12))^2
  Tn_q.TUK_12[i] = 2*pi/n0*sum(sp.disc.TUK_12^2)/(2*pi/n0*sum(sp.disc.TUK_12))^2
  
  sp.disc.pow.TUK_24 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_TUK_24.fft.test*2*pi/n0
  sp.disc.TUK_24 = (x.prd.ar1/sp.ar1.est)%*%w_TUK_24.fft.test*2*pi/n0
  
  Tn.TUK_24[i] = 2*pi/n0*sum(sp.disc.pow.TUK_24^2)/(2*pi/n0*sum(sp.disc.pow.TUK_24))^2
  Tn_q.TUK_24[i] = 2*pi/n0*sum(sp.disc.TUK_24^2)/(2*pi/n0*sum(sp.disc.TUK_24))^2
  
  sp.disc.pow.TUK_48 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_TUK_48.fft.test*2*pi/n0
  sp.disc.TUK_48 = (x.prd.ar1/sp.ar1.est)%*%w_TUK_48.fft.test*2*pi/n0
  
  Tn.TUK_48[i] = 2*pi/n0*sum(sp.disc.pow.TUK_48^2)/(2*pi/n0*sum(sp.disc.pow.TUK_48))^2
  Tn_q.TUK_48[i] = 2*pi/n0*sum(sp.disc.TUK_48^2)/(2*pi/n0*sum(sp.disc.TUK_48))^2
  #########
  
  ###QS###
  sp.disc.pow.QS_12 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_QS_12.fft.test*2*pi/n0
  sp.disc.QS_12 = (x.prd.ar1/sp.ar1.est)%*%w_QS_12.fft.test*2*pi/n0
  
  Tn.QS_12[i] = 2*pi/n0*sum(sp.disc.pow.QS_12^2)/(2*pi/n0*sum(sp.disc.pow.QS_12))^2
  Tn_q.QS_12[i] = 2*pi/n0*sum(sp.disc.QS_12^2)/(2*pi/n0*sum(sp.disc.QS_12))^2
  
  sp.disc.pow.QS_24 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_QS_24.fft.test*2*pi/n0
  sp.disc.QS_24 = (x.prd.ar1/sp.ar1.est)%*%w_QS_24.fft.test*2*pi/n0
  
  Tn.QS_24[i] = 2*pi/n0*sum(sp.disc.pow.QS_24^2)/(2*pi/n0*sum(sp.disc.pow.QS_24))^2
  Tn_q.QS_24[i] = 2*pi/n0*sum(sp.disc.QS_24^2)/(2*pi/n0*sum(sp.disc.QS_24))^2
  
  sp.disc.pow.QS_48 = (x.prd.ar1.pow/sp.ar1.est.pow)%*%w_QS_48.fft.test*2*pi/n0
  sp.disc.QS_48 = (x.prd.ar1/sp.ar1.est)%*%w_QS_48.fft.test*2*pi/n0
  
  Tn.QS_48[i] = 2*pi/n0*sum(sp.disc.pow.QS_48^2)/(2*pi/n0*sum(sp.disc.pow.QS_48))^2
  Tn_q.QS_48[i] = 2*pi/n0*sum(sp.disc.QS_48^2)/(2*pi/n0*sum(sp.disc.QS_48))^2
  #########
  progress(i, max.value = N)
  if (i == N) cat("Done!\n")
}

###BAR###
Tn.standard.BAR_48 = n0*(Tn.BAR_48 - Cn_BAR_48)/sqrt(Dn_BAR_48)
Tn_q.standard.BAR_48 = n0*(Tn_q.BAR_48 - Cn_BAR_48)/sqrt(Dn_BAR_48)

Tn.standard.BAR_24 = n0*(Tn.BAR_24 - Cn_BAR_24)/sqrt(Dn_BAR_24)
Tn_q.standard.BAR_24 = n0*(Tn_q.BAR_24 - Cn_BAR_24)/sqrt(Dn_BAR_24)

Tn.standard.BAR_12 = n0*(Tn.BAR_12 - Cn_BAR_12)/sqrt(Dn_BAR_12)
Tn_q.standard.BAR_12 = n0*(Tn_q.BAR_12 - Cn_BAR_12)/sqrt(Dn_BAR_12)
#########

###TUK###
Tn.standard.TUK_48 = n0*(Tn.TUK_48 - Cn_TUK_48)/sqrt(Dn_TUK_48)
Tn_q.standard.TUK_48 = n0*(Tn_q.TUK_48 - Cn_TUK_48)/sqrt(Dn_TUK_48)

Tn.standard.TUK_24 = n0*(Tn.TUK_24 - Cn_TUK_24)/sqrt(Dn_TUK_24)
Tn_q.standard.TUK_24 = n0*(Tn_q.TUK_24 - Cn_TUK_24)/sqrt(Dn_TUK_24)

Tn.standard.TUK_12 = n0*(Tn.TUK_12 - Cn_TUK_12)/sqrt(Dn_TUK_12)
Tn_q.standard.TUK_12 = n0*(Tn_q.TUK_12 - Cn_TUK_12)/sqrt(Dn_TUK_12)
#########

###QS###
Tn.standard.QS_48 = n0*(Tn.QS_48 - Cn_QS_48)/sqrt(Dn_QS_48)
Tn_q.standard.QS_48 = n0*(Tn_q.QS_48 - Cn_QS_48)/sqrt(Dn_QS_48)

Tn.standard.QS_24 = n0*(Tn.QS_24 - Cn_QS_24)/sqrt(Dn_QS_24)
Tn_q.standard.QS_24 = n0*(Tn_q.QS_24 - Cn_QS_24)/sqrt(Dn_QS_24)

Tn.standard.QS_12 = n0*(Tn.QS_12 - Cn_QS_12)/sqrt(Dn_QS_12)
Tn_q.standard.QS_12 = n0*(Tn_q.QS_12 - Cn_QS_12)/sqrt(Dn_QS_12)
#########

#95th and 90th Quantile
###BAR###
sig_05.BAR_12 = quantile(Tn_q.standard.BAR_12, 1 - alpha.1, names = FALSE)
sig_10.BAR_12 = quantile(Tn_q.standard.BAR_12, 1 - alpha.2, names = FALSE)

r.rate.1.emp.BAR_12 = 100*length(Tn.standard.BAR_12[Tn.standard.BAR_12 > sig_05.BAR_12])/N
r.rate.2.emp.BAR_12 = 100*length(Tn.standard.BAR_12[Tn.standard.BAR_12 > sig_10.BAR_12])/N

sig_05.BAR_24 = quantile(Tn_q.standard.BAR_24, 1 - alpha.1, names = FALSE)
sig_10.BAR_24 = quantile(Tn_q.standard.BAR_24, 1 - alpha.2, names = FALSE)

r.rate.1.emp.BAR_24 = 100*length(Tn.standard.BAR_24[Tn.standard.BAR_24 > sig_05.BAR_24])/N
r.rate.2.emp.BAR_24 = 100*length(Tn.standard.BAR_24[Tn.standard.BAR_24 > sig_10.BAR_24])/N

sig_05.BAR_48 = quantile(Tn_q.standard.BAR_48, 1 - alpha.1, names = FALSE)
sig_10.BAR_48 = quantile(Tn_q.standard.BAR_48, 1 - alpha.2, names = FALSE)

r.rate.1.emp.BAR_48 = 100*length(Tn.standard.BAR_48[Tn.standard.BAR_48 > sig_05.BAR_48])/N
r.rate.2.emp.BAR_48 = 100*length(Tn.standard.BAR_48[Tn.standard.BAR_48 > sig_10.BAR_48])/N
#########

###TUK###
sig_05.TUK_12 = quantile(Tn_q.standard.TUK_12, 1 - alpha.1, names = FALSE)
sig_10.TUK_12 = quantile(Tn_q.standard.TUK_12, 1 - alpha.2, names = FALSE)

r.rate.1.emp.TUK_12 = 100*length(Tn.standard.TUK_12[Tn.standard.TUK_12 > sig_05.TUK_12])/N
r.rate.2.emp.TUK_12 = 100*length(Tn.standard.TUK_12[Tn.standard.TUK_12 > sig_10.TUK_12])/N

sig_05.TUK_24 = quantile(Tn_q.standard.TUK_24, 1 - alpha.1, names = FALSE)
sig_10.TUK_24 = quantile(Tn_q.standard.TUK_24, 1 - alpha.2, names = FALSE)

r.rate.1.emp.TUK_24 = 100*length(Tn.standard.TUK_24[Tn.standard.TUK_24 > sig_05.TUK_24])/N
r.rate.2.emp.TUK_24 = 100*length(Tn.standard.TUK_24[Tn.standard.TUK_24 > sig_10.TUK_24])/N

sig_05.TUK_48 = quantile(Tn_q.standard.TUK_48, 1 - alpha.1, names = FALSE)
sig_10.TUK_48 = quantile(Tn_q.standard.TUK_48, 1 - alpha.2, names = FALSE)

r.rate.1.emp.TUK_48 = 100*length(Tn.standard.TUK_48[Tn.standard.TUK_48 > sig_05.TUK_48])/N
r.rate.2.emp.TUK_48 = 100*length(Tn.standard.TUK_48[Tn.standard.TUK_48 > sig_10.TUK_48])/N
#########

###QS###
sig_05.QS_12 = quantile(Tn_q.standard.QS_12, 1 - alpha.1, names = FALSE)
sig_10.QS_12 = quantile(Tn_q.standard.QS_12, 1 - alpha.2, names = FALSE)

r.rate.1.emp.QS_12 = 100*length(Tn.standard.QS_12[Tn.standard.QS_12 > sig_05.QS_12])/N
r.rate.2.emp.QS_12 = 100*length(Tn.standard.QS_12[Tn.standard.QS_12 > sig_10.QS_12])/N

sig_05.QS_24 = quantile(Tn_q.standard.QS_24, 1 - alpha.1, names = FALSE)
sig_10.QS_24 = quantile(Tn_q.standard.QS_24, 1 - alpha.2, names = FALSE)

r.rate.1.emp.QS_24 = 100*length(Tn.standard.QS_24[Tn.standard.QS_24 > sig_05.QS_24])/N
r.rate.2.emp.QS_24 = 100*length(Tn.standard.QS_24[Tn.standard.QS_24 > sig_10.QS_24])/N

sig_05.QS_48 = quantile(Tn_q.standard.QS_48, 1 - alpha.1, names = FALSE)
sig_10.QS_48 = quantile(Tn_q.standard.QS_48, 1 - alpha.2, names = FALSE)

r.rate.1.emp.QS_48 = 100*length(Tn.standard.QS_48[Tn.standard.QS_48 > sig_05.QS_48])/N
r.rate.2.emp.QS_48 = 100*length(Tn.standard.QS_48[Tn.standard.QS_48 > sig_10.QS_48])/N
#########

#95th and 90th Nromal Quantile
sig_05.true = qnorm(alpha.1, lower.tail = FALSE)
sig_10.true = qnorm(alpha.2, lower.tail = FALSE)

###BAR###
r.rate.1.BAR_12 = 100*length(Tn.standard.BAR_12[Tn.standard.BAR_12 > sig_05.true])/N
r.rate.2.BAR_12 = 100*length(Tn.standard.BAR_12[Tn.standard.BAR_12 > sig_10.true])/N

r.rate.1.BAR_24 = 100*length(Tn.standard.BAR_24[Tn.standard.BAR_24 > sig_05.true])/N
r.rate.2.BAR_24 = 100*length(Tn.standard.BAR_24[Tn.standard.BAR_24 > sig_10.true])/N

r.rate.1.BAR_48 = 100*length(Tn.standard.BAR_48[Tn.standard.BAR_48 > sig_05.true])/N
r.rate.2.BAR_48 = 100*length(Tn.standard.BAR_48[Tn.standard.BAR_48 > sig_10.true])/N
#########

###TUK###
r.rate.1.TUK_12 = 100*length(Tn.standard.TUK_12[Tn.standard.TUK_12 > sig_05.true])/N
r.rate.2.TUK_12 = 100*length(Tn.standard.TUK_12[Tn.standard.TUK_12 > sig_10.true])/N

r.rate.1.TUK_24 = 100*length(Tn.standard.TUK_24[Tn.standard.TUK_24 > sig_05.true])/N
r.rate.2.TUK_24 = 100*length(Tn.standard.TUK_24[Tn.standard.TUK_24 > sig_10.true])/N

r.rate.1.TUK_48 = 100*length(Tn.standard.TUK_48[Tn.standard.TUK_48 > sig_05.true])/N
r.rate.2.TUK_48 = 100*length(Tn.standard.TUK_48[Tn.standard.TUK_48 > sig_10.true])/N
#########

###QS###
r.rate.1.QS_12 = 100*length(Tn.standard.QS_12[Tn.standard.QS_12 > sig_05.true])/N
r.rate.2.QS_12 = 100*length(Tn.standard.QS_12[Tn.standard.QS_12 > sig_10.true])/N

r.rate.1.QS_24 = 100*length(Tn.standard.QS_24[Tn.standard.QS_24 > sig_05.true])/N
r.rate.2.QS_24 = 100*length(Tn.standard.QS_24[Tn.standard.QS_24 > sig_10.true])/N

r.rate.1.QS_48 = 100*length(Tn.standard.QS_48[Tn.standard.QS_48 > sig_05.true])/N
r.rate.2.QS_48 = 100*length(Tn.standard.QS_48[Tn.standard.QS_48 > sig_10.true])/N
#########

size.emp.table = data.frame(
  B12_5perct = c(r.rate.1.emp.BAR_12, r.rate.1.emp.TUK_12, r.rate.1.emp.QS_12),
  B12_10perct = c(r.rate.2.emp.BAR_12, r.rate.2.emp.TUK_12, r.rate.2.emp.QS_12),
  B24_5perct = c(r.rate.1.emp.BAR_24, r.rate.1.emp.TUK_24, r.rate.1.emp.QS_24),
  B24_10perct = c(r.rate.2.emp.BAR_24, r.rate.2.emp.TUK_24, r.rate.2.emp.QS_24),
  B48_5perct = c(r.rate.1.emp.BAR_48, r.rate.1.emp.TUK_48, r.rate.1.emp.QS_48),
  B48_10perct = c(r.rate.2.emp.BAR_48, r.rate.2.emp.TUK_48, r.rate.2.emp.QS_48)
)

size.table = data.frame(
  B12_5perct = c(r.rate.1.BAR_12, r.rate.1.TUK_12, r.rate.1.QS_12),
  B12_10perct = c(r.rate.2.BAR_12, r.rate.2.TUK_12, r.rate.2.QS_12),
  B24_5perct = c(r.rate.1.BAR_24, r.rate.1.TUK_24, r.rate.1.QS_24),
  B24_10perct = c(r.rate.2.BAR_24, r.rate.2.TUK_24, r.rate.2.QS_24),
  B48_5perct = c(r.rate.1.BAR_48, r.rate.1.TUK_48, r.rate.1.QS_48),
  B48_10perct = c(r.rate.2.BAR_48, r.rate.2.TUK_48, r.rate.2.QS_48)
)