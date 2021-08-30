########################################
##The following codes are only for simulation,
##they are not used to build package.
########################################

#############
#Table of L2#
#############
###Simulation on glkerns and lokerns###
library(lokern)
n0 = 1024
i = seq(1,n0/2,by = 1)
w = 2*pi*i/n0
L2 = 0
for(times in c(1:500)){
  y.sim = arima.sim(list(order = c(2,0,2),ar = c(-1/5,-9/10),ma = c(0,1)),n = n0)
  x.sim = y.sim + 1/2*rnorm(n = n0)
  
  I = (Mod(fft(x.sim))^2/(2*pi*n0))[2:(n0/2 + 1)]
  lok = lokerns(w,I)
  glk = glkerns(w,I)
  sp.density = ((2 + 2*cos(2*w))/(1.85 + 0.76*cos(w) + 1.8*cos(2*w)) + 1/4)/(2*pi)
  
  L2 = L2 + 2*pi*sum((lok$est-sp.density)^2)/512
}

###Simulation on Haar, Franklin and Meyer###
library(lokern)
setwd("set your own work directory")
n0 = 1024
j0 = 3 #j0 = 4
j1 = 5

h.global = c(1:60)
i = seq(1,n0/2,by = 1)
w = 2*pi*i/n0

sp.ds = function(w){((2 + 2*cos(2*w))/(1.85 + 0.76*cos(w) + 1.8*cos(2*w)) + 1/4)/(2*pi)}

r.sample = function(x,h){
  n = length(x)
  x.avg = mean(x)
  x.1 = x[(abs(h)+1):n]-x.avg
  x.2 = x[1:(n-abs(h))]-x.avg
  return((t(x.1)%*%x.2)[1,1]*1/n)
}

#####################################
#We have the following Haar, Meyer and 
#Franklin three wavelets. Please only  
#run one of them for your each simulation, 
#otherwise codes would be overlapped.
#####################################

###Haar Wavelets###
w.x = function(x){sin(pi*x)/(pi*x)}

beta.jk = function(j,k,r.v){
  m = length(r.v)
  h = c(1:m)
  beta.true = 1/sqrt(2*pi)*sum(r.v*w.x(h/m)*2^(-j/2)*(sin(2*pi*h/2^(j+2)))^2/(2*pi*h/2^(j+2))*(-2*sin(2*pi*h*(2*k+1)/2^(j+1))))
  
  return(beta.true)
}

phi.jk = function(j,k,w.unit,h = 1,result = 0,st = 1){
  if(h == 1)
  {
    result.int = 2^(-j/2)*(sin(2*pi*h/2^(j+2)))^2/(2*pi*h/2^(j+2))*(-2*sin(2*pi*h*(2*k+1)/2^(j+1)-w.unit*h))
    phi.jk(j,k,w.unit,h + 1,result = result.int, st + 1)
  }
  else if(h > 1)
  {
    result.current = result + 2^(-j/2)*(sin(2*pi*h/2^(j+2)))^2/(2*pi*h/2^(j+2))*(-2*sin(2*pi*h*(2*k+1)/2^(j+1)-w.unit*h))

    if(st == 250)
    {return(result.current)}
    else
    {phi.jk(j,k,w.unit,h + 1,result.current, st + 1)}

  }
}

##Meyer Wavelets##
v.x = function(x){
  if(x <= 0)
    return(0)
  else if(x > 0 && x < 1)
    return(x^2*(3-2*x))
  else
    return(1)
}

B.x = function(h,j){
  if(abs(h) <= 2^(j+1)/3 && abs(h) >= 2^j/3)
    return(sin(pi/2*v.x(3/2^j*abs(h)-1)))
  else if(abs(h) > 2^(j+1)/3 && abs(h) <= 2^(j+2)/3)
    return(cos(pi/2*v.x(3/2^(j+1)*abs(h)-1)))
  else
    return(0)
}

w.x = function(x){sin(pi*x)/(pi*x)}

beta.jk = function(j,k,r.v){
  m = length(r.v)
  h = c(1:m)
  beta.true = 1/sqrt(2*pi)*sum(r.v*w.x(h/m)*2^(-(j-2)/2)*cos(2*pi*h*(2*k-1)/2^(j+1))*sapply(h,B.x, j = j))
  
  return(beta.true)
}

phi.jk = function(j,k,w.unit,h = 1,result = 0,st = 1){
  if(h == 1)
  {
    result.int = 2^(-(j-2)/2)*cos(2*pi*h*(2*k-1)/2^(j+1)-w.unit*h)*B.x(h,j)
    phi.jk(j,k,w.unit,h + 1,result = result.int, st + 1)
  }
  else if(h > 1)
  {
    result.current = result + 2^(-(j-2)/2)*cos(2*pi*h*(2*k-1)/2^(j+1)-w.unit*h)*B.x(h,j)
    
    if(st == 150)
    {return(result.current)}
    else
    {phi.jk(j,k,w.unit,h + 1,result.current, st + 1)}
    
  }
}

##Franklin Wavelets##
w.x = function(x){sin(pi*x)/(pi*x)}

beta.jk = function(j,k,r.v){
  m = length(r.v)
  h = c(1:m)
  beta.true = 1/sqrt(2*pi)*sum(r.v*w.x(h/m)*2^(3/2*j + 5)*cos(2*pi*h*(2*k-1)/2^(j+1))*(sin(2*pi*h/2^(j+2)))^4/(2*pi*h)^2*((1-2/3*(cos(2*pi*h/2^(j+2)))^2)/((1-2/3*(sin(2*pi*h/2^(j+1)))^2)*(1-2/3*(sin(2*pi*h/2^(j+2)))^2)))^(1/2))
  
  return(beta.true)
}

phi.jk = function(j,k,w.unit,h = 1,result = 0,st = 1){
  if(h == 1)
  {
    result.int = 2^(3/2*j + 5)*cos(2*pi*h*(2*k-1)/2^(j+1) - w.unit*h)*(sin(2*pi*h/2^(j+2)))^4/(2*pi*h)^2*((1-2/3*(cos(2*pi*h/2^(j+2)))^2)/((1-2/3*(sin(2*pi*h/2^(j+1)))^2)*(1-2/3*(sin(2*pi*h/2^(j+2)))^2)))^(1/2)
    phi.jk(j,k,w.unit,h + 1,result = result.int, st + 1)
  }
  else if(h > 1)
  {
    result.current = result + 2^(3/2*j + 5)*cos(2*pi*h*(2*k-1)/2^(j+1) - w.unit*h)*(sin(2*pi*h/2^(j+2)))^4/(2*pi*h)^2*((1-2/3*(cos(2*pi*h/2^(j+2)))^2)/((1-2/3*(sin(2*pi*h/2^(j+1)))^2)*(1-2/3*(sin(2*pi*h/2^(j+2)))^2)))^(1/2)
    
    if(st == 200)
    {return(result.current)}
    else
    {phi.jk(j,k,w.unit,h + 1,result.current, st + 1)}
    
  }
}

#####Theoretic Threshold#####
#Based on simulation, we find in terms of simulation performance,
#there's no big difference between true and estimated threshold,
#in consideration of speed issue, the following theoretic threshold
#is highly recommended.
#############################
delta.sq.jk.integrand = function(w.unit, j, k){
  (phi.jk(j = j,k = k,w.unit = w.unit) + phi.jk(j = j,k = k,w.unit = -w.unit))^2*(sp.ds(w.unit))^2
}

delta.sq.jk = function(j,k){
  integrate(delta.sq.jk.integrand,-pi,pi,j = j,k = k)$value/(2*n0)
}

lambda = 2^(j1+1)-2^j0

#Sigma1
sigma.sq.jk = function(j,k){
  2*log(lambda)*delta.sq.jk(j = j,k = k)
}

#Sigma2
# sigma.sq.jk = function(j,k){
#   2*log(lambda/(2^j0))*delta.sq.jk(j = j,k = k)
# }

##Simulation Starts##
L2 = 0
for(times in c(1:500)){
  y.sim.meyer = arima.sim(list(order = c(2,0,2),ar = c(-1/5,-9/10),ma = c(0,1)),n = n0)
  x.sim.meyer = y.sim.meyer + 1/2*rnorm(n = n0)
  
  #############################
  #The following commented codes are
  #for estimated threshold, but it
  #takes long time to run. It is 
  #not recommended since there's
  #no difference between estimated 
  #and theoretic threshold.
  #############################
  #######Sigma2 Est########
  # I = (Mod(fft(x.sim.meyer))^2/(2*pi*n0))[2:(n0/2 + 1)]
  # lok = lokerns(w,I)
  # lok.est = lok$est
  # lok.all = c(lok.est[1],lok.est)
  # w.all = c(0,w)
  # 
  # lok.final = midpoints.add(lok.all)
  # w.final = midpoints.add(w.all)
  # 
  # 
  # sp.ds.est = function(w.unit){
  #   index = which(abs(w.unit - w.final) <= pi/(2^8*n0))[1]
  #   lok.final[index]
  # }
  # 
  # delta.sq.jk.integrand = function(w.unit, j, k){
  #   sp.ds.est.v = Vectorize(sp.ds.est,vectorize.args = "w.unit")
  # 
  #   (phi.jk(j = j,k = k,w.unit = w.unit) + phi.jk(j = j,k = k,w.unit = -w.unit))^2*(sp.ds.est.v(w.unit))^2
  # }
  # 
  # delta.sq.jk.est = function(j,k){
  # 
  #   withCallingHandlers(
  #     withRestarts({
  #       withCallingHandlers(
  #         withRestarts({
  #           withCallingHandlers(
  #             withRestarts({
  #               return(integrate(delta.sq.jk.integrand,0,pi,j = j,k = k,subdivisions = 500)$value/n0)
  #             },
  #             change_two_limits = function(){
  #               return((integrate(delta.sq.jk.integrand,0,pi/2,j = j,k = k,subdivisions = 500)$value
  #                       +integrate(delta.sq.jk.integrand,pi/2,pi,j = j,k = k,subdivisions = 500)$value)/n0)
  #             }
  #             ),error = function(e) {invokeRestart("change_two_limits")})},
  #           change_three_limits = function(){
  #             return((integrate(delta.sq.jk.integrand,0,pi/3,j = j,k = k,subdivisions = 500)$value
  #                     +integrate(delta.sq.jk.integrand,pi/3,2*pi/3,j = j,k = k,subdivisions = 500)$value
  #                     +integrate(delta.sq.jk.integrand,2*pi/3,pi,j = j,k = k,subdivisions = 500)$value)/n0)
  #           }
  #         ),error = function(e) {invokeRestart("change_three_limits")})},
  #       change_four_limits = function(){
  #         return((integrate(delta.sq.jk.integrand,0,pi/4,j = j,k = k,subdivisions = 500)$value
  #                 +integrate(delta.sq.jk.integrand,pi/4,pi/2,j = j,k = k,subdivisions = 500)$value
  #                 +integrate(delta.sq.jk.integrand,pi/2,3*pi/4,j = j,k = k,subdivisions = 500)$value
  #                 +integrate(delta.sq.jk.integrand,3*pi/4,pi,j = j,k = k,subdivisions = 500)$value)/n0)
  #       }
  #     ),error = function(e){invokeRestart("change_four_limits")})
  # }
  # #Sigma1
  # sigma.sq.jk = function(j,k){
  #   2*log(lambda)*delta.sq.jk.est(j = j,k = k)
  # }
  # #Sigma2
  # # sigma.sq.jk = function(j,k){
  # #   2*log(lambda/2^j0)*delta.sq.jk.est(j = j,k = k)
  # # }
  #########################
  
  Alpha = 1/sqrt(2*pi)*r.sample(x.sim.meyer,h = 0)
  
  r.x.sample = sapply(h.global, r.sample, x = x.sim.meyer)
  
  psum1.haar = beta.jk(0,0,r.x.sample)*sapply(w, phi.jk, j = 0, k = 0)/sqrt(2*pi)
  
  for (j in c(1:(j0-1))) {
    
    k = c(0:(2^j-1))
    ####Sample Model####
    beta.result = matrix(sapply(k,beta.jk,j = j, r.v = r.x.sample), nrow = n0/2, ncol = length(k),byrow = TRUE)

    phi.result = sapply(k, function(t){sapply(w,phi.jk,j = j,k = t)})
    psum1.haar = psum1.haar + colSums(t(beta.result)*t(phi.result))*1/sqrt(2*pi)
    
  }
  
  psum2.haar = rep(0,n0/2)
  result.in.loop = rep(0,n0/2)

  for (j in c(j0:j1)) {
    for (k in c(0:(2^j-1))) {
      if(abs(beta.jk(j = j,k = k,r.x.sample)) > sigma.sq.jk(j = j,k = k)){
        result.in.loop = beta.jk(j = j,k = k,r.x.sample)*sapply(w,phi.jk,j = j,k = k)/sqrt(2*pi)
       } else {
        result.in.loop = rep(0,n0/2)
      }
      psum2.haar = psum2.haar + result.in.loop
    }
  }
  
  f.haar = 1/sqrt(2*pi)*Alpha + psum1.haar + psum2.haar
  sp.density = sp.ds(w)
  ###Codes are used to generate plot for each simulation###
  name = paste("replicateFranklin.esth60",times,".png",sep="")
  png(name)
  plot(sp.density,type = "l",xlab = "Index of Fourier frequency",lty = 3,ylim = c(-.15,1))
  lines(f.haar,type = "l",lty = 1,col = 'red')
  dev.off()
  #########################################################
  L2 = L2 + 2*pi*sum((f.haar-sp.density)^2)/512
}

####DAUB 8 Wavelets####
library(fourierin)
library(lokern)
library(wavethresh)

setwd("set your own work directory")

#This is used to save the results of true threshold
sigma1.vector = rep(0,56)#length is 48 if j0 = 4
# sigma2.vector = rep(0,56)

n0 = 1024
h.global = c(1:60)
j0 = 3#j0 = 4
j1 = 5

i = seq(1,n0/2,by = 1)
w = 2*pi*i/n0

sp.ds = function(w.unit){((2 + 2*cos(2*w.unit))/(1.85 + 0.76*cos(w.unit) + 1.8*cos(2*w.unit)) + 1/4)/(2*pi)}

r.sample = function(x,h){
  n = length(x)
  x.avg = mean(x)
  x.1 = x[(abs(h)+1):n]-x.avg
  x.2 = x[1:(n-abs(h))]-x.avg
  return((t(x.1)%*%x.2)[1,1]*1/n)
}

##D phase Wavelets##
filters = filter.select(filter.number = 8,family = "DaubExPhase")
h.filter = filters$H
L = length(h.filter)
matrix.filter.scale = function(h.index){
  if (h.index %in% c(0:(L-1))){
    return(sqrt(2)*h.filter[h.index+1])
  } else {
    return(0)
  }
}
vector.filter.wavle = function(h.index,x){
  if ((h.index+1-floor(2*x)) %in% c(0:(L-1))){
    return((-1)^(floor(2*x)-h.index)*h.filter[h.index+1-floor(2*x)+1])
  } else {
    return(0)
  }
}

T.row = matrix(c(1:(L-1)),nrow = (L-1),ncol = (L-1))
T.col = t(T.row)
T.0 = 2*T.row - T.col - 1
T.1 = 2*T.row - T.col

T.0[] = vapply(T.0,matrix.filter.scale,numeric(1))
T.1[] = vapply(T.1,matrix.filter.scale,numeric(1))


dyad = function(x,n){
  dyad.repst = rep(-1,n)
  x = x - floor(x)
  for(i in c(1:n)){
    if(2*x >= 1){
      dyad.repst[i] = 1
      x = 2*x - 1
    } else {
      dyad.repst[i] = 0
      x = 2*x
    }
  }
  return(dyad.repst)
}

u.x = function(x){
  miu = c(0:(L-2))
  return(sapply(miu,vector.filter.wavle,x = x))
}

v.x = function(x,n){
  T.f = diag(L-1)
  repst = dyad(2*x,n)
  for(i in repst){
    if(i == 1){
      T.f = T.f %*% T.1
    } else {
      T.f = T.f %*% T.0
    }
  }
  return(sqrt(2)/(L-1)*T.f%*%rep(1,L-1))
}

D.wal = function(x){
  u.x(x) %*% v.x(x,30)
}
###################

rol.increase = function(r){
  if(r < 5000){
    rol.increase(2*r - 1)
  } else {
    return(r)
  }
}

phi.s.fourier = function(w.v){
  rol = rol.increase(length(w.v))
  D.wal.v = Vectorize(D.wal)
  #D4:-3  4, D5:-4  5, D6:-5  6, D7:-6  7,D8:-7  8,D9:-8  9,D10:-9 10
  return(-fourierin(f = D.wal.v, lower_int = -7,upper_int = 8,eval_grid = w.v,const_adj = 0,freq_adj = -1,resolution = rol))
}

phi.jk.s.fourier = function(t,j,k){
  (cos(t*k/2^j) - sin(t*k/2^j)*1i)*2^(-j/2)*phi.s.fourier(t/2^j)
}

phi.jk.L.fourier = function(h,j,k){
  sqrt(2*pi)*phi.jk.s.fourier(2*pi*h,j = j,k = k)
}

w.x = function(x){sin(pi*x)/(pi*x)}

beta.jk = function(j,k,r.v){
  m = length(r.v)
  h = c(1:m)
  beta.true = 1/sqrt(2*pi)*sum(r.v*w.x(h/m)*(Re(phi.jk.L.fourier(h,j = j,k = k)+phi.jk.L.fourier(-h,j = j,k = k))))
  
  return(beta.true)
}

comp.exp = function(w.unit,h.unit){
  cos(w.unit*h.unit) + sin(w.unit*h.unit)*1i
}

phi.jk = function(j,k,w.v){
  h = c(1:500)
  
  result = (t(phi.jk.L.fourier(h,j = j,k = k)) %*% sapply(w.v, comp.exp,h.unit = h)) + (t(phi.jk.L.fourier(-h,j = j,k = k)) %*% sapply(w.v, comp.exp,h.unit = -h)) 
  
  return(t(Re(result)))
}

phi.jk.integrand = function(j,k,w.unit,phi.1,phi.2,h.0){
  result = sum(Re(phi.1*comp.exp(w.unit = w.unit,h.unit = h.0) + phi.2*comp.exp(w.unit = w.unit,h.unit = -h.0)))
  
  return(result)
}

lambda = 2^(j1+1)-2^j0
############Theoretic Sigma1 & Sigma2##################
delta.sq.jk.integrand = function(w.unit, j, k, phi.1, phi.2, h.0){
  phi.jk.integrand.v = Vectorize(phi.jk.integrand,vectorize.args = "w.unit")
  
  (phi.jk.integrand.v(j = j,k = k,w.unit = w.unit,phi.1 = phi.1,phi.2 = phi.2,h.0 = h.0) + phi.jk.integrand.v(j = j,k = k,w.unit = -w.unit,phi.1 = phi.1,phi.2 = phi.2,h.0 = h.0))^2*(sp.ds(w.unit))^2
}

delta.sq.jk = function(j,k){
  h0 = c(1:500)
  phi1 = phi.jk.L.fourier(h0,j = j,k = k)
  phi2 = phi.jk.L.fourier(-h0,j = j,k = k)
  
  2*integrate(delta.sq.jk.integrand,0,pi,j = j,k = k,phi.1 = phi1,phi.2 = phi2,h.0 = h0)$value/(2*n0)
}

sigma.sq.jk = function(j,k){
  2*log(lambda)*delta.sq.jk(j = j,k = k)
}

# sigma.sq.jk = function(j,k){
#   2*log(lambda/(2^j0))*delta.sq.jk(j = j,k = k)
# }
##############################################

L2 = 0

for(times in c(1:500)){
  y.sim.meyer = arima.sim(list(order = c(2,0,2),ar = c(-1/5,-9/10),ma = c(0,1)),n = n0)
  x.sim.meyer = y.sim.meyer + 1/2*rnorm(n = n0)
  
  Alpha = 1/sqrt(2*pi)*r.sample(x.sim.meyer,h = 0)
  
  r.x.sample = sapply(h.global, r.sample, x = x.sim.meyer)
  
  psum1.D = beta.jk(0,0,r.x.sample)*phi.jk(0,0,w)/sqrt(2*pi)
  
  for (j in c(1:(j0-1))) {
    k = c(0:(2^j-1))
    beta.result = matrix(sapply(k,beta.jk,j = j, r.v = r.x.sample), nrow = n0/2, ncol = length(k),byrow = TRUE)
    phi.result = sapply(k, function(t){phi.jk(j = j,k = t,w.v = w)})
    psum1.D = psum1.D + colSums(t(beta.result)*t(phi.result))*1/sqrt(2*pi)
  }
  
  psum2.D = rep(0,n0/2)
  result.in.loop = rep(0,n0/2)
  col.index = 0
  for (j in c(j0:j1)) {
    for (k in c(0:(2^j-1))) {
      
      ####Sample Model####
      col.index = col.index + 1
      if(times == 1){
        sigma1.vector[col.index] = sigma.sq.jk(j = j,k = k)
      }
      if(abs(beta.jk(j = j,k = k,r.x.sample)) > sigma1.vector[col.index]){
        result.in.loop = beta.jk(j = j,k = k,r.x.sample)*phi.jk(j = j,k = k,w.v = w)/sqrt(2*pi)
      } else {
        result.in.loop = rep(0,n0/2)
      }
      
      psum2.D = psum2.D + result.in.loop
    }
  }
  
  f.D = 1/sqrt(2*pi)*Alpha + psum1.D + psum2.D
  sp.density = sp.ds(w)
  ###Codes are used to generate plot for each simulation###
  name = paste("replicateD8w(60)",times,".png",sep="")
  png(name)
  plot(sp.density,type = "l",ylab = "",cex.lab = .75,xlab = "Index of Fourier frequency",lty = 3,ylim = c(-.15,1))
  lines(f.D,type = "l",lty = 1)
  dev.off()
  #########################################################
  L2 = L2 + 2*pi*sum((f.D-sp.density)^2)/512
}

###Plot of Simulation###
####DAUB 8 Wavelets####
library(fourierin)
library(lokern)
library(wavethresh)

setwd("H:/Nonparametric/PlotOfSDE/Temporary")

n0 = 1024
h.global = c(1:60)
j0 = 3#j0 = 4
j1 = 5

i = seq(1,n0/2,by = 1)
w = 2*pi*i/n0

sp.ds = function(w.unit){((2 + 2*cos(2*w.unit))/(1.85 + 0.76*cos(w.unit) + 1.8*cos(2*w.unit)) + 1/4)/(2*pi)}

r.sample = function(x,h){
  n = length(x)
  x.avg = mean(x)
  x.1 = x[(abs(h)+1):n]-x.avg
  x.2 = x[1:(n-abs(h))]-x.avg
  return((t(x.1)%*%x.2)[1,1]*1/n)
}

##D phase Wavelets##
filters = filter.select(filter.number = 8,family = "DaubExPhase")
h.filter = filters$H
L = length(h.filter)
matrix.filter.scale = function(h.index){
  if (h.index %in% c(0:(L-1))){
    return(sqrt(2)*h.filter[h.index+1])
  } else {
    return(0)
  }
}
vector.filter.wavle = function(h.index,x){
  if ((h.index+1-floor(2*x)) %in% c(0:(L-1))){
    return((-1)^(floor(2*x)-h.index)*h.filter[h.index+1-floor(2*x)+1])
  } else {
    return(0)
  }
}

T.row = matrix(c(1:(L-1)),nrow = (L-1),ncol = (L-1))
T.col = t(T.row)
T.0 = 2*T.row - T.col - 1
T.1 = 2*T.row - T.col

T.0[] = vapply(T.0,matrix.filter.scale,numeric(1))
T.1[] = vapply(T.1,matrix.filter.scale,numeric(1))


dyad = function(x,n){
  dyad.repst = rep(-1,n)
  x = x - floor(x)
  for(i in c(1:n)){
    if(2*x >= 1){
      dyad.repst[i] = 1
      x = 2*x - 1
    } else {
      dyad.repst[i] = 0
      x = 2*x
    }
  }
  return(dyad.repst)
}

u.x = function(x){
  miu = c(0:(L-2))
  return(sapply(miu,vector.filter.wavle,x = x))
}

v.x = function(x,n){
  T.f = diag(L-1)
  repst = dyad(2*x,n)
  for(i in repst){
    if(i == 1){
      T.f = T.f %*% T.1
    } else {
      T.f = T.f %*% T.0
    }
  }
  return(sqrt(2)/(L-1)*T.f%*%rep(1,L-1))
}

D.wal = function(x){
  u.x(x) %*% v.x(x,30)
}
###################

rol.increase = function(r){
  if(r < 5000){
    rol.increase(2*r - 1)
  } else {
    return(r)
  }
}

phi.s.fourier = function(w.v){
  rol = rol.increase(length(w.v))
  D.wal.v = Vectorize(D.wal)
  #D4:-3  4, D5:-4  5, D6:-5  6, D7:-6  7,D8:-7  8,D9:-8  9,D10:-9 10
  return(-fourierin(f = D.wal.v, lower_int = -7,upper_int = 8,eval_grid = w.v,const_adj = 0,freq_adj = -1,resolution = rol))
}

phi.jk.s.fourier = function(t,j,k){
  (cos(t*k/2^j) - sin(t*k/2^j)*1i)*2^(-j/2)*phi.s.fourier(t/2^j)
}

phi.jk.L.fourier = function(h,j,k){
  sqrt(2*pi)*phi.jk.s.fourier(2*pi*h,j = j,k = k)
}

w.x = function(x){sin(pi*x)/(pi*x)}

beta.jk = function(j,k,r.v){
  m = length(r.v)
  h = c(1:m)
  beta.true = 1/sqrt(2*pi)*sum(r.v*w.x(h/m)*(Re(phi.jk.L.fourier(h,j = j,k = k)+phi.jk.L.fourier(-h,j = j,k = k))))
  
  return(beta.true)
}

comp.exp = function(w.unit,h.unit){
  cos(w.unit*h.unit) + sin(w.unit*h.unit)*1i
}

phi.jk = function(j,k,w.v){
  h = c(1:500)
  
  result = (t(phi.jk.L.fourier(h,j = j,k = k)) %*% sapply(w.v, comp.exp,h.unit = h)) + (t(phi.jk.L.fourier(-h,j = j,k = k)) %*% sapply(w.v, comp.exp,h.unit = -h)) 
  
  return(t(Re(result)))
}

phi.jk.integrand = function(j,k,w.unit,phi.1,phi.2,h.0){
  result = sum(Re(phi.1*comp.exp(w.unit = w.unit,h.unit = h.0) + phi.2*comp.exp(w.unit = w.unit,h.unit = -h.0)))
  
  return(result)
}

lambda = 2^(j1+1)-2^j0
############Theoretic Sigma2##################
delta.sq.jk.integrand = function(w.unit, j, k, phi.1, phi.2, h.0){
  phi.jk.integrand.v = Vectorize(phi.jk.integrand,vectorize.args = "w.unit")
  
  (phi.jk.integrand.v(j = j,k = k,w.unit = w.unit,phi.1 = phi.1,phi.2 = phi.2,h.0 = h.0) + phi.jk.integrand.v(j = j,k = k,w.unit = -w.unit,phi.1 = phi.1,phi.2 = phi.2,h.0 = h.0))^2*(sp.ds(w.unit))^2
}

delta.sq.jk = function(j,k){
  h0 = c(1:500)
  phi1 = phi.jk.L.fourier(h0,j = j,k = k)
  phi2 = phi.jk.L.fourier(-h0,j = j,k = k)
  
  2*integrate(delta.sq.jk.integrand,0,pi,j = j,k = k,phi.1 = phi1,phi.2 = phi2,h.0 = h0)$value/(2*n0)
}

sigma.sq.jk = function(j,k){
  2*log(lambda/(2^j0))*delta.sq.jk(j = j,k = k)
}
##############################################
#Read the data to replicate results
times = 24
replicate.x.meyer.table = read.table("result(24).txt",sep = "\t")
x.sim.meyer = as.matrix(replicate.x.meyer.table[times,])

###Plot of glkern and lokern###
I = (Mod(fft(x.sim.meyer))^2/(2*pi*n0))[2:(n0/2 + 1)]
lok = lokerns(w,I)
glk = glkerns(w,I)
  
Alpha = 1/sqrt(2*pi)*r.sample(x.sim.meyer,h = 0)
  
r.x.sample = sapply(h.global, r.sample, x = x.sim.meyer)
  
psum1.D = beta.jk(0,0,r.x.sample)*phi.jk(0,0,w)/sqrt(2*pi)
  
for (j in c(1:(j0-1))) {
  k = c(0:(2^j-1))
  beta.result = matrix(sapply(k,beta.jk,j = j, r.v = r.x.sample), nrow = n0/2, ncol = length(k),byrow = TRUE)
  phi.result = sapply(k, function(t){phi.jk(j = j,k = t,w.v = w)})
  psum1.D = psum1.D + colSums(t(beta.result)*t(phi.result))*1/sqrt(2*pi)
}
  
psum2.D = rep(0,n0/2)
result.in.loop = rep(0,n0/2)
col.index = 0
for (j in c(j0:j1)) {
  for (k in c(0:(2^j-1))) {
      
    ####Sample Model####
    col.index = col.index + 1
    if(abs(beta.jk(j = j,k = k,r.x.sample)) > sigma.sq.jk(j = j,k = k)){
      result.in.loop = beta.jk(j = j,k = k,r.x.sample)*phi.jk(j = j,k = k,w.v = w)/sqrt(2*pi)
    } else {
      result.in.loop = rep(0,n0/2)
    }
      
    psum2.D = psum2.D + result.in.loop
  }
}
  
f.D = 1/sqrt(2*pi)*Alpha + psum1.D + psum2.D
sp.density = sp.ds(w)

par(mfrow = c(1,2),mar=c(5,2.5,1.75,1.25))
plot(sp.density,type = "l",ylab = "",cex.lab = .75,xlab = "Index of Fourier frequency",lty = 3,ylim = c(-.15,1))
lines(lok$est,type = "l",lty = 1)
lines(glk$est,type = "l",lty = 2)
plot(sp.density,type = "l",ylab = "",cex.lab = .75,xlab = "Index of Fourier frequency",lty = 3,ylim = c(-.15,1))
lines(f.D,type = "l",lty = 1)
