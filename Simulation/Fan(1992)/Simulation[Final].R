#######################
##Table of Statistics##
#######################

# percent function is used, therefore we need to run 
# this package in front of codes.
library(scales)

# Warn: You should only input integer 1, 2 or 3. Since
# we only have those three simulations in total.
Simulation = 3

# The value you input is not 1, 2, or 3, the code stops
# and reports you an error !
if(!(Simulation %in% c(1,2,3)))
{stop('Simulation must be 1, 2, or 3 !')}

###Common Function###
Kernel.normal.s = function(x,miu=0,sd=1,l=1){(exp(-(x-miu)^2/(2*sd^2))/(sqrt(2*pi)*sd))^(l)}


################
##Simulation 1##
################
if(Simulation == 1)
{
integrand.K.var1.s = function(x){Kernel.normal.s(x,0,1,2)}
integrand.K.var2.s = function(x,deltax=.09){deltax/(Kernel.normal.s(x,0,delta.s,1))}
integrand.K.bias1.s = function(x){x^2*Kernel.normal.s(x)}
integrand.K.bias2.s = function(x){(-(.75)^2*sin(.75*x))^2}
integrand.K.bias2.s.nw = function(x){(-(.75)^2*sin(.75*x)+2*.75*cos(.75*x)*(-x/(delta.s)^2))^2}
h.formula = function(n,del,delta.x){(1/n)^(1/5)*(integrate(integrand.K.var1.s,-Inf,Inf)$value*integrate(integrand.K.var2.s,-2*del,2*del,deltax=delta.x)$value/((integrate(integrand.K.bias1.s,-Inf,Inf)$value)^2*integrate(integrand.K.bias2.s,-2*del,2*del)$value))^(1/5)}
h.formula.nw = function(n,del){(1/n)^(1/5)*(integrate(integrand.K.var1.s,-Inf,Inf)$value*integrate(integrand.K.var2.s,-2*del,2*del)$value/((integrate(integrand.K.bias1.s,-Inf,Inf)$value)^2*integrate(integrand.K.bias2.s.nw,-2*del,2*del)$value))^(1/5)}
}

################
##Simulation 2##
################
if(Simulation == 2)
{
integrand.K.var1.s = function(x){Kernel.normal.s(x,0,1,2)}
integrand.K.var2.s = function(x,deltax=.16){deltax/(0.5*Kernel.normal.s(x,-1,1,1)+0.5*Kernel.normal.s(x,1.75,.25,1))}
integrand.K.bias1.s = function(x){x^2*Kernel.normal.s(x)}
integrand.K.bias2.s = function(x){(-.5*(2.5)^2*sin(2.5*x))^2}
integrand.K.bias2.s.nw = function(x){(-.5*(2.5)^2*sin(2.5*x)+2.5*cos(2.5*x)*(-0.5*Kernel.normal.s(x,-1,1,1)*(x+1)-.5*Kernel.normal.s(x,1.75,.25,1)*(x-1.75)/(.25^2))/(0.5*Kernel.normal.s(x,-1,1,1)+0.5*Kernel.normal.s(x,1.75,.25,1)))^2}
h.formula = function(n,del,delta.x){(1/n)^(1/5)*(integrate(integrand.K.var1.s,-Inf,Inf)$value*integrate(integrand.K.var2.s,-2*del,2*del,deltax=delta.x)$value/(4*(integrate(integrand.K.bias1.s,-Inf,Inf)$value)^2*integrate(integrand.K.bias2.s,-2*del,2*del)$value))^(1/5)}
h.formula.nw = function(n,del){(1/n)^(1/5)*(integrate(integrand.K.var1.s,-Inf,Inf)$value*integrate(integrand.K.var2.s,-2*del,2*del)$value/(4*(integrate(integrand.K.bias1.s,-Inf,Inf)$value)^2*integrate(integrand.K.bias2.s.nw,-2*del,2*del)$value))^(1/5)}
}

################
##Simulation 3##
################
if(Simulation == 3)
{
integrand.K.var1.s = function(x){Kernel.normal.s(x,0,1,2)}
integrand.K.bias1.s = function(x){x^2*Kernel.normal.s(x)}
integrand.K.bias2.s = function(x){(-.5*(2.5)^2*sin(2.5*x))^2}
h.formula = function(n,del,delta.x){(1/n)^(1/5)*(integrate(integrand.K.var1.s,-Inf,Inf)$value*(4*del*delta.x/.2)/(4*(integrate(integrand.K.bias1.s,-Inf,Inf)$value)^2*integrate(integrand.K.bias2.s,-2*del,2*del)$value))^(1/5)}
}

##local linear smoother##
mise.lls.s = rep(0,300)
s.n.l.s = function(x,X,h,l){sum(Kernel.normal.s((x-X)/h)*(x-X)^(l))}
est.lls = function(x,X=X.r.lls.s,h=h.ls.s,y=y.r.lls.s){
  w = Kernel.normal.s((x-X)/h)*(s.n.l.s(x,X,h,2)-(x-X)*s.n.l.s(x,X,h,1))
  result = sum(w*y)/sum(w)
  #Simulation 1
  if(Simulation == 1)
  {return((result-sin(.75*x))^2)}
  
  #Simulation 2 & 3
  if(Simulation == 2 || Simulation == 3)
  {return((result-sin(2.5*x))^2)}
}

##Gasser-Muller##
#library(pracma)
mise.gm.s = rep(0,300)
est.gm = function(x,y=y.r.gm.s,h=h.gm.s,tt=t.gm.s){
  result = sapply(x,function(z){y[1]*integrate(function(t){1/h*Kernel.normal.s((z-t)/h)},tt[1],tt[2])$value})
  for(j in 3:(n0+1)){
    result = result + sapply(x,function(z){y[j-1]*integrate(function(t){1/h*Kernel.normal.s((z-t)/h)},tt[j-1],tt[j])$value})
  }
  #Simulation 1
  if(Simulation == 1)
  {return((result - sin(.75*x))^2)}
  #Simulation 2 & 3
  if(Simulation == 2 || Simulation == 3)
  {return((result - sin(2.5*x))^2)}
}

################################################
####Calculation has been done by Riemman Sum####
################################################
# x.s1 = seq(-2*delta.gm.s1,2*delta.gm.s1,by = .01)
# integral.K=function(x,j,h,tt){
#   sapply(x,function(z){integrate(function(t){1/h*exp(-((z-t)/h)^2/2)/sqrt(2*pi)},tt[j-1],tt[j])$value})
# }
# 
# m.est = function(x,y,h,tt){
#   result = y[1]*integral.K(x,2,h,tt)
#   for(j in 3:(n2+1)){
#     result = result + y[j-1]*integral.K(x,j,h,tt)
#   }
#   return(result)
# }
# 
# diff.m = function(x,y,h,tt){(m.est(x,y,h,tt)-sin(.75*x))^2}

##Nadaraya-Waston##
mise.nw.s = rep(0,300)
#Simulation 1
if(Simulation == 1)
{est.nw = function(x,X=X.r.nw.s,h=h.nw.s){(sum(Kernel.normal.s((x-X)/h)*y.r.nw.s)/sum(Kernel.normal.s((x-X)/h))-sin(.75*x))^2}}
#Simulation 2 & 3
if(Simulation == 2 || Simulation == 3)
{est.nw = function(x,X=X.r.nw.s,h=h.nw.s){(sum(Kernel.normal.s((x-X)/h)*y.r.nw.s)/sum(Kernel.normal.s((x-X)/h))-sin(2.5*x))^2}}


###Table###
#Simulation 1
if(Simulation == 1)
{output=data.frame(matrix(0,ncol = 7,nrow = 6))}
#Simulation 2 & 3
if(Simulation == 2 || Simulation == 3)
{output=data.frame(matrix(0,ncol = 7,nrow = 3))}
names(output)=c("variance","Sample.Size","LLs.MISE","GM.MISE","Eff.GM","NW.MISE","Eff.NW")

if(Simulation == 1)
{n.rows = 6}
if(Simulation == 2 || Simulation == 3)
{n.rows = 3}


for(s in 1:n.rows){

###Common Parameter###
  
#Simulation 1
if(Simulation == 1)
{
if(s>=1 & s<=3)
{
  delta.s = .25
  if(s == 1)  {n0 = 100}
  else if(s == 2) {n0 = 200}
  else {n0 = 400}
}
else{
  delta.s = 1
  if(s == 4)  {n0 = 100}
  else if(s == 5) {n0 = 200}
  else {n0 = 400}
}
} 
#Simulation 2 & 3
if(Simulation == 2 || Simulation ==3)
{
  delta.s = 1
  if(s == 1)  {n0 = 100}
  else if(s == 2) {n0 = 200}
  else {n0 = 400}
}
for(i in 1:300){
  ###Local linear Smoother Simulation###
  #Simulation 1
  if(Simulation == 1)
  {
  X.r.lls.s = rnorm(n0,0,delta.s)
  error.lls.s = rnorm(n0)
  y.r.lls.s = sin(.75*X.r.lls.s)+.3*error.lls.s
  h.ls.s = h.formula(n0,delta.s,.09)
  }
  

  if(Simulation == 2 || Simulation == 3)
  {
  #Simulation 2
  if(Simulation == 2)
  {
  coeff = rbinom(n0,1,.5)
  X.r.lls.s = coeff*rnorm(n0,-1,1)+(1-coeff)*rnorm(n0,1.75,.25)
  }
  #Simulation 3
  else{
  X.r.lls.s = runif(n0,min = -2.5,max = 2.5)
  }
    
  error.lls.s = rnorm(n0)
  y.r.lls.s = sin(2.5*X.r.lls.s)+.4*error.lls.s
  h.ls.s = h.formula(n0,delta.s,.16)
  }

  
  mise.lls.s[i]=integrate(Vectorize(est.lls,vectorize.args = 'x'),-2*delta.s,2*delta.s)$value
  
  ###Gasser-Muller Simulation###
  X.r.gm.s = sort(X.r.lls.s,decreasing = FALSE)
  t.gm.s = c(-Inf,rep(0,n0-1),Inf)
  for(k in 2:n0)
  {
    t.gm.s[k] = (X.r.gm.s[k-1]+X.r.gm.s[k])/2
  }
  error.gm.s = error.lls.s
  #Simulation 1
  if(Simulation == 1)
  {
  y.r.gm.s = sin(.75*X.r.gm.s)+.3*error.gm.s
  h.gm.s = h.formula(n0,delta.s,.09*1.5)
  }
  #Simulation 2 & 3
  if(Simulation == 2 || Simulation ==3)
  {
  y.r.gm.s = sin(2.5*X.r.gm.s)+.4*error.gm.s
  h.gm.s = h.formula(n0,delta.s,.16*1.5)
  }
  mise.gm.s[i] = integrate(est.gm,-2*delta.s,2*delta.s)$value
  
  
  ###Nadaraya-Waston Simulation###
  X.r.nw.s = X.r.lls.s
  error.nw.s = error.lls.s
  #Simulation 1
  if(Simulation == 1)
  {
    y.r.nw.s = sin(.75*X.r.nw.s)+.3*error.nw.s
    h.nw.s = h.formula.nw(n0,delta.s)
  }
  #Simulation 2 & 3
  if(Simulation == 2 || Simulation == 3)
  {
    y.r.nw.s = sin(2.5*X.r.nw.s)+.4*error.nw.s
    if(Simulation == 2)
    {h.nw.s = h.formula.nw(n0,delta.s)}
    else
    {h.nw.s = h.formula(n0,delta.s,.16)}
  }
  
  mise.nw.s[i] = integrate(Vectorize(est.nw,vectorize.args = 'x'),-2*delta.s,2*delta.s)$value

}
###Calcualte MISE and Efficiency###
mise.lls.s.mean = mean(mise.lls.s)

mise.gm.s.mean = mean(mise.gm.s)
eff.gm.s = (mise.lls.s.mean/mise.gm.s.mean)^(5/4)

mise.nw.s.mean = mean(mise.nw.s)
eff.nw.s = (mise.lls.s.mean/mise.nw.s.mean)^(5/4)

###Store Results###
output[s,]=c(delta.s,n0,round(mise.lls.s.mean,digits = 4),round(mise.gm.s.mean,digits = 4),percent(eff.gm.s),round(mise.nw.s.mean,digits = 4),percent(eff.nw.s))
}

print(output)


#########
##Plots##
#########

###Common Variables###
x = seq(-0.6,0.6,by=0.01)
n = 100
delta = .25

###Common Functions###
Kernel.normal = function(x,miu=0,sd=1,l=1){(exp(-(x-miu)^2/(2*sd^2))/(sqrt(2*pi)*sd))^(l)}
integrand.K.var1 = function(x){Kernel.normal(x,0,1,2)}
integrand.K.var2 = function(x,deltax=.09){deltax/(Kernel.normal(x,0,delta,1))}
integrand.K.bias1 = function(x){x^2*Kernel.normal(x)}
integrand.K.bias2 = function(x){(-(.75)^2*sin(.75*x))^2}
integrand.K.bias2.nw = function(x){(-(.75)^2*sin(.75*x)+2*.75*cos(.75*x)*(-x/(delta)^2))^2}
h.formula.plot = function(n,del,delta.x){(1/n)^(1/5)*(integrate(integrand.K.var1,-Inf,Inf)$value*integrate(integrand.K.var2,-2*del,2*del,deltax=delta.x)$value/((integrate(integrand.K.bias1,-Inf,Inf)$value)^2*integrate(integrand.K.bias2,-2*del,2*del)$value))^(1/5)}
h.formula.nw.plot = function(n,del){(1/n)^(1/5)*(integrate(integrand.K.var1,-Inf,Inf)$value*integrate(integrand.K.var2,-2*del,2*del)$value/((integrate(integrand.K.bias1,-Inf,Inf)$value)^2*integrate(integrand.K.bias2.nw,-2*del,2*del)$value))^(1/5)}

##Nadaraya-Waston##
h.nw = h.formula.nw.plot(n,delta)

for(k in 1:5){

X.r.nw = rnorm(n,0,delta)
error.nw = rnorm(n)
y.r.nw = sin(.75*X.r.nw)+.3*error.nw
m.est.nw = rep(0,length(x))

for(i in 1:length(x)){
     m.est.nw[i] = sum(exp(-((x[i]-X.r.nw)/h.nw)^2/2)*y.r.nw)/sum(exp(-((x[i]-X.r.nw)/h.nw)^2/2))
}
if(k == 1){
plot(x,m.est.nw,type = "l",ylim = c(-0.6,0.6),ylab = "y")
}
else{
  lines(x,m.est.nw,lty = k)
}
}

##Gasser-Muller##
h.gm = h.formula.plot(n,delta,.09*1.5)

for(l in 1:5){

X.r.gm = sort(rnorm(n,0,delta),decreasing = FALSE)
t.gm = c(-Inf,rep(0,n-1),Inf)
for(k in 2:n)
{
  t.gm[k] = (X.r.gm[k-1]+X.r.gm[k])/2
}
error.gm = rnorm(n)
y.r.gm = sin(.75*X.r.gm)+.3*error.gm
m.est.gm = rep(0,length(x))

for(i in 1:length(x)){
  for(j in 2:(n+1)){
  m.est.gm[i] = m.est.gm[i] + y.r.gm[j-1]*integrate(function(t){1/(h.gm*sqrt(2*pi))*exp(-((x[i]-t)/h.gm)^2/2)},t.gm[j-1],t.gm[j])$value
  }
}
if(l == 1){
  plot(x,m.est.gm,type = "l",ylim = c(-0.6,0.6),ylab = "y")
}
else{
  lines(x,m.est.gm,lty = l)
}
}

##Local linear smoother##
h.ls = h.formula.plot(n,delta,.09)

for(k in 1:5){

X.r.lls = rnorm(n,0,delta)
error.lls = rnorm(n)
y.r.lls = sin(.75*X.r.lls)+.3*error.lls
m.est.ls = rep(0,length(x))

for(i in 1:length(x)){
  s.n.1 = sum(1/sqrt(2*pi)*exp(-((x[i]-X.r.lls)/h.ls)^2/2)*(x[i]-X.r.lls))
  s.n.2 = sum(1/sqrt(2*pi)*exp(-((x[i]-X.r.lls)/h.ls)^2/2)*(x[i]-X.r.lls)^2)
  w = 1/sqrt(2*pi)*exp(-((x[i]-X.r.lls)/h.ls)^2/2)*(s.n.2-(x[i]-X.r.lls)*s.n.1)
  m.est.ls[i] = sum(w*y.r.lls)/sum(w)
}
if(k == 1){
  plot(x,m.est.ls,type = "l",ylim = c(-0.6,0.6),ylab = "y")
}
else{
  lines(x,m.est.ls,lty = k)
}
}
