library(PCDSpline)
library(HI)
load("C:/Users/wangl/Dropbox/project_lin/realdata/skinTumor.rda")

n<-max(skinTumor$id)
#record the number of observations for all patients
k<-as.numeric(table(skinTumor$id))
K<-max(k)
tz<-matrix(,n,2*K)
x1<-c();x2<-c();x3<-c();x4<-c();
for (r in 1:n){
  rownum<-which(skinTumor$id==r)
  #record all observation times
  tz[r,1:(2*k[r])]<-c(skinTumor[rownum,]$time, skinTumor[rownum,]$count)
  x1[r]<-skinTumor[which(skinTumor$id==r),]$age[1]
  x2[r]<-skinTumor[which(skinTumor$id==r),]$male[1]
  x3[r]<-skinTumor[which(skinTumor$id==r),]$dfmo[1];
  x4[r]<-skinTumor[which(skinTumor$id==r),]$priorTumor[1]
}
X<-cbind((x1-mean(x1))/sd(x1),x2,x3,(x4-mean(x4))/sd(x4))

G<-cbind(X,tz)

##############################################################################
# Ispline and arms preparation 
##############################################################################
#Ispline function; n(spline function) = n(knots) + order -2; return results of n spline function for each value in x.
Ispline<-function(x,order,knots){
  # M Spline function with order k=order+1. or I spline with order
  # x is a row vector of recurrent time
  # k is the order of I spline
  # knots are a sequence of increasing points
  # the number of free parameters in M spline is the length of knots plus 1.
  k=order+1
  m=length(knots)
  n=m-2+k # number of parameters
  t=c(rep(1,k)*knots[1], knots[2:(m-1)], rep(1,k)*knots[m]) # newknots
  yy1=array(rep(0,(n+k-1)*length(x)),dim=c(n+k-1, length(x)))
  for (l in k:n){
    yy1[l,]=(x>=t[l] & x<t[l+1])/(t[l+1]-t[l])
  }
  
  yytem1=yy1
  for (ii in 1:order){
    yytem2=array(rep(0,(n+k-1-ii)*length(x)),dim=c(n+k-1-ii, length(x)))
    for (i in (k-ii):n){
      yytem2[i,]=(ii+1)*((x-t[i])*yytem1[i,]+(t[i+ii+1]-x)*yytem1[i+1,])/(t[i+ii+1]-t[i])/ii
    }
    yytem1=yytem2
  }
  index=rep(0,length(x))
  for (i in 1:length(x)){
    index[i]=sum(t<=x[i])
  }
  yy=array(rep(0,(n-1)*length(x)),dim=c(n-1,length(x)))
  if (order==1){
    for (i in 2:n){
      yy[i-1,]=(i<index-order+1)+(i==index)*(t[i+order+1]-t[i])*yytem2[i,]/(order+1)
    }
  }else{
    for (j in 1:length(x)){
      for (i in 2:n){
        if (i<(index[j]-order+1)){
          yy[i-1,j]=1
        }else if ((i<=index[j]) && (i>=(index[j]-order+1))){
          yy[i-1,j]=(t[(i+order+1):(index[j]+order+1)]-t[i:index[j]])%*%yytem2[i:index[j],j]/(order+1)
        }else{
          yy[i-1,j]=0
        }
      }
    }
  }
  return(yy)
}


#prepare for sampling beta (arms)
beta_fun=function(x, j, beta,  zi, rlj) 
{
  beta[j]=x
  tt = sum(xcov[,j]*x*zi)-sum(w*exp(xcov%*%beta)*rlj)-x^2/(2*sig0^2) # prior(beta)=N(0,sig0^2) 
  return(tt)     
}
beta_ind_fun=function(x,j,beta, zi, rlj) (x > -10)*(x < 5)


#prepare for sampling theta (arms)
#theta0 is the parameter for G0.
theta_fun = function(x, w, nh, theta0){
  nh*(x*log(x) - log(gamma(x))) +(x - 1)*sum(log(w)) - x*sum(w) - theta0*x + (theta0 -1)*log(x)
}
theta_ind_fun = function(x,w,nh, theta0)(x > 1e-4)*(x<70)

#prepare for sampling eta
eta_fun = function(x)
{
  lfc = (x-1)*sum(log(w)) - x*(sum(w)+b_eta) + (n*x+a_eta-1)*log(x) -n*log(gamma(x))
  return(lfc)
}
eta_ind_fun = function(x) (x>0)*(x<70)


###############################################################################
#dirichlet simulation
###############################################################################
total=10000
burnin=3000

#the number of subjects
n=dim(G)[1]
#the number of knots and order
klength = 20
order = 3
#the number of spline functions
kk = klength-2+order
#the number of parameters(beta)
pars = 4
#the number of columns in the dataset
col_l = dim(G)[2]
col = (col_l - pars)/2


####################################################
#set initial values and store vectors
######################################################
# observed times for each subject Ji
Ji=(rowSums(!is.na(G)) - pars)/2

# x covariates
xcov=G[, 1:pars]

#observed time
TT= c()
for (i in 1:n){
  temp = G[i,(pars+1):(Ji[i]+pars)]
  TT=unlist(c(TT,temp))
}

# N(t) of two adjacent times for the first event
#ZZ11 is the total recurrent times of event 1 for each subject.
Z1 = c()
ZZ11 = c()
for (i in 1:n){
  temp = G[i,(pars+1+Ji[i]):(pars + 2*Ji[i])]
  temp2 = sum(temp)
  Z1=unlist(c(Z1,temp))
  ZZ11 = c(ZZ11, temp2)
}


# the range of time considered    
Range=max(TT)
knots=seq(0,Range,by=(Range-0)/(klength-1))

# results from baseline function for the observation times
bl1=Ispline(TT, order, knots)

# results from baseline function for the fixed times
seqt = seq(0.1, 1880, by = 1)
blseqt=Ispline(seqt,order,knots)

#points for test frailty density
t_s <- seq(1e-2,6,0.05)

#bl1_enk is the spline for the last value of each subject.
bl1_enk=cbind(bl1[,cumsum(Ji)])

#get all  I(ti,j) - I(ti,(j-1)); N is the number of observation times of all the subjects
N = sum(Ji)
df=array(,dim=c(kk,N))
df[,1]=0
for (j in 2:N)
{
  df[,j]=bl1[,j]-bl1[,j-1];
}
df = round(df, 20)

#df1 = I(t(i1))
sum1=1;
for (k in 1:n)
{    
  df[,sum1]=bl1[,sum1];	
  sum1=sum1+Ji[k]
}

# inititial values
rl1 = matrix(rgamma(kk, 1, 1), ncol=kk)
#beta
beta1=matrix(rep(0, pars), ncol=1)
sig0 = 100
# lambda for rl
a_lam=1
b_lam=1
lambda = rgamma(1, a_lam,b_lam)
#alpha
aa = 1
ba = 1
a = rgamma(1, aa, ba)
#Nc the number of clusters
Nc =2
#p
p <- rep(1/Nc, Nc)
#theta
theta <-rgamma(Nc, 1, 1)
#ki, the membership of each individual
ki <- t(rmultinom(n,1,p))
#nh, the counts belongs to each cluster
nh <- colSums(ki)
#w, the initial frailty parameters
v <- ki%*%theta
w <- rgamma(n,v,v)
#theta0
theta0 = 0.01

#store MCMC results
pij =  array(0,dim=c(kk,N))
parRl =  array(0,dim=c(total,kk))
parLam = array(0,dim=c(total))
parBeta =  array(0,dim=c(total,pars))
parW =  array(0,dim=c(total,n))
parDenW = array(0,dim=c(total,length(t_s)))

################################################################
#Simulation 
################################################################
for (iter in 2:total)
{ 
  pij=as.vector(rl1)*df
  for (j in 1:N)
  {
    if (sum(pij[,j])==0) {pij[,j]=10^(-10)}
  }
  
  #sample Zijl1 and Zijl2
  ZZ1 =NULL 
  for (j in 1:N)
  {
    ZZZ1=rmultinom(1,Z1[j], pij[,j])
    ZZ1=cbind(ZZ1, ZZZ1)
  }
  
  
  ##sample rl1 
  sum_zij1 = apply (ZZ1, 1, sum)
  rl1 = rgamma(kk, sum_zij1+1, bl1_enk%*%(exp(xcov%*%beta1)*w) + lambda)
  
  parRl[iter,]=rl1    
  
  ###sample lambda
  lambda = rgamma(1, kk+a_lam, b_lam+sum(rl1))
  
  parLam[iter]=lambda
  
  ###commonly used in the arms function########
  rlj1 = as.vector(rl1%*%bl1_enk)
  ########################################
  
  #adaptive rejection metropolis sampling sample beta
  for (j in 1:pars){
    beta1[j]=arms(beta1[j],beta_fun, beta_ind_fun, 1, j = j, beta=beta1,  zi = ZZ11, rlj = rlj1)
  }
  parBeta[iter,]=beta1
  
  #sample wi
  eta <- ki%*%theta
  w= rgamma(n, ZZ11 + eta, exp(xcov%*%beta1)*rlj1 + eta)   
  w = sapply(w, max, 1e-20)
  parW[iter,]=w
  
  
  #sample theta, only only update those theta that generate frailty wi.
  for (d in which(nh != 0)){
    theta[d] <- arms(theta[d], theta_fun, theta_ind_fun, 1, w[ki[,d]==1], nh[d], theta0)
  }
  
  #update pi
  ga_i <- lapply(w, dgamma, theta, theta)
  ga_i <- do.call(rbind, ga_i)
  pi <- p*t(ga_i)
  
  #update p
  # a is alpha in the dirichlet distribution.
  v_a <- 1+nh[1:Nc]
  v_b <- a + n -cumsum(nh)[1:Nc]
  V <- rbeta(Nc, v_a, v_b)
  IV <- 1-V
  p[1] <- V[1]
  for (i in 2:Nc){
    p[i] = prod(IV[1:(i-1)])*V[i]
  }
  
  #update ki 
  add_to = 0
  ui <- runif(n,rep(0,n),ki%*%p)
  while(sum(p) < 1 - min(ui)){
    Nc = Nc + 1
    V <- c(V, rbeta(1,1,a))
    IV <- 1-V
    p[1] <- V[1]
    for (i in 2:Nc){
      p[i] = prod(IV[1:(i-1)])*V[i]
    }
    theta = c(theta,rgamma(1,1,1))
    #update pi
    ga_Nc <- dgamma(w, theta[Nc], theta[Nc])
    ga_i <- cbind(ga_i,ga_Nc)
    pi <- p*t(ga_i)
    add_to = add_to+1
  }
  pi <- p*t(ga_i)
  
  ki = cbind(ki,matrix(NA,n,add_to))
  for (i in 1:n){
    ki[i,] <-rmultinom(1, 1, pi[, i])
  }
  #update nh
  nh <- colSums(ki)
  
  a <- rgamma(1,(aa+sum(V != 1)),(ba-sum(log(IV[V != 1]))))
  
  dphi <- lapply(t_s,dgamma, theta, theta)
  dphi <- do.call(rbind,dphi)
  density_phi <- colSums(p*t(dphi))
  parDenW[iter,]=density_phi 

}


estbeta = apply(parBeta[seq((burnin+1), total),], 2, mean)
estbetaq = apply(parBeta[seq((burnin+1), total),], 2, quantile, c(0.025, 0.975))
estRl = apply(parRl[seq((burnin+1), total),], 2, mean)
estRlq = apply(parRl[seq((burnin+1), total),], 2, quantile, c(0.025, 0.975))
estW = apply(parW[seq((burnin+1), total),], 2, mean)
estDenW = apply(parDenW[seq((burnin+1), total),], 2, mean)
estbeta1std1 = apply(parBeta[seq((burnin+1), total),], 2, sd)
estblseqt = estRl%*%blseqt
estblseqtq = estRlq %*%blseqt

round(cbind(estbeta, estbeta1std1,t(estbetaq)),4)

estDens_2 = read.table("e:/density.txt")
plot(t_s,c(t(as.matrix(estDens_2))),lty=3,type = "l",col="blue",xlab = expression(phi),ylab = "Estimated Frailty Density" )
points(t_s,estDenW,lty=1,type = "l")
legend(4,0.7,legend = c('DPGM-PM','GFPM'), col = c('black','blue'),lty = c(1,3),cex = 1)

plot(seqt,estblseqt,type = "l", ylab = "Estimated Mean Function",xlab = "Observation Time (Days)")

# draw the posterior distribution for beta_2 AND beta_3
#install.packages("ggplot2")
#install.packages("patchwork")
library(ggplot2)
library(patchwork)
beta_m <- data.frame(parBeta[seq((burnin+1), total),])
hist1 = ggplot(data = beta_m, aes(x=X2)) +
  geom_histogram(mapping = aes(y=after_stat(density), alpha= 0.5)) + # Adjust binwidth as needed
  geom_density() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = expression(beta[2])) # Using expression() for Greek letter

# For hist2 with x-axis label as beta_3 in Greek
hist2 = ggplot(data = beta_m, aes(x=X3)) +
  geom_histogram(mapping = aes(y=after_stat(density), alpha= 0.5)) + # Adjust binwidth as needed
  geom_density() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = expression(beta[3])) # Using expression() for Greek letter

hist1 + hist2


