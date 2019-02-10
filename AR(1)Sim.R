###Simulation AR(1)
alpha <- 0.3
e_t <- rnorm(100,0,1)
mu <- 5

y_t <- as.numeric(length(e_t))
y_t[1] <- 0

for(i in 2:length(e_t)){
  y_t[i] <- mu + alpha*y_t[i-1] + e_t[i]
}

plot.ts(y_t)

#######################################################################

updateAlpha <- function(y,mu,sig2){
  alpha <- NULL
  y1 <- y[-c(1)]
  y2 <- y[-c(length(y))]
  
  m.post <- (sum(y1*y2) - mu*sum(y2))/sum(y2^2)
  sig.post <- sig2/sum(y2^2)
  
  I <- 2
  while(I > 0){ #This truncates the support of alpha in the interval (-1,1),
                #because we are trying to model a stationary series.
    alpha <- rnorm(1,m.post,sqrt(sig.post))
    if(alpha >-1 & alpha < 1){
      I <- 0
    }
  }
  
  return(alpha)
}

updateSig2 <- function(a,b,y,alpha,mu,n){
  y1 <- y[-c(1)]
  y2 <- y[-c(length(y))]
  
  a.post <- (n-1)/2 + a
  b.post <- sum((y1 - mu - alpha*y2)^2)/2 + b
  
  return(1/rgamma(1,a.post,b.post))
}

updateMu <- function(sig2mu,alpha,y,sig2,n){
  y1 <- y[-c(1)]
  y2 <- y[-c(length(y))]
  
  sig2.post <- 1/((n-1)/sig2 + 1/sig2mu)
  m.post <- sig2.post * (sum(y1) - alpha*sum(y2))/sig2
  return(rnorm(1,m.post,sqrt(sig2.post)))
}

######################################################################

Niter <- 10000
Alpha.out <- array(NA, dim = Niter)
Sig2.out <- array(NA, dim = Niter)
Mu.out <- array(NA, dim = Niter)

Mu.out[1] <- 1
Alpha.out[1] <- alpha
Sig2.out[1] <- 1
Y <- y_t
N <- 100

for(i in 2:Niter){
  Alpha.out[i] <- updateAlpha(Y,Mu.out[i-1],Sig2.out[i-1])
  Sig2.out[i] <- updateSig2(0.01,0.01,Y,Alpha.out[i],Mu.out[i-1],N)
  Mu.out[i] <- updateMu(100,Alpha.out[i],Y,Sig2.out[i],N)
  print(i)
} 

plot(Alpha.out,type='l')
hist(Alpha.out)
mean(Alpha.out)
quantile(Alpha.out,probs = c(0.025,0.975))

plot(Sig2.out,type='l')
hist(Sig2.out)
mean(Sig2.out)
quantile(Sig2.out,probs = c(0.025,0.975))

plot(Mu.out,type='l')
hist(Mu.out)
mean(Mu.out)
quantile(Mu.out,probs = c(0.025,0.975))

#############################################################################
data(airquality, package = "datasets")
Wind = airquality$Wind  # wind speed
Temp = airquality$Temp  # air temperature
N = dim(airquality)[1]  # number of data points

#############################################################################
#Modeling AR(1)
Niter <- 20000
Alpha.out <- array(NA, dim = Niter)
Sig2.out <- array(NA, dim = Niter)
Mu.out <- array(NA, dim = Niter)

Alpha.out[1] <- 0
Mu.out[1] <- 1
Sig2.out[1] <- 1
Y <- Wind
N <- 153

for(i in 2:Niter){
  Alpha.out[i] <- updateAlpha(Y,Mu.out[i-1],Sig2.out[i-1])
  Sig2.out[i] <- updateSig2(0.01,0.01,Y,Alpha.out[i],Mu.out[i-1],N)
  Mu.out[i] <- updateMu(100,Alpha.out[i],Y,Sig2.out[i],N)
  print(i)
} 

plot(Alpha.out,type='l')
hist(Alpha.out)
mean(Alpha.out)
quantile(Alpha.out,probs = c(0.025,0.975))

plot(Sig2.out,type='l')
hist(Sig2.out)
mean(Sig2.out)
quantile(Sig2.out,probs = c(0.025,0.975))

plot(Mu.out,type='l')
hist(Mu.out)
mean(Mu.out)
quantile(Mu.out,probs = c(0.025,0.975))

#################################################################################
Nburn <- 5000
Thinning <- 10
Alpha <- Alpha.out[seq(Nburn+1,Niter,by=Thinning)]
Sig2 <- Sig2.out[seq(Nburn+1,Niter,by=Thinning)]
Mu <- Mu.out[seq(Nburn+1,Niter,by=Thinning)]

################################################################################
#Residuals
Y2 <- Wind[-N]

Residuals <- array(NA, dim = c(length(Alpha), length(Y2)))

for(i in 1:dim(Residuals)[1]){
    for(j in 1:dim(Residuals)[2]){
        Residuals[i,j] <- Mu[i] + Alpha[i] * Y2[j] + Sig2[i]*rnorm(1,0,1)
    }
}

medians <- apply(Residuals,2,median)
plot(medians,type='l')
library(tseries)
adf.test(medians)

#Prediction
Prediction <- array(NA, dim = c(length(Alpha), 10))
Prediction[,1] <- Mu + Alpha * Y[N] + Sig2*rnorm(length(Alpha),0,1)
for(i in 1:dim(Prediction)[1]){
    for(j in 2:dim(Prediction)[2]){
      Prediction[i,j] <- Mu[i] + Alpha[i] * Prediction[i,j-1] + Sig2[i]*rnorm(1,0,1)
    }
}

means <- colMeans(Prediction)
plot.ts(Y,xlim=c(0,163))
lines(x = 154:163,means,type='l',col='red')

cred1 <- NULL; cred2 <- NULL
for(i in 1:dim(Prediction)[2]){
    cred1[i] <- quantile(Prediction[,i], probs = 0.025)
    cred2[i] <- quantile(Prediction[,i], probs = 0.975)
}

