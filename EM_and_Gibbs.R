library(mvtnorm)
library(ggplot2)
library(Rlab)
library(MCMCpack)
library(extraDistr)
x = faithful
alpha = c(0.5,0.5)
mu1 = c(2,55)
mu2 = c(4.4,80)
sigma1 = matrix(c(0.8,7,7,70),ncol = 2)
sigma2 = matrix(c(0.8,7,7,70),ncol = 2)
n = length(x[,1])
#EM----
# mixture normal
p_y1 = (alpha[1]*dmvnorm(x,mu1,sigma1))/(alpha[1]*dmvnorm(x,mu1,sigma1)+alpha[2]*dmvnorm(x,mu2,sigma2))
p_y2 = (alpha[2]*dmvnorm(x,mu2,sigma2))/(alpha[1]*dmvnorm(x,mu1,sigma1)+alpha[2]*dmvnorm(x,mu2,sigma2))
#iter 9 times
alpha = c(mean(p_y1),mean(p_y2))
sigma1 = 0
for (i in 1:n) {
  a = 0
  a = (as.matrix(t(x[i,]-mu1))%*%as.matrix(x[i,]-mu1)*p_y1[i])/sum(p_y1)
  sigma1 = sigma1+a
}
sigma2 = 0
for (i in 1:n) {
  a = 0
  a = (as.matrix(t(x[i,]-mu2))%*%as.matrix(x[i,]-mu2)*p_y2[i])/sum(p_y2)
  sigma2 = sigma2+a
}
mu1 = as.vector(apply(x*p_y1/sum(p_y1),2,sum))
mu2 = as.vector(apply(x*p_y2/sum(p_y2),2,sum))
p_y1 = (alpha[1]*dmvnorm(x,mu1,sigma1))/(alpha[1]*dmvnorm(x,mu1,sigma1)+alpha[2]*dmvnorm(x,mu2,sigma2))
p_y2 = (alpha[2]*dmvnorm(x,mu2,sigma2))/(alpha[1]*dmvnorm(x,mu1,sigma1)+alpha[2]*dmvnorm(x,mu2,sigma2))
Q = function(x){
  a = sum(log(alpha[1]*dmvnorm(x,mu1,sigma1))*p_y1)
  b = sum(log(alpha[2]*dmvnorm(x,mu2,sigma2))*p_y2)
  return(a+b)
}
Q(x)
#cluster
group = p_y1<0.5
x$group = group
x$group = as.numeric(x$group)
x$group = as.factor(x$group)
ggplot(x)+geom_point(aes(x = eruptions,y = waiting,col = group))+
  labs(title = "EM algorithm for different clusters",hjust = 0.5)+
  theme(plot.title = element_text(hjust = 0.5))
#Gibbs----
#generate p
p_y = function(alpha,mu1,mu2,tau1,tau2){
  p1 = alpha[1]*dnorm(x$eruptions,mu1,tau1)/(alpha[1]*dnorm(x$eruptions,mu1,tau1)+alpha[2]*dnorm(x$eruptions,mu2,tau2))
  p2 = alpha[2]*dnorm(x$eruptions,mu2,tau2)/(alpha[1]*dnorm(x$eruptions,mu1,tau1)+alpha[2]*dnorm(x$eruptions,mu2,tau2))
  p = cbind(p1,p2)
  y = apply(p,1,function(p)rcat(1,p))
  return(y)
}
p_mu = function(alpha,mu1,mu2,tau1,tau2,n1,n2,xbar1,xbar2){
  mubar1 = (1.25*2+n1*tau1*xbar1)/(1.25+n1*tau1)
  mubar2 = (1.25*4.4+n2*tau2*xbar2)/(1.25+n2*tau2)
  taubar1 = 1.25+n1*tau1
  taubar2 = 1.25+n2*tau2
  mu = c(rnorm(1,mubar1,(taubar1)^-1),rnorm(1,mubar2,(taubar2)^-1))
  return(mu)
}
p_tau = function(alpha,mu1,mu2,tau1,tau2,n1,n2,x1,x2){
  alpha1 = a + n1/2
  alpha2 = a + n2/2
  beta1 = b + sum((x1[,1]-mu1)^2)/2
  beta2 = b + sum((x2[,1]-mu2)^2)/2
  tau = c((rgamma(1,alpha1,beta1))^(-0.5),(rgamma(1,alpha2,beta2))^(-0.5))
  return(tau)
}
p_alpha = function(n1,n2){
  alpha = rdirichlet(1,c(a1+n1,a2+n2))
  return(alpha)
} 
updata = function(iter){
  ITER = NULL
  mu_1 = mu1;mu_2 = mu2;tau_1 = tau1;tau_2 = tau2
  for (i in 1:iter) {
    C = NULL
    y = p_y(alpha,mu_1,mu_2,tau_1,tau_2)
    sub_x = split(x,y)
    x1 = sub_x$`1`
    x2 = sub_x$`2`
    n1 = length(x1[,1])
    n2 = length(x2[,1])
    xbar1 = mean(x1[,1])
    xbar2 = mean(x2[,1])
    mu = p_mu(alpha,mu_1,mu_2,tau_1,tau_2,n1,n2,xbar1,xbar2)
    tau = p_tau(alpha,mu_1,mu_2,tau_1,tau_2,n1,n2,x1,x2)
    alpha = p_alpha(n1,n2)
    mu_1 = mu[1]
    mu_2 = mu[2]
    tau_1 = tau[1]
    tau_2 = tau[2]
    C = c(mu_1,mu_2,tau_1,tau_2,alpha)
    ITER = rbind(ITER,C)
  }
  y = y<1.5
  x$y = y
  x$y = as.numeric(x$y)
  x$y = as.factor(x$y)
  plot = ggplot(x)+geom_point(aes(x = eruptions,y = waiting,col = y))+
    labs(title = "Gibbs algorithm for different clusters",hjust = 0.5)+
    theme(plot.title = element_text(hjust = 0.5))
  return(list(plot,ITER))
}
#iter
mu1 = 2;mu2 = 4.4
a = b = a1 = a2 = 1
tau1 = 1/(sigma1[1,1])
tau2 = 1/(sigma2[1,1])
plot_iter = updata(1000)
plot_iter[[1]]
hyper_theta = plot_iter[[2]]
