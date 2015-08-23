-- 18.05.2007

library("POT")
library(MASS)
library(mvtnorm)
library(scatterplot3d)
library(sn)
library(mnormt)
library(copula)
a2_load3 <- function (p_path) {
		source(paste(p_path,"x1.r",sep=""))
		source(paste(p_path,"x2.r",sep=""))
		source(paste(p_path,"x3.r",sep=""))		
		return(cbind(X1,X2,X3))
}
a2_mult <- function (p_m, p_dim) {
		for (j in 1:p_dim) {
				for (i in 1:(length(X)/p_dim)) {
					X[i,j] = X[i,j]*p_m
				}
		}
		return (X)
}
a2_pnorm_w_point <- function(x0, F0, mu, precision) 
{
 step = 1
 sigma = 1
 dif = pnorm(c(x0),mu,sigma ) - F0 
 while (abs(dif) > precision) {	
	if (dif > precision)
		if (step < 0) step = -step/2 else step = step
	else if (dif < -precision)
		if (step > 0) step = -step/2 else step = step
	else
		1=1;	
	sigma = sigma + step
	dif = pnorm(c(x0),mu,sigma ) - F0 
 }
 return(sigma)
}
ptail_gpd <- function(t, p_X, p_mu, p_j, p_threshold_pct)
{
return ((1-p_threshold_pct)*pgpd(c(t-quantile(p_X[,p_j], p_threshold_pct)), p_mu[p_j], mle[,p_j]$fitted.values[1], mle[,p_j]$fitted.values[2]) + p_threshold_pct)
}
pcomb_distr <- function(t, p_X, p_mu, p_j, p_sigma, p_threshold_pct)
{
 if(t > quantile(p_X[,p_j], p_threshold_pct)) return (ptail_gpd(t, p_X, p_mu, p_j, p_threshold_pct))
 else return (pnorm(c(t), p_mu[p_j], p_sigma[p_j]))
}
a2_pjoint <- function(t, p_X, p_mu, p_sigma, p_threshold_pct ) {
		  pcopula(a2_cop, c(pcomb_distr(t[1], p_X, p_mu, 1, p_sigma, p_threshold_pct),pcomb_distr(t[2], p_X, p_mu, 2, p_sigma, p_threshold_pct),pcomb_distr(t[3], p_X, p_mu, 3, p_sigma, p_threshold_pct)))
}
inv_pcomb_distr <- function(t_f, p_X, p_mu, p_j, p_sigma, p_threshold_pct, precision) 
{
 step = 1
 v_x = 0
 dif = pcomb_distr(v_x, p_X, p_mu, p_j, p_sigma, p_threshold_pct) - t_f
 while (abs(dif) > precision) {	
	if (dif > precision)
		if (step > 0) step = -step/2 else step = step
	else if (dif < -precision)
		if (step < 0) step = -step/2 else step = step
	else
		step = step;	
	v_x = v_x + step
	dif = pcomb_distr(v_x, p_X, p_mu, p_j, p_sigma, p_threshold_pct) - t_f
 }
 return(v_x)
}
