# ESFDRC
An R package for two dimensional marginal screening for the ultra-high dimensional survival data.

Head and neck cancer is the 6th most common cancer worldwide with an expected 1.08 million new cases each year. 
Such cancer data are ultra-high dimensional with thousands of clinical features and gene expressions, 
making it challenging for the traditional analytical tools to extract the potential biomarker for the cancer survival and control false discoveries.
In addition, presence of heavy censoring can affect the screening procedures based on Kaplan-Meier (K-M) survival estimates. 
In this package, we develope a model-free screening procedure for the ultra-high dimensional right censored survival data without the involvement 
of the K-M estimates. Our method can capture both linear and non-linear correlation between the predictor and outcome with the help of **ECCFIC**-correlation proposed by Ke and Yin (2019) and 
does not require any complex estimation or optimization. Furthermore, we tackled the issue of false discoveries in the screening procedure by adapting the 
**KnockOff** methods by Candes et al. (2018). 

# Description

A nice vignettee demonstrates the example of [TCGA-HNSC](https://portal.gdc.cancer.gov/projects/TCGA-HNSC) data is available 
here: [Vignette](http://htmlpreview.github.io/?https://github.com/urmiaf/ESFDRC/blob/master/vignettes/Introduction.html)

# Usage
### Data preparation
p=5000             #number of covariates

n=700              #sample size

library(MASS)

rho=.5

#function for AR(1) covariance structure (pxp)

ar1_cov <- function(p, rho) {
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}


k<-ar1_cov(p,rho)                                 #generating AR(1) covariance

x<-mvrnorm(n,mu=rep(0,p),k)                     #covariates from multivariate normal

e<- rnorm(n,0,1)                                #error distribution


t=exp(5*(x[,1])+8*x[,10]^2+3*abs(x[,10])+e)     #generating survival time


delta<-as.numeric(t<=cens)                      #generating censoring time

censoring_prop<-round(1-(sum(delta)/n),2)       #true censoring rate in the data


data<-data.frame(t,cens,delta)

time<-ifelse(delta==1,t,cens)                   #observed time=min(t,cens)

data<-data.frame(time,delta,x)                  #simulated data 

### ECCFIC-Screening
library(foreach)

library(ESFDRC)

eccfic<-ECCFIC_screen(time=data$time, data$delta, x_mat=data[,-c(1:2)], kernel = "gaussian")

### ESFDRC (Screening with FDR control)
library(knockoff)

rand_num=3

two_datasets<-split_data(data, n1=100, n2=200, rand_num=rand_num)

data_n1<-two_datasets[[1]]

data_n2<-two_datasets[[2]]

final_covariates<-ESFDRC_func(
  data_n1,
  data_n2,
  rand_num,q=.05
)

                                                                                                                                                   
## References
<a id="1">[1]</a> 
Ke, C., & Yin, X. (2020). 
Expected conditional characteristic function-based measures for
testing independence. 
Journal of the American Statistical Association, 115:530, 985-996.

<a id="2">[2]</a> 
Candes, E., Fan, Y., Janson, L., & Lv, J. (2018). 
Panning for gold:‘model-x’knockoffs for high
dimensional controlled variable selection. 
Journal of the Royal Statistical Society: Series B
(Statistical Methodology), 80 (3), 551–577.
