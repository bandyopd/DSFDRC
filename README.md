# DSFDRC
An R package for two dimensional marginal screening for the ultra-high dimensional survival data controlling the false discoveries.

Head and neck cancer is the 6th most common cancer worldwide with an expected 1.08 million new cases each year. 
Such cancer data are ultra-high dimensional with thousands of clinical features and gene expressions, 
making it challenging for the traditional analytical tools to extract the potential biomarker for the cancer survival and control false discoveries.
In addition, presence of heavy censoring can affect the screening procedures based on Kaplan-Meier (K-M) survival estimates. 
In this package, we develope a model-free dual screening procedure for the ultra-high dimensional right censored survival data without the involvement 
of the K-M estimates. Our method can capture both linear and non-linear correlation between the predictor and outcome with the help of **ECCFIC**-correlation proposed by Ke and Yin (2019) and 
does not require any complex estimation or optimization. Furthermore, we tackled the issue of false discoveries in the screening procedure by adapting the 
**KnockOff** methods by Candes et al. (2018). 

# Description

A nice vignettee demonstrates the example of [TCGA-HNSC](https://portal.gdc.cancer.gov/projects/TCGA-HNSC) data is available 
here: [Vignette](http://htmlpreview.github.io/?https://github.com/urmiaf/DSFDRC/blob/master/vignettes/Introduction.html)

# Installation
install.packages("devtools")

devtools::install_github("urmiaf/DSFDRC")

# Usage
### Data preparation

ar1_cov <- function(p, rho) { #AR(1) covariance structure
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - 
                    (1:p - 1))
  rho^exponent
}

rand_num<- floor(runif(1, min=100, max=9999))    #generating random number

set.seed(rand_num) 

p=5000                       #number of covariate

n=200                           #sample size

d=round(n/log(n),0)              # model size

rho=.5                           #auto regressive correlation rho

sigma2<-ar1_cov(p,rho)             #generating AR(1) cov

x<-mvrnorm(n,mu=rep(0,p),sigma2)                     #covariates from multivariate normal:mean=0,var=1,cov=k

e<- rnorm(n,0,1)

t=exp(x[,1]+.8*x[,2]+2*x[,7]^2+e)

ul<-6.4

c=rexp(n,1/(ul*exp(x[,1])))                         #depending censoring with 50% censoring rate.

delta<-as.numeric(t<=c)

censoring_prop<-round(1-mean(delta),2)                #true censoring rate in the data

censoring_prop

y=pmin(t,c)

#x,y,delta created

data<-data.frame(y,delta,x)

### Dual-Screnning: KIDS

library(DSFDRC)

screen_result<-KIDS(x,y,delta,swap=F)

### Screening with FDR control: a_KIDS

library(knockoff)

fdr=c(.1,.15,.2,.25,.3)  #desired false discovery rate

akids_result<-aKIDS(x,y,delta,n1=300,d=100,fdr=fdr,rand_number=123,swap=F)

                                                                                                                                                   
# References
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
