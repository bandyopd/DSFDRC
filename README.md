# ESFDRC
R package for two dimensional marginal screening for the ultra-high dimensional survival data.

Head and neck cancer is the 6th most common cancer worldwide with an expected 1.08 million new cases each year. 
Such cancer data are ultra-high dimensional with thousands of clinical features and gene expressions, 
making it challenging for the traditional analytical tools to extract the potential biomarker for the cancer survival and control false discoveries.
In addition, presence of heavy censoring can affect the screening procedures based on Kaplan-Meier (K-M) survival estimates. 
In this package, we develope a model-free screening procedure for the ultra-high dimensional right censored survival data without the involvement 
of the K-M estimates. Our method can capture both linear and non-linear correlation between the predictor and outcome with the help of **ECCFIC**-correlation proposed by Ke and Yin (2019) and 
does not require any complex estimation or optimization. Furthermore, we tackled the issue of false discoveries in the screening procedure by adapting the 
**KnockOff** methods by Candes et al. (2018). 

## Description

A nice vignettee demonstrates the example of [TCGA-HNSC](https://portal.gdc.cancer.gov/projects/TCGA-HNSC) data is available 
here: [Vignette](http://htmlpreview.github.io/?https://github.com/urmiaf/ESFDRC/blob/master/vignettes/Introduction.html)

## Usage

ECCFIC_screen(time, delta, x_mat, kernel = "gaussian")  #screening function

ESFDRC_func(      #screening with FDR control

  data_n1,
  
  data_n2,
  
  rand_num,
  
  q,
  
  s = round(nrow(data_n1)/log(nrow(data_n1)), 0)
)

#where,

#time	=a numeric vector of survival time
#delta	=a numeric vector of censoring indicator
#x_mat =a matrix/dataframe of continious covariates
#kernel	=a kernel to use for x, 'gaussian' or 'distance',default gaussian. 

  #data_n1=a data set with column 1= 'time',column2='delta' and rest are 'covariates' (can be obtained from the function split_data())
  
  #data_n2=a data set with column 1= 'time',column2='delta' and rest are 'covariates' Note: n1<n2
                                                                                                                                                   
  #rand_num=a random seed to reproduce the result
                                                                                                                                                   
  #q=a prespecified false discovery rate (usually .05 or .10)                                                                                                                                                  
 #s =the number of covariates to be screened in the 1st step. Default is (n/log(n)) where 'n' i the number of rows in 'data_n1                                                                                                                               


                                                                                                                                                   
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
