#' @title TCGA head and neck cancer data with standardized mRNA gene expression and some demographic information
#'
#' @description  518 primary-solid tumor samples from the TCGA-HNSCC data with 15,887 normalized mRNA gene expression and 6 clinical features.Patients were subjected to right censoring (lost to follow up or no event until end of the study) with 57% censoring rate. 
#' Each patients were followed up from entry into the study until event/censorship.
#'
#' @format A data frame
#' \describe{
#'   \item{time}{observed survival time of the patients}
#'   \item{status}{Patient's censoring indicator (1=death,0=censored)}
#' }
#' @source <https://github.com/urmiaf/DSFDRC>
"hnsc_data"
