rm(list=ls()) #clean up memory before start
options(dplyr.summarise.inform = FALSE) #Do not print info from dplyr package 
# ---------------------------------------------
#User parameters
# ---------------------------------------------
#path to data file
path='data/EMODnet/EMD_seafloorbaselines_EEA_waters_filtered_20230317.csv'

#which dataset do you want to import? Choose amongst EMOD, BITS, 
#BTS, DYFS, EVHOE, IE-IGFS, NS-IBTS, SCOWCGFS
type='EMOD'
CSV_separator='\t' # '\t' of tab separetor other options are ',' (comma separated), ';' (semi column), ' ' (space separated)
# ---------------------------------------------
# Make the Overview Excel
# ---------------------------------------------
source("Litter_summary/Calculate_Litter_Parameters.R")

# ---------------------------------------------
# Make the GAM models
# ---------------------------------------------
source('GAM/runLitter.R')
