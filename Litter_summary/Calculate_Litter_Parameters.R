# ---------------------------------------------
#Define and install libraries
# ---------------------------------------------
list.of.packages <-c('dplyr', 'lubridate', 'xlsx')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# ---------------------------------------------
#load libraries
# ---------------------------------------------
lapply(list.of.packages, require, character.only = TRUE)

# ---------------------------------------------
#load funtions
# ---------------------------------------------
source('Litter_summary/Functions/Import_Data.R')
source('Litter_summary/Functions/Make_Matrix.R')
source('Litter_summary/Functions/Translation_LitterCat.R')

# ---------------------------------------------
#Import Data
# ---------------------------------------------
Data<-Import.Data(path,type)
#translate other cat keys form other refs
Data<-Translation.LitterCat(Data)  
#remove litteritems == -9 --> no count data available
Data<-Data%>%filter( LT_Items !=-9 )


# ---------------------------------------------
#Make summary table
# ---------------------------------------------

SummaryTable<-Data%>%
  select(Country, Year, Gear, Date_Haul)%>%
  group_by(Country, Year, Gear)%>%
  dplyr::summarise(NumberOfTracks=length(unique(Date_Haul)))

# ---------------------------------------------
#Make Category Table
# ---------------------------------------------
#isolate data with litter counts
data<-Data%>%filter(LTREF!='RECO-LT')%>%filter(LTREF!='RECO_LT')

data$LitterCat<-substr(data$PARAM,1,1)
#Summarize info in matrix
CatTable<-Make.matrix(data, type='cat')

# ---------------------------------------------
#Make Subcategory Table
# ---------------------------------------------
data<-Data%>%filter(LTREF!='RECO-LT')%>%filter(LTREF!='RECO_LT')
#Summarize info in matrix
SubcatTable<-Make.matrix(data, type='sub')

# ---------------------------------------------
#Make Excel
# ---------------------------------------------
write.xlsx(as.data.frame(SummaryTable), file=paste0('output/Summary_',type,  '_Data.xlsx'), sheetName = "Litter_summary", 
           col.names = TRUE, row.names = F, append = FALSE)
write.xlsx(as.data.frame(CatTable), file=paste0('output/Summary_',type,  '_Data.xlsx'), sheetName = "Litter_Cat_summary", 
           col.names = TRUE, row.names = TRUE, append = TRUE)
write.xlsx(SubcatTable, file=paste0('output/Summary_',type,  '_Data.xlsx'), sheetName = "Litter_Subcat_summary", 
           col.names = TRUE, row.names = TRUE, append = TRUE)