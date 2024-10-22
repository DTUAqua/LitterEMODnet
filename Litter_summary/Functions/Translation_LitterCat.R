Translation.LitterCat<-function(Data){
  #Read excel + define colnames
  litterCat<-as.data.frame(read.xlsx('Input/LitterCat.xlsx',1))
  rownames(litterCat)<-litterCat[,1]
  litterCat<-litterCat[,-1]
  
  # Check that al categories were defined in the excel
  id<-which(is.na(litterCat[Data$PARAM[Data$LTREF=='TSG-ML'],1]))
  if (length(unique( Data$PARAM[Data$LTREF=='TSG-ML'][id]))>0){
    print(paste0('Missing categories voor TSG-ML list: ', unique( Data$PARAM[Data$LTREF=='MEDITS'][id]), '. Please complete the LitterCat.xlsx.'))
  }
  
  #Change the TSG-ML to C_TS_REV
  Data$PARAM[Data$LTREF=='TSG-ML']<-
    litterCat[Data$PARAM[Data$LTREF=='TSG-ML'],1]
  
  # Check that al categories were defined in the excel
  id<-which(is.na(litterCat[Data$PARAM[Data$LTREF=='MEDITS'],1]))
  if (length(unique( Data$PARAM[Data$LTREF=='MEDITS'][id]))>0){
    print(paste0('Missing categories voor MEDITS list: ', unique( Data$PARAM[Data$LTREF=='MEDITS'][id]), '. Please complete the LitterCat.xlsx.'))
  }
  

  #Change the MEDITS to C_TS_REV
  Data$PARAM[Data$LTREF=='MEDITS']<-
    litterCat[Data$PARAM[Data$LTREF=='MEDITS'],1]

  #Delete the cats that were marked as nolitter (including total unspecified counts)
  todelete<-which(Data$PARAM=='nolitter')
  if (length(todelete)>0){
  Data<-Data[-1*todelete,]
  }
  return(Data)
  
}