Import.Data<-function(path='Data/EMD_seafloorbaselines_EEA_waters_filtered_20230317.csv',type='EMOD'){
  
  #Import Data  
  data<-read.csv(path)
  
  if (type == 'EMOD'){
    #isolate the variables that will be used for assessment
    small_data<-data%>%select(country_waters, SurveyName, Gear, hal_id, Date, LTREF, PARAM, 
                LTSZC, LT_Items, UnitItem, LT_Weigth, UnitWgt)
    #define parameter names for further use
    colnames(small_data)<-c('Country', 'Survey', 'Gear', 'HaulID', 'Date', 'LTREF', 'PARAM', 
                          'LTSZC', 'LT_Items', 'UnitItem', 'LT_Weigth', 'UnitWgt')
    #convert string to numeric value
    small_data$LT_Items <-as.numeric(small_data$LT_Items)
    #Isolate the Year
    small_data$Year<-lubridate::year(small_data$Date)
    #Data_Haul allows to find the unique number of hauls. A haul containing different litter items will have the same name.
    small_data$Date_Haul<-paste0( small_data$Date, small_data$HaulID)
    #return dataset
    return(small_data)
  }else if(type %in% c('BITS', 'BTS', 'DYFS', 'EVHOE', 'IE-IGFS', 'NS-IBTS', 'SCOWCGFS')){
    
    dim1<-which(data$RecordType=='RecordType')
    
    #Isolate litter data from the dataset
    cnames<-data[dim1,]
    dim2<-min(which(is.na(cnames)))-1
    data<-data[seq(dim1+1, dim(data)[1]), seq(dim2)]
    colnames(data)<-cnames[seq(dim2)]
    
    #isolate the variables that will be used for assessment
    small_data<-data%>%select('Country', 'Survey', 'Gear', 'StNo', 'LTREF', 'PARAM', 'LTSZC', 
                              'LT_Items', 'UnitItem', 'LT_Weight', 'UnitWgt', 'Year','HaulNo')
    #convert string to numeric value
    small_data$LT_Items <-as.numeric(small_data$LT_Items)
    #Data_Haul allows to find the unique number of hauls. A haul containing different litter items will have the same name.  
    colnames(small_data)<-c('Country', 'Survey', 'Gear', 'StNo', 'LTREF', 'PARAM', 
                            'LTSZC', 'LT_Items', 'UnitItem', 'LT_Weigth', 'UnitWgt', 'Year', 'Date_Haul')  
    #return dataset
    return(small_data)
  }else{
    stop('Chosen type is not supported.')
  }
}
