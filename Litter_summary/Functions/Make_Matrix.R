Make.matrix<-function(Data=data, type='cat'){
  #remove no litter records
  data<-Data%>%filter(LTREF!= 'RECO-LT')%>%filter(LTREF!= 'RECO_LT')
  
  if (type=='cat'){
    #isolate the category
    data$Cat<-substr(data$PARAM, 1,1)
    category<-sort(unique(data$Cat))
    
  }else if(type=='sub'){
    data$Cat<-data$PARAM
    category<-sort(unique(data$Cat))
    
  }else{
    stop('type not supported')
  }
  
  #create empty matrix
  m<-matrix( 0,0  ,length(category)+2)
  #countries present in dataset
  Countries<-unique(data$Country)
  #years reported in the dataset
  Years<-sort(unique(data$Year))
  #define the variables 
  colnames(m)<-c('Country', 'Year', category)
  
  #check that all values are positive
  if (any(data$LT_Items<0)){
    stop('Negative number of litter count')
  }
  
  #define counter
  i=1
  for (c in sort(unique(Countries))){
    for (y in sort(unique(data$Year))){
      
      #isolate data for country c in year y
      temp<-data%>%filter(Year==y)%>%filter(Country==c)%>%
              group_by(Country, Year,Cat)%>%
              filter(!is.na(Country)& !is.na(Year) & !is.na(Cat))%>%
              dplyr::summarise(N=sum(LT_Items, na.rm=T))
      #if there is litter reported
      if (dim(temp)[1]>0){
        #add line to matrix with country year
        m<-rbind(m, c(c,y,rep(0,length(category))))
        #fill in the numer of litter items per cat
        m[i, temp$Cat]<-temp$N
        #increase the counter
        i<-i+1
      }
    }
  }
  #return the matrix
  return(m)
}
