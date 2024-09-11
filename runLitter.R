############################################################################################
## Example script for calculating standardized indices of litter abundance
## using EMODnet litter data.
## Author: Casper W. Berg, DTU Aqua.
## Sep 2024
#############################################################################################

library(DATRAS)
library(maps); library(mapdata)
library(surveyIndex)
library(marmap)


## data frame to DATRASraw object - surveyIndex package expects this structure
df2dr<-function(x){
    x$haul.id = 1:nrow(x)
    dd = list()
    dd[[1]] = data.frame()
    dd[[2]] = x
    dd[[3]] = data.frame()
    class(dd)<-"DATRASraw"
    dd
}


datafiles = list.files("data/EMODnet",pattern="*.csv$",full.names=TRUE,recursive=TRUE)

dat <- read.csv(datafiles[1])
dat$LT_Items = as.numeric(dat$LT_Items)

## Turn these columns into factors
factorcols = c("country_waters","subregion","region","ite_id","Survey_id","hal_id","SurveyName",
               "Project_id","sur_policy","Ship","Gear","Country","Originator","Collator",
               "CoordRefSys","LTREF","PARAM","LTSZC","LTSRC","TYPPL","LTPRP","UnitWgt",
               "UnitItem")
for(cc in factorcols) dat[,cc] <- factor(dat[,cc])

## Assume zero count when LT_ITEMS is NA (Important to avoid dropping hauls with zero litter)
## However, note danger that for some hauls numbers may not be recorded (only weight)
## in which case they are NOT zero hauls, and should rather be left out of the analysis...
dat$LT_Items[is.na(dat$LT_Items)] = 0

## keep only row line per haul and merge with total numbers
d <- dat[!duplicated(dat$hal_id),!grepl("LT",names(dat))]
totallitter <- aggregate(LT_Items ~ hal_id,data=dat,FUN=sum)
colnames(totallitter)[2]<-"TotalLitter"
dim(d)
d <- merge(d,totallitter,by="hal_id")
dim(d)

## Add/rename some variables (SurveyIndex package assumes some specific column names)
d$lon = d$HaulLong
d$lat = d$HaulLat
d$Year = sapply(d$Date,substr,start=0,stop=4)
d$Year = factor(d$Year,levels=as.character(sort(as.numeric(unique(d$Year)))))
d$Month = as.numeric(sapply(d$Date,substr,start=6,stop=7))
d$Day = as.numeric(sapply(d$Date,substr,start=9,stop=10))
d$timeOfYear = (d$Month-1)/12 + (d$Day-1)/365
d$DoorSpread = as.numeric(d$DoorSpread)
d$WingSpread = as.numeric(d$WingSpread)
d$GroundSpeed = as.numeric(d$GroundSpeed)
d$Distance = as.numeric(d$Distance) 
d$TimeShotHour = 12 ## fake, not used

## Simplify into two gear categories: Otter trawls and beam trawls.
d$Gear2 = "OTT"
d$Gear2[ grep("^BT",d$Gear) ] <- "BT"
d$Gear2 = factor(d$Gear2)


## Subsetting:
## In this example we only consider Northern part of data set - Medits data should probably be analyzed separately
d <- subset(d, !Gear %in% c("GOC73","BOT","BMT","NCT"))
d = subset(d, HaulDur < 120) ## exclude weirdly long hauls 
d = subset(d,! Year %in% as.character(2006:2011)) ## Too few data from 2011 and earlier to do anything meaningful
d$Year = factor(d$Year)


addSweptAreaSimple<-function(d,minSpeed=1,minDist=500,maxDistDev=0.2){
    d$GroundSpeed[ d$GroundSpeed < minSpeed ] <- NA
    d$Distance[ d$GroundSpeed < minDist ] <- NA
    d$WingSpread[ d$WingSpread<=0 ] <- NA
    d$WingSpread[ d$Gear=="GOV" & ( d$WingSpread<5 | d$WingSpread >40 ) ] <- NA

    ## Remove some distance outliers
    Distance2 <- (d$HaulDur / 60 * 1852 * d$GroundSpeed)
    badDist <- abs((d$Distance - Distance2)/d$Distance)>maxDistDev
    d$Distance[ badDist ] <- NA

    ## Impute missing wing spreads (median within gear)
    wtab <- aggregate(WingSpread~Gear,data=d,FUN=median,na.rm=TRUE)
    for(i in 1:nrow(wtab)){
        sel <- which(d$Gear==wtab$Gear[i] & is.na(d$WingSpread))
        d$WingSpread[ sel ] <- wtab$WingSpread[i]
    }


    d$WingSpread[d$Gear == "BT8"]    <- 8
    d$WingSpread[d$Gear == "BT4A"]   <- 4
    d$WingSpread[d$Gear == "BT7"]    <- 7
    d$WingSpread[d$Gear == "BT4AI"]  <- 4
    d$WingSpread[d$Gear == "BT4S"]   <- 4
    d$WingSpread[d$Gear == "BT4P"]   <- 4
    d$WingSpread[d$Gear == "BT6"]    <- 6
    d$WingSpread[d$Gear == "BT3"]    <- 3
    
    ## For Beam trawls, GearEx == "DB" means double beam, i.e. catches from two beam trawls added.
    ## But GearEx variable is not included, only one "DB" value in database from 2014.
    ## Ignore for now...
    ##d$BeamWidth[ !is.na(d$GearEx) & d$GearEx=="DB" ] <- d$BeamWidth[ !is.na(d$GearEx) & d$GearEx=="DB" ]*2

    ## There are probably errors in recorded in wing spread and unclear how CPUE correlates with wing spread. It might therefore be a more robust assumption to use a fixed (median) Wingspread for each gear type. 
    d$WingSpread.median <- d$WingSpread
    for(i in 1:nrow(wtab)){
        sel <- which(d$Gear==wtab$Gear[i])
        d$WingSpread.median[ sel ] <- wtab$WingSpread[i]   
    }

    
    ## Impute missing ground speeds
    stab <- aggregate(GroundSpeed~Gear,data=d,FUN=median,na.rm=TRUE)
    for(i in 1:nrow(stab)){
        if(!is.na(stab$GroundSpeed[i])){
            d$GroundSpeed[ is.na(d$GroundSpeed) & d$Gear==stab$Gear[i]  ] <- stab$GroundSpeed[i]
        }
    }
    for(i in 1:nrow(stab)){
        sel <- which(d$Gear==stab$Gear[i] & is.na(d$GroundSpeed))
        d$GroundSpeed[ sel ] <- stab$GroundSpeed[i]
    }

    ## Impute missing distances
    noDist <- which( is.na(d$Distance) )
    d$Distance[ noDist ] <- (d$HaulDur[noDist] / 60 * 1852 * d$GroundSpeed[noDist])
       
    d$SweptArea <- d$WingSpread * d$Distance
    d$SweptArea.median <- d$WingSpread.median * d$Distance
    d
}

d = addSweptAreaSimple(d)

d$EFFORT = d$SweptArea.median 

## Get prediction grid
bgrid <- getBathyGrid(df2dr(d),minDepth=1,maxDepth=1000,maxDist=0.3,resolution=10,shapefile="shapefiles/ICES/ICES_areas.shp",select=c("ICES_SUB"))
bgrid$Gear2 = "OTT" ## Use "OTT" as reference gear

###############################
## Add EEZ from shapefile
###############################

## name clash with EEZ column country - rename
d$Country.orig = d$Country
d$Country = NULL
d = addSpatialData(df2dr(d),shape="shapefiles/EEZshape/EMODnet_HA_OtherManagementAreas_EEZ_v11_20210506.shp")[[2]]
d$Territory = factor(d$Territory,levels=unique(d$Territory)) ## drop empty factor levels

## Add to grid
tmp = addSpatialData(df2dr(bgrid),shape="shapefiles/EEZshape/EMODnet_HA_OtherManagementAreas_EEZ_v11_20210506.shp")
bgrid = tmp[[2]]
bgrid = subset(bgrid,Territory %in% unique(d$Territory) ) ## Only include territories with data
bgrid = subset(bgrid,!Territory %in% c("Portugal","United Kingdom")) ## Drop UK and Portugal 
bgrid$Territory = factor(bgrid$Territory,levels=unique(d$Territory)) ## drop empty factor levels


## Plot grid
my.palette<-colorRampPalette(c("darkblue","mediumblue","lightblue1"))
my.palette.vec=my.palette(100);
png("output/bathygrid.png",width=1200,height=800)
plot(bgrid$lon,bgrid$lat,col=rev(my.palette.vec)[cut(bgrid$Depth,100)],pch=15,cex=1.3)
maps::map("worldHires", fill = TRUE, plot = TRUE, 
                add = TRUE, col = grey(0.5))
points(d$lon,d$lat,col=2,pch=".",cex=3)
dev.off()

## Plot EEZ map
png("output/EEZmap.png",width=1200,height=800)
eezcols = rainbow(nlevels(bgrid$Territory))
eezcols[1] = "#000000"
plot(bgrid$lon,bgrid$lat,col=eezcols[bgrid$Territory],pch=15,cex=1.3)
maps::map("worldHires", fill = TRUE, plot = TRUE, 
                add = TRUE, col = grey(0.9))
legend("bottomright",legend=levels(bgrid$Territory),col=eezcols,pch=15,bg="white",cex=1.3)
dev.off()


### Standardized effort: 1 km^2 (1e6 m^2).
## Index is sum of grid points, so divide by number of grid points, such that the index (sum) is
## the mean litter abundance.
## NB : Remember to update this if grid is redefined!
StdEffort = 1e6 / nrow(bgrid) 

##########################################
## Define and fit model
##########################################

formul = "Year + Gear2 + s(lon,lat, bs='ds',m=c(1,0.5),k=512) + s(lon,lat,bs=c('ds'),k=64,m=c(1,0.5),by=Year,id=1) + offset(log(EFFORT))"

## Data needs to be in DATRASraw format, and with a matrix called "Nage" for the response:
drd = df2dr(d)
drd$Nage = matrix(d[,"TotalLitter"], nrow=nrow(d),ncol=1)
colnames(drd$Nage)<-1

system.time( model  <- getSurveyIdx(drd,ages=1,predD=bgrid,cutOff=0,fam="negbin",mc.cores=1,modelP=formul,nBoot=1000, predfix=list(EFFORT=StdEffort),control=list(trace=TRUE,maxit=20,nthreads=2) ) )

summary(model$pModels[[1]])
AIC.surveyIdx(model)


mycolors = c(rev(heat.colors(6)),"darkred")
mycolors.alt = rev(c("red","orange","yellow","green"))

lt="Total litter"

mycex = 0.4
png("output/maps.png",width=1200,height=1000,pointsize=24)
surveyIdxPlots(model,drd,myids=NULL,predD=bgrid,select="absolutemap",colors=mycolors,legend=TRUE,legend.signif=2,map.cex=mycex,par=list(mfrow=n2mfrow(nlevels(drd$Year)),mar=c(0,0,2,0)),main=lt,year=levels(drd$Year),legend.pos="bottomright")
##mapLegend(model,levels(drd$Year),mycolors) ## This can be used to get legend in separate plot 
dev.off()

png("output/maps4col.png",width=1200,height=1000,pointsize=24)
surveyIdxPlots(model,drd,myids=NULL,predD=bgrid,select="absolutemap",colors=mycolors.alt,legend=TRUE,legend.signif=2,map.cex=mycex,par=list(mfrow=n2mfrow(nlevels(drd$Year)),mar=c(0,0,2,0)),main=lt,year=levels(drd$Year),legend.pos="bottomright")
dev.off()

png("output/mapsBubbles.png",width=1200,height=1000,pointsize=24)
surveyIdxPlots(model,drd,myids=NULL,predD=bgrid,select="absolutemap",colors=mycolors,legend=TRUE,legend.signif=2,map.cex=mycex,par=list(mfrow=n2mfrow(nlevels(drd$Year)),mar=c(0,0,2,0)),main=lt,year=levels(drd$Year),legend.pos="bottomright",mapBubbles=TRUE)
dev.off()

## CV maps of estimated litter abundance
surveyIdxPlots(model,drd,myids=NULL,predD=bgrid,select="CVmap",year=levels(drd$Year),colors=cm.colors(6),par=list(mfrow=n2mfrow(nlevels(drd$Year)),mar=c(0,0,2,0),cex=0.6),legend=TRUE,legend.signif=2,map.cex=0.7,cutp=c(0,0.1,0.25,0.5,0.75,1,Inf))


pdf("output/overalltrend.pdf",width=8,height=6)
par(mfrow=c(1,1),mar=c(4,4,4,4))
surveyIndex:::plot.SIlist(list(model),main="Total litter mean trend (numbers/km^2 for OTT)")
dev.off()

#########################################################################
## Calculate litter trends for subareas using the grand model
#########################################################################
subgrids = list()
subidx = list()


tertab <- table(d$Territory)
goodEEZs <- names(tertab)[tertab>200] ## Only EEZs with at least 200 hauls
for(eez in goodEEZs){
    subgrids[[ eez ]] <- subset(bgrid,Territory==eez)
}

## Calculate indices and uncertainties by subgrid
for(i in 1:length(subgrids)){
    stdefi = 1e6 / nrow(subgrids[[i]])
    subidx[[ names(subgrids)[i] ]] <- redoSurveyIndex(drd,model,predD=subgrids[[i]],predfix=list(EFFORT=stdefi))
}

pdf("output/trends.pdf",width=8,height=6)
par(mfrow=c(1,1))
surveyIndex:::plot.SIlist( subidx , rescale=TRUE,allCI=TRUE,main="Litter trend by subarea")
surveyIndex:::plot.SIlist( subidx , rescale=FALSE,allCI=TRUE,main="Litter density by subarea")
dev.off()

sink("output/modelSummary.txt")
print(summary(model$pModels[[1]]))
sink()

## QQ-plot
resid = residuals(model)
png("output/qqplot.png",width=1200,height=800,pointsize=24)
qqnorm(resid)
abline(0,1,col=2,lwd=2)
dev.off()

trendAnalysis<-function(x,nyears=6){
    tmp = data.frame(litter = x$idx[,1],
                     Year = as.numeric(as.character(rownames(x$idx))),
                     sig2 = ((x$up[,1] - x$lo[,1])/4)^2)
    if(nyears>nrow(tmp)) stop("you asked for too many years")
    tm = lm(litter ~ Year, weights = tail(1/tmp$sig2,nyears), data=tail(tmp,nyears))
    summa <- summary(tm)
    pval = summa$coefficients[2,4]
    trend = summa$coefficients[2,1]
    sigma = summa$coefficients[2,2]
    citrend = confint(tm)[2,] 
    return( list( model = tm, pvalue = pval, trend = trend, sigma = sigma,citrend = citrend))
}


## as abline, but clip to xlim
myabline <- function(x,xlim,...){
    usr <- par("usr")
    clip(min(xlim),max(xlim),usr[3],usr[4])
    abline(x,...)
    clip(usr[1],usr[2],usr[3],usr[4])
}

## Global trend
trendAnalysis(model)
trendAnalysis(model,nyears=10)

trends <- list()
pdf("output/trendsbyEEZ.pdf",width=10,height=8,pointsize=10)
par(mfrow=n2mfrow(length(subidx)),mar=c(4,3,4,1))
for(i in 1:length(subidx)){
    surveyIndex:::plot.SIlist(list(subidx[[i]]),main=names(subidx)[i])
    ta = trendAnalysis(subidx[[i]])
    trends[[ names(subidx)[i] ]] <- ta 
    yrange <- range(as.numeric( tail( rownames(subidx[[i]]$idx),6)))
    myabline(ta$model,col=3,lwd=2,xlim=yrange)
    
    ta2 = trendAnalysis(subidx[[i]],10)
    yrange <- range(as.numeric( tail( rownames(subidx[[i]]$idx),10)))
    myabline(ta2$model,col=4,lwd=2,xlim=yrange)
    
    legend("topleft",lty=1,col=c(4,3),legend=c(paste("Trend last 10 years",round(ta2$trend,1),"(+/-",round(ta2$sigma*2,1),") p=",round(ta2$pvalue,3)),
                                               paste("Trend last 6 years",round(ta$trend,1),"(+/-",round(diff(ta$citrend)/2,1),") p=",round(ta$pvalue,3))),lwd=2)
}
dev.off()

