## Extra stuff: Compare sample means with model for selected subareas (assumes main litter script has been run first)

###################
## Baltic exmaple
###################

subgrid = subset(bgrid,lon>13)
subdat = subset(drd,lon>13)[[2]]
StdEffort = 1e6 / nrow(subgrid) 

idx2 = redoSurveyIndex(drd,model,predD=subgrid,predfix=list(EFFORT=StdEffort))

subdat$CPUE = subdat$TotalLitter/subdat$EFFORT * mean(subdat$EFFORT)
subdat.y = split(subdat,subdat$Year)
png("output/BalticExample.png",width=1200,height=800,pointsize=16)
par(mfrow=n2mfrow(nlevels(d$Year)+6),mar=c(2,2,2,2),oma=c(2,2,2,2))

xlims = range(subdat$lon); ylims = range(subdat$lat)
lapply(subdat.y,function(x){  bubblePlot(df2dr(x),response="CPUE",xlim=xlims,ylim=ylims,rim=TRUE,scale=1/3); title(x$Year[1]) })

surveyIdxPlots(idx2,df2dr(subdat),predD=subgrid,select="absolutemap",par=list(),year=c(2012,2019,2021),legend=FALSE,mapBubbles=TRUE,colors=mycolors)

surveyIndex:::plot.SIlist(list(idx2),main="TotalLitter trend")
ta = trendAnalysis(idx2)
abline(ta$model,col=3)
ta2 = trendAnalysis(idx2,10)
abline(ta2$model,col=2)
legend("topleft",lty=1,col=c(2,3),legend=c(paste("Trend all years p=",round(ta2$pvalue,3)),
                                           paste("Trend last 6 years p=",round(ta$pvalue,3))))

## Function to calculate mean + bootstrapped confidence intervals
do.mean<-function(x, g2, y, minHauls=3, minSum=1 , nboot = 1000, fun=mean){
    tmp = subset(x[[2]], Gear2==g2 & Year==as.character(y))
    if(nrow(tmp)<minHauls) return(c(NA,NA,NA))
    med = fun(tmp$Nage[,1] / (tmp$EFFORT/1e6))
    if(sum(tmp$Nage[,1])<minSum) return(c(med,NA,NA))
    bsmed = sapply(1:nboot,function(xx){
        set.seed(xx)
        sel = sample.int(nrow(tmp),nrow(tmp),replace=TRUE)
        bsdata = tmp$Nage[sel, 1 ] / (tmp$EFFORT[sel]/1e6)
        fun(bsdata)
    })
    return(c(med,quantile(bsmed,probs=c(0.025,0.975))))
}

out=list()
for(yy in levels(subdat$Year))
    out[[yy]] = do.mean(df2dr(subdat),g2="OTT",y=yy)
out.df = do.call("rbind",out)

plot(rownames(out.df),out.df[,1],ylim=c(0,200),type="b",main="Sample mean")
lines(rownames(out.df),out.df[,2],lty=2)
lines(rownames(out.df),out.df[,3],lty=2)

dev.off()
#################################
## Repeat Celtic Sea / Biscay
#################################
subgrid = subset(bgrid,lon<0 & lat <52)
subdat = subset(drd,lon<0 & lat <52)[[2]]

StdEffort = 1e6 / nrow(subgrid) 
idx2 = redoSurveyIndex(drd,model,predD=subgrid,predfix=list(EFFORT=StdEffort))

subdat$CPUE = subdat$TotalLitter/subdat$EFFORT * mean(subdat$EFFORT)
subdat.y = split(subdat,subdat$Year)

png("output/CelticExample.png",width=1200,height=800,pointsize=16)
par(mfrow=n2mfrow(nlevels(d$Year)+7),mar=c(2,2,2,2),oma=c(2,2,2,2))

xlims = range(subdat$lon); ylims = range(subdat$lat)
lapply(subdat.y,function(x){  bubblePlot(df2dr(x),response="CPUE",xlim=xlims,ylim=ylims,rim=TRUE,scale=1/3); title(x$Year[1]) })

surveyIdxPlots(idx2,df2dr(subdat),predD=subgrid,select="absolutemap",par=list(),year=2014:2018,legend=FALSE,mapBubbles=TRUE,colors=mycolors)

surveyIndex:::plot.SIlist(list(idx2),main="TotalLitter trend")
ta = trendAnalysis(idx2)
abline(ta$model,col=3)
ta2 = trendAnalysis(idx2,10)
abline(ta2$model,col=2)
legend("topleft",lty=1,col=c(2,3),legend=c(paste("Trend all years p=",round(ta2$pvalue,3)),
                                           paste("Trend last 6 years p=",round(ta$pvalue,3))))

out=list()
for(yy in levels(subdat$Year))
    out[[yy]] = do.mean(df2dr(subdat),g2="OTT",y=yy)
out.df = do.call("rbind",out)

plot(rownames(out.df),out.df[,1],ylim=c(0,200),type="b",main="Sample mean")
lines(rownames(out.df),out.df[,2],lty=2)
lines(rownames(out.df),out.df[,3],lty=2)

dev.off()

