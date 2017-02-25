setwd("D:/TEXAS")

# Load required libraries in R
packages<-function(){
  library(sp)
  library(ggplot2)
  library(foreign)
  library(stpp)
  library(dplyr)
  library(spatstat)
  library(rgdal)
  library(maptools)
  library(raster)
  library(data.table)
  library(plyr)
  library(gridExtra)
  library(spdep)
  library(INLA)
  library(goftest)
  library(spacetime)
}
packages()

##Load all tornados
TornALL <- readOGR(dsn = "./tornado/torn", layer = "torn", stringsAsFactors = FALSE)
Data_correction<-function(i=TornALL){
  i$yr <- as.integer(i$yr)
  i$mo <- as.integer(i$mo)
  i$EF <- as.integer(i$mag) 
  i$Date <- as.Date(i$date, format="%Y-%m-%d")
  i$Length <- as.numeric(i$len) * 1609.34
  i$Width <- as.numeric(i$wid) * 0.9144
  i$fat <- as.integer(i$fat)
  i$slon <- as.numeric(i$slon)
  i$slat <- as.numeric(i$slat)
  i$elon <- as.numeric(i$elon)
  i$elat <- as.numeric(i$elat)
  i$inj <- as.numeric(i$inj)
  i$Ref <- 1:nrow(i)
  return(i)
}

TornALL <- Data_correction(TornALL)
CRS.new <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
TornALL <- spTransform(TornALL, CRS.new)

#Load Boundaries
US.sp <- readOGR(dsn = "./tmp", layer = "cb_2013_us_county_5m", 
                 stringsAsFactors = FALSE)
TX.sp <- US.sp[US.sp$STATEFP == 48, ]
county <- TX.sp$GEOID
county2 <- geometry(spChFIDs(TX.sp, county)) 
counties <- spTransform(county2, CRS.new)
county<-as.numeric(county)

#Load pop
Pop <- read.csv("Population_final2.csv", header=T, sep=";", dec=".")
Pop <- Pop[,-2]
Pop<- Pop[,-48]
Pop.df = melt(Pop, id.vars = "FIP")
Pop.df$Year = as.numeric(substring(Pop.df$variable, first = 4, last = 7))
names(Pop.df)[3:4] = c("pop", "YearPop")
Pop.df$lpop = log10(Pop.df$pop)
Pop.df$ID <- "" 
Pop.df$ID<-match(Pop.df$FIP, county)



###POP Change
##pop changes by county #http://pages.uoregon.edu/rgp/PPPM613/class8a.htm
PC <- Pop.df %>% group_by(ID) %>%
  summarize(Change = (pop[YearPop == max(YearPop)] - pop[YearPop == min(YearPop)])/pop[YearPop == min(YearPop)] * 100)
PC.df = as.data.frame(PC)
row.names(PC.df) = county


####Preparatrion for inla
#1. Subset Texas from big shape
TornTexas<-subset(TornALL, TornALL$st=="TX")
TornTexas<-subset(TornTexas, TornTexas$Ref != 2949 &
                    TornTexas$Ref != 5732 &
                    TornTexas$Ref != 9752 &
                    TornTexas$Ref != 10216 &
                    TornTexas$Ref != 10641 &
                    TornTexas$Ref != 57954 & yr>=1970)
#2. Return Number of tornados, first by state, per year, starting in 1970
ct = over(counties, TornTexas, returnList = TRUE)
names(ct) = county
TornAll <- ldply(ct, data.frame)
nTornados <- TornAll %>% filter(!duplicated(Ref)) %>%
  group_by(yr, .id) %>%
  dplyr::summarize(numberTorn = n())

#3. Creation of Dataframe Number counts/county/year
years<-c(1970:2015)
yearTorn<-rep(years, each=length(county))

countyTorn<-rep(sort(county), times=length(years))

Random.df<- data.frame(County=countyTorn, Year=yearTorn)
colnames(nTornados) [1:2] <- c("Year", "County")
TornInla <-  merge(Random.df,nTornados,by=c("Year","County"), all=TRUE)
TornInla[is.na(TornInla)] <- 0

#4. Prep. of INLA graph
spdf = SpatialPolygonsDataFrame(counties, PC.df)
spdf$ID<-seq(1:254)
View(spdf)
spdf$area = round((rgeos::gArea(counties, byid = TRUE)/10^6), 5)
spdf$Name = TX.sp$NAME
spdf$FIP = county
a <- aggregate(TornInla$numberTorn, by=list(TornInla$County), FUN=sum)
colnames(a)<-c("FIP", "nT")
spdf<-merge(spdf, a, by="FIP")
View(spdf)

nb = poly2nb(spdf)
nb2INLA("tornb.inla", nb)
tornb.inla = inla.read.graph("tornb.inla")
image(inla.graph2matrix(tornb.inla),xlab="",ylab="")

##Generate spatial ID in INLAtORN. DF
TornInla$ID <- "" 
spdf234<-as.data.frame(spdf)
spdf234<-arrange(spdf234, FIP)
ID<-rep(spdf234$ID, 46)
TornInla$ID=ID

##Now ID in INLAtable and ID in INLA graph match!!!
TornInla$ID2<-TornInla$ID
TornInla$Year2<-TornInla$Year

##several controls:
control <- list(
  predictor = list(compute = TRUE),
  results = list(return.marginals.random = TRUE, return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
                 po = TRUE, waic=TRUE, graph=TRUE, gdensity=TRUE, openmp.strategy="huge"), 
  group = list(model="rw2"))

#1 Model linear trend (or non-linear) only of tornados/year
b<- aggregate(TornInla$numberTorn, by=list(TornInla$Year), FUN=sum)
colnames(b)<-c("Year", "nT")
b<-b[-47,]
spdf<-merge(spdf, a, by="FIP")

formula<-nT ~ Year
modelt0 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt0)

b$Year2<-b$Year

#non-linear trend
formula<-nT ~ Year + f(Year2, model="iid")
modelt1 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt1)

###non -linear trend better fit!! But which model? - by DIC, CRW2
formula<-nT ~ Year + f(Year2, model="rw2")
modelt2 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt2)

formula<-nT ~ Year + f(Year2, model="rw1")
modelt3 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt3)

formula<-nT ~ Year + f(Year2, model="crw2")
modelt4 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt4)

formula<-nT ~ Year + f(Year2, model="mec")
modelt5 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt5)

formula<-nT ~ Year + f(Year2, model="meb")
modelt6 = inla(formula = formula, 
               family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = b,
               control.compute = control$compute,
               control.predictor = control$predictor)
summary(modelt6)

#Prepare data for spatial component
nTornadospatial <- TornInla %>% 
  group_by(County) %>%
  dplyr::summarize(numberTorn = sum(numberTorn))
nTornadospatial$ID<-spdf234$ID

############################Frailty model
#simple random effect (spatially unstructured heterogeneity model)

frailtyformula<-numberTorn~f(ID, model="iid")
E <- mean(nTornadospatial$numberTorn)
frailtymodel=inla(formula=frailtyformula, family = "poisson", 
                  data=nTornadospatial, E=E,
                  quantiles = c(.05, .5, .95),
                  control.compute = control$compute, 
                  control.results = control$results, 
                  control.predictor = control$predictor)
summary(frailtymodel)

brier.score <- function(x, m){
  with(m, {mean(x^2) - 2 * mean(x * mean) + mean(mean^2 + sd^2)})
}
-mean(log(frailtymodel$cpo$cpo))
brier.score(nTornadospatial[["numberTorn"]], frailtymodel[["summary.fitted.values"]])
goftest::cvm.test(frailtymodel$cpo$pit, null="punif")

exp(frailtymodel$summary.fixed)

####analyse random effects
refm<-exp(frailtymodel$summary.random$ID[,2])
a<-data.frame(refm)
a <- ggplot(a, aes(refm))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Density Plot - Spatial Random Effect Distribution \n Frailty Model")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")
a
a<-data.frame(refm)
#plot random effects
#1.merge dataframes
spdf_img<-spdf[,3]
spdf_img$re1<-a[,1]
#3. spplot
range(spdf_img$re1)
rng = c(seq(0, 4, length=5), 9)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-900000,-1100000), 
             scale = 300000, fill=c("transparent","black"))
text1 = list("sp.text", c(-900000,-1150000), "0")
text2 = list("sp.text", c(-550000,-1150000), "300 Km")
text4<-list("sp.text", c( -730000, -1270000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-900000, -400000), scale = 200000)
spplot(spdf_img, "re1", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Random Effects", 
       main="Spatial Unstructured Heterogeneity \n Occurence of Tornados in Texas")


#####Convolution model
convolformula <- numberTorn ~ f(ID, model = "bym", graph=tornb.inla)
convolmodel <- inla(formula = convolformula, family = "poisson", 
                    quantiles = c(.05, .5, .95),
                    data = nTornadospatial, E=E,
                    control.compute = control$compute,
                    control.predictor = control$predictor)


summary(convolmodel)

exp(convolmodel$summary.fixed)
-mean(log(convolmodel$cpo$cpo))
brier.score(nTornadospatial[["numberTorn"]], convolmodel[["summary.fitted.values"]])
goftest::cvm.test(convolmodel$cpo$pit, null="punif")


####analyse random effects
refm1<-exp(convolmodel$summary.random$ID[1:254,2])
plot(density(refm))
a<-data.frame(refm)
a <- ggplot(a, aes(refm))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Density Plot - Spatially Unstructured Effects Distribution \n Convolution Model")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")
a
refm<-exp(convolmodel$summary.random$ID[255:508,2])
plot(density(refm))
a<-data.frame(refm)
a <- ggplot(a, aes(refm))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Density Plot - Spatially Structured Effects Distribution \n Convolution Model")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")
a

#plot random effects
#1.merge dataframes
re2<-refm1
spdf_img$re2<-re2
#3. spplot
range(spdf_img$re2)
rng = c(seq(0, 4, length=5), 9)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-900000,-1100000), 
             scale = 300000, fill=c("transparent","black"))
text1 = list("sp.text", c(-900000,-1150000), "0")
text2 = list("sp.text", c(-550000,-1150000), "300 Km")
text4<-list("sp.text", c( -730000, -1270000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-900000, -400000), scale = 200000)
spplot(spdf_img, "re2", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Random Effects", 
       main="Spatial Unstructured Heterogeneity \n Occurence of Tornados in Texas \n Convolution Model")

re3<-refm
spdf_img$re3<-re3
#3. spplot
range(spdf_img$re3)
rng = c(seq(0, 4, length=5), 7)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-900000,-1100000), 
             scale = 300000, fill=c("transparent","black"))
text1 = list("sp.text", c(-900000,-1150000), "0")
text2 = list("sp.text", c(-550000,-1150000), "300 Km")
text4<-list("sp.text", c( -730000, -1270000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-900000, -400000), scale = 200000)
spplot(spdf_img, "re3", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Spatial Effects", 
       main="Spatial Structured Heterogeneity \n Occurence of Tornados in Texas \n Convolution Model")


#fitted
CARfit<-NULL
CARfit<-convolmodel$summary.fitted.values[,1]
CARfit<-data.frame(CARfit=CARfit, ID=nTornadospatial$ID)
View(spdf_img)
spdf_img<-spdf[,3]

spdf_img<-merge(spdf_img, CARfit, by="ID")
range(spdf_img$CARfit)
rng = c(seq(0, 4, length=5), 8)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-900000,-1100000), 
             scale = 300000, fill=c("transparent","black"))
text1 = list("sp.text", c(-900000,-1150000), "0")
text2 = list("sp.text", c(-550000,-1150000), "300 Km")
text4<-list("sp.text", c( -730000, -1270000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-900000, -400000), scale = 200000)
spplot(spdf_img, "CARfit", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Fitted Effects", 
       main="Fitted Effects \n Occurence of Tornados in Texas \n Convolution Model")



#spatial risk
CARmarginals<-convolmodel$marginals.random$ID[1:254]
CARzeta<-lapply(CARmarginals, function (x) inla.emarginal(exp, x))
risk<-data.frame(CARzeta=unlist(CARzeta), ID=seq(1, 254, 1))
risk<-merge(nTornadospatial, risk, by="ID")
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, risk, by="ID")
View(spdf_img)
rng = c(seq(0, 4, length=5), 9)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
spplot(spdf_img, "CARzeta", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Marginal Effects", 
       main="Spatial Risk \n Occurence of Tornados in Texas \n Convolution Model")


#######spatio-temporal
#SP+TIME uncorrelated
E=mean(TornInla$numberTorn)
convolformulat <- numberTorn ~ f(ID, model = "bym", graph=tornb.inla)+ f(Year, model="iid")
convolmodelt <- inla(formula = convolformulat, family = "poisson", 
                     quantiles = c(.05, .5, .95),
                     data = TornInla, E=E,
                     control.compute = control$compute,
                     control.predictor = control$predictor)
summary(convolmodelt)
exp(convolmodelt$summary.fixed)
-mean(log(convolmodelt$cpo$cpo))
brier.score(TornInla[["numberTorn"]], convolmodelt[["summary.fitted.values"]])
goftest::cvm.test(convolmodelt$cpo$pit, null="punif")

re<-convolmodelt$summary.random$ID[1:254,2]
plot(density(re))
a<-data.frame(re)
a <- ggplot(a, aes(re))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Spatially Unstructured Effects")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")

care<-convolmodelt$summary.random$ID[255:508, 2]
a<-data.frame(care)
a <- ggplot(a, aes(care))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Spatially Structured Effects")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")

tre<-convolmodelt$summary.random$Year[,2]
a<-data.frame(tre)
a <- ggplot(a, aes(tre))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Temporal Unstructured Effects")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")


##sp-time correlated
convolformulat2<-numberTorn ~f (ID, model="bym", graph = tornb.inla)+
  f(Year, model="rw1")
convolmodelt2<-inla(formula = convolformulat2, family = "poisson", 
                    quantiles = c(.05, .5, .95),
                    data = TornInla, E=E,
                    control.compute = control$compute,
                    control.predictor = control$predictor)

summary(convolmodelt2)
exp(convolmodelt2$summary.fixed)
-mean(log(convolmodelt2$cpo$cpo))
brier.score(TornInla[["numberTorn"]], convolmodelt2[["summary.fitted.values"]])
goftest::cvm.test(convolmodelt2$cpo$pit, null="punif")


re<-convolmodelt2$summary.random$ID[1:254,2]
a<-data.frame(re)
a <- ggplot(a, aes(re))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Spatially Unstructured Effects")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")

care<-convolmodelt2$summary.random$ID[255:508, 2]
a<-data.frame(care)
a <- ggplot(a, aes(care))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Spatially Structured Effects")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")

tre<-convolmodelt2$summary.random$Year[,2]
a<-data.frame(tre)
a <- ggplot(a, aes(tre))
a+geom_density(fill="cadetblue4", alpha=0.2, colour="azure4")+
  ggtitle("Temporal Structured Effects")+
  theme(plot.title = element_text(hjust = 0.5))+
  xlab("N=254")+ylab("Density")


#map fitted effects
carfit<-data.frame(FIT=convolmodelt2$summary.fitted.values[,1], ID=TornInla$ID, Year=TornInla$Year)
meanfit<-carfit %>% 
  group_by(ID) %>%
  dplyr::summarize(meannumberTorn = mean(FIT))
spdf_img<-merge(spdf_img, meanfit, by="ID")

rng = c(seq(0, 4, length=5), 7)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
spplot(spdf_img, "meannumberTorn", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Mean Fitted Effects", 
       main="Mean Spatio-temporal Fitted Effects")


#Marginal Effects
CARmarginals<-convolmodelt2$marginals.random$ID[1:254]
CARzeta<-lapply(CARmarginals, function (x) inla.emarginal(exp, x))
risk<-data.frame(CARzeta=unlist(CARzeta), ID=seq(1, 254, 1))
risk<-merge(nTornadospatial, risk, by="ID")
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, risk, by="ID")

rng = c(seq(0, 4, length=5), 9)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
spplot(spdf_img, "CARzeta", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Marginal Effects", 
       main="Spatial Risk \n Occurence of Tornados in Texas")


df = convolmodelt2$summary.fitted.values
names(df) = c("mean", "sd", "QL", "QM", "QH", "mode")
df$ID<-TornInla$ID
df$Year<-TornInla$Year
df = df %>% group_by(ID) %>% 
  dplyr::summarize(mean=mean(mean), sd=mean(sd), 
                   QL=mean(QL), QM=mean(QM), 
                   QH= mean(QH))

df<- df%>% mutate(QL = QL - 1,
                  QH = QH - 1,
                  Sig = sign(QL) == sign(QH),
                  sd = sd,
                  ctyPerState = (mean - 1)*100) 

sum(df$Sig)

spdfR = spdf
spdfR@data = df

range(spdfR$ctyPerState)
rng = c(seq(-100, 100, length=8), 200, 600)
rnq = c(rev(RColorBrewer::brewer.pal(8, "RdYlGn")), "#8c510a", "#543005")
spplot(spdfR, "ctyPerState", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Percentage of difference from Statewide Average Rate", 
       main="Occurence rate of Tornados in Texas")

range(spdfR$sd)
rng<-seq(0, 0.80, length=9)
rnq<-c("#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081")
spplot(spdfR, "sd", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 2))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Standard Error")



###################################ADD COVARIATES 
TornInla$Pop<-Pop.df$pop
TornInla$Lpop<-Pop.df$lpop
a<-data.frame(area=round(spdf$area, 5), County=spdf$FIP)

TornInla<-merge(TornInla, a, by="County")
TornInla <- dplyr::arrange(TornInla, Year)
TornInla$DPop<-""
TornInla$DPop<-TornInla$Pop/TornInla$area

convolformulatc2<-numberTorn ~f (ID, model="bym", graph = tornb.inla)+
  f(Year, model="rw1") + DPop
convolmodeltc2<-inla(formula = convolformulatc2, family = "poisson", 
                    quantiles = c(.05, .5, .95),
                    data = TornInla, E=E,
                    control.compute = control$compute,
                    control.predictor = control$predictor)

summary(convolmodeltc2)
exp(convolmodeltc2$summary.fixed)
-mean(log(convolmodeltc2$cpo$cpo))
brier.score(TornInla[["numberTorn"]], convolmodeltc2[["summary.fitted.values"]])
goftest::cvm.test(convolmodeltc2$cpo$pit, null="punif")

###Insert Roughness Index
wcounty<-unionSpatialPolygons(counties, ID = rep("1", length(row.names(counties))))
Ind<-raster("Index_Value1.tif")
Ind<-as(crop(Ind, extent(wcounty)), "SpatialGridDataFrame")
proj4string(Ind) = proj4string(wcounty) #same projection, different datum & ellipsoid
RI.data<-over(counties, Ind, returnList = TRUE)

Elev.df = data.frame(county = rep(county, sapply(RI.data, nrow)),
                     Elev = unlist(RI.data), 
                     ID = rep(spdf$ID, sapply(RI.data, nrow)),
                     stringsAsFactors = FALSE)

CE.df = Elev.df %>% group_by(ID) %>%  
  dplyr::summarize(elev = mean(Elev, na.rm = TRUE),
                   elevS = sd(Elev, na.rm = TRUE),
                   elevCV = elevS/elev)
all(spdf$ID == CE.df$ID)
TornInla = merge(TornInla, CE.df, by = "ID")


####Models with st of RI
formulamodeltc3 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + f(Year, model="rw1") +elevS
modeltc3 <- inla(formula = formulamodeltc3, family = "poisson", 
               quantiles = c(.05, .5, .95), 
               data = TornInla, E=E,
               control.compute = control$compute,
               control.predictor = control$predictor, 
               control.results=control$results)

summary(modeltc3)
exp(modeltc3$summary.fixed)
-mean(log(modeltc3$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc3[["summary.fitted.values"]])
goftest::cvm.test(modeltc3$cpo$pit, null="punif")


####INSERT LAND-COVER
###Adding Land-Cover
lc1992<-read.dbf("./Area_Per_County_sqM/1992.dbf")
lc2001<-read.dbf("./Area_Per_County_sqM/2001.dbf")
lc2006<-read.dbf("./Area_Per_County_sqM/2006.dbf")
lc2011<-read.dbf("./Area_Per_County_sqM/2011.dbf")

# percentages:
tarea<-data.frame(Area= spdf$area*10^6, GEOID=spdf$FIP)
lc1992<-merge(lc1992, tarea, by="GEOID")
lc1992$perc11<-lc1992$VALUE_11/lc1992$Area*100
lc1992$perc21<-lc1992$VALUE_21/lc1992$Area*100
lc1992$perc31<-lc1992$VALUE_31/lc1992$Area*100
lc1992$perc41<-lc1992$VALUE_41/lc1992$Area*100
lc1992$perc51<-lc1992$VALUE_51/lc1992$Area*100
lc1992$perc91<-lc1992$VALUE_91/lc1992$Area*100


lc2001<-merge(lc2001, tarea, by="GEOID")
lc2001$perc11<-lc2001$VALUE_11/lc2001$Area*100
lc2001$perc21<-lc2001$VALUE_21/lc2001$Area*100
lc2001$perc31<-lc2001$VALUE_31/lc2001$Area*100
lc2001$perc41<-lc2001$VALUE_41/lc2001$Area*100
lc2001$perc51<-lc2001$VALUE_51/lc2001$Area*100
lc2001$perc91<-lc2001$VALUE_91/lc2001$Area*100

lc2006<-merge(lc2006, tarea, by="GEOID")
lc2006$perc11<-lc2006$VALUE_11/lc2006$Area*100
lc2006$perc21<-lc2006$VALUE_21/lc2006$Area*100
lc2006$perc31<-lc2006$VALUE_31/lc2006$Area*100
lc2006$perc41<-lc2006$VALUE_41/lc2006$Area*100
lc2006$perc51<-lc2006$VALUE_51/lc2006$Area*100
lc2006$perc91<-lc2006$VALUE_91/lc2006$Area*100

lc2011<-merge(lc2011, tarea, by="GEOID")
lc2011$perc11<-lc2011$VALUE_11/lc2011$Area*100
lc2011$perc21<-lc2011$VALUE_21/lc2011$Area*100
lc2011$perc31<-lc2011$VALUE_31/lc2011$Area*100
lc2011$perc41<-lc2011$VALUE_41/lc2011$Area*100
lc2011$perc51<-lc2011$VALUE_51/lc2011$Area*100
lc2011$perc91<-lc2011$VALUE_91/lc2011$Area*100

a<-rep(lc1992$perc11, times= 23)
b<-rep(lc2001$perc11, times=10)
c<-rep(lc2006$perc11, times=5)
d<-rep(lc2011$perc11, times=8)

asd<-c(a,b,c,d)

TornInla$perc11<-asd

a<-rep(lc1992$perc21, times= 23)
b<-rep(lc2001$perc21, times=10)
c<-rep(lc2006$perc21, times=5)
d<-rep(lc2011$perc21, times=8)

asd<-c(a,b,c,d)

TornInla$perc21<-asd

a<-rep(lc1992$perc41, times=23)
b<-rep(lc2001$perc41, times=10)
c<-rep(lc2006$perc41, times=5)
d<-rep(lc2011$perc41, times=8)

asd<-c(a,b,c,d)
TornInla$perc41<-asd

a<-rep(lc1992$perc51, times=23)
b<-rep(lc2001$perc51, times=10)
c<-rep(lc2006$perc51, times=5)
d<-rep(lc2011$perc51, times=8)

asd<-c(a,b,c,d)

TornInla$perc51<-asd

a<-rep(lc1992$perc31, times= 23)
b<-rep(lc2001$perc31, times=10)
c<-rep(lc2006$perc31, times=5)
d<-rep(lc2011$perc31, times=8)

asd<-c(a,b,c,d)

TornInla$perc31<-asd

a<-rep(lc1992$perc91, times= 23)
b<-rep(lc2001$perc91, times=10)
c<-rep(lc2006$perc91, times=5)
d<-rep(lc2011$perc91, times=8)

asd<-c(a,b,c,d)

TornInla$perc91<-asd


formulamodeltc4 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + f(Year, model="rw1") + perc11
modeltc4 <- inla(formula = formulamodeltc4, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)

summary(modeltc4)
exp(modeltc4$summary.fixed)
-mean(log(modeltc4$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc4[["summary.fitted.values"]])
goftest::cvm.test(modeltc4$cpo$pit, null="punif")


formulamodeltc5 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + f(Year, model="rw1") + perc21
modeltc5 <- inla(formula = formulamodeltc5, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)

summary(modeltc5)
exp(modeltc5$summary.fixed)
-mean(log(modeltc5$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc5[["summary.fitted.values"]])
goftest::cvm.test(modeltc5$cpo$pit, null="punif")


formulamodeltc6 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="rw1") + perc31
modeltc6 <- inla(formula = formulamodeltc6, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)

summary(modeltc6)
exp(modeltc6$summary.fixed)
-mean(log(modeltc6$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc6[["summary.fitted.values"]])
goftest::cvm.test(modeltc6$cpo$pit, null="punif")


formulamodeltc7 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="rw1") + perc41
modeltc7 <- inla(formula = formulamodeltc7, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)

summary(modeltc7)
exp(modeltc7$summary.fixed)
-mean(log(modeltc7$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc7[["summary.fitted.values"]])
goftest::cvm.test(modeltc7$cpo$pit, null="punif")


formulamodeltc8 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="rw1") + perc51
modeltc8 <- inla(formula = formulamodeltc8, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)
summary(modeltc8)
exp(modeltc8$summary.fixed)
-mean(log(modeltc8$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc8[["summary.fitted.values"]])
goftest::cvm.test(modeltc8$cpo$pit, null="punif")


formulamodeltc9 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="rw1") + perc91
modeltc9 <- inla(formula = formulamodeltc9, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)

summary(modeltc9)
exp(modeltc9$summary.fixed)
-mean(log(modeltc9$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc9[["summary.fitted.values"]])
goftest::cvm.test(modeltc9$cpo$pit, null="punif")


##ALL covariates
formulamodeltc5 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="rw1") + perc21+perc31+perc41+perc11+perc51+perc91+elevS+DPop
modeltc5 <- inla(formula = formulamodeltc5, family = "poisson", 
                 quantiles = c(.05, .5, .95), 
                 data = TornInla, E=E,
                 control.compute = control$compute,
                 control.predictor = control$predictor, 
                 control.results=control$results)

MOD=modeltc5
par(mfrow=c(2,4))

plot(MOD$marginals.fixed$elevS, main="SDTPI", xlab="", ylab="", xlim=c(-5,5),  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc21, main="(%) Residential", xlim=c(-0.1,0.1),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc51, main="(%) Low-Grass", xlim=c(-0.1,0.1), xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc31, main=" (%) Barren", xlim=c(-0.1,0.1),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc41, main=" (%) Forest",xlim=c(-0.1,0.1),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc11, main= "(%) Water", xlim=c(-0.1,0.1),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$DPop, main="Pop Density", xlim=c(-0.1,0.1),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$elevS, main="(%) Wetlands", xlim=c(-5,5),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
dev.off()

exp(MOD$summary.fixed)



#####Bernardinelli
formulamodeltc11 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="iid") + Year2+
  elevS + perc31 + perc41 + perc11 + perc51 + perc91 + perc21

modeltc11 <- inla(formula = formulamodeltc11, family = "poisson", 
                  quantiles = c(.05, .5, .95), 
                  data = TornInla, E=E,
                  control.compute = control$compute, 
                  control.results=control$results)

summary(modeltc11)
exp(modeltc11$summary.fixed)
-mean(log(modeltc11$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc11[["summary.fitted.values"]])
goftest::cvm.test(modeltc11$cpo$pit, null="punif")


##Knorr held
formulamodeltc12 <- numberTorn ~f(ID, model = "bym", graph=tornb.inla) + 
  f(Year, model="rw1")+f(Year2, model = "iid") +
  elevS + perc11+perc21+perc31 + perc41 +perc51+perc91
modeltc12 <- inla(formula = formulamodeltc12, family = "poisson", 
                  quantiles = c(.05, .5, .95), 
                  data = TornInla, E=E,
                  control.compute = control$compute, 
                  control.results=control$results)

summary(modeltc12)
exp(modeltc12$summary.fixed)
-mean(log(modeltc11$cpo$cpo))
brier.score(TornInla[["numberTorn"]], modeltc11[["summary.fitted.values"]])
goftest::cvm.test(modeltc11$cpo$pit, null="punif")


#type I interaction
TornInla$area.year <- seq(1,length(countyTorn))

formTypeI <- numberTorn~ + f(ID, model="bym", graph=tornb.inla)+
  f(Year, model="rw1") + f(Year2, model="iid")+
  f(area.year, model="iid")+
  perc11+perc21+perc31+perc41+perc51+perc91+elevS

mod.intI <- inla(formTypeI,family="poisson",data=TornInla,
                 control.predictor=control$predictor,
                 control.compute=control$compute)

summary(mod.intI)
exp(mod.intI$summary.fixed)
-mean(log(mod.intI$cpo$cpo))
brier.score(TornInla[["numberTorn"]], mod.intI[["summary.fitted.values"]])
goftest::cvm.test(mod.intI$cpo$pit, null="punif")

save.image("hopefully.RData")

#exploring final model:
MOD<-mod.intI
exp(MOD$summary.fixed)

summary(MOD)
#map fitted effects
carfit<-data.frame(FIT=MOD$summary.fitted.values[,1], ID=TornInla$ID, Year=TornInla$Year)
meanfit<-carfit %>% 
  group_by(ID) %>%
  dplyr::summarize(meannumberTorn = mean(FIT))
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, meanfit, by="ID")

range(spdf_img$meannumberTorn)
rng = c(seq(0, 2, length=5), 5)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
spplot(spdf_img, "meannumberTorn", col = "white", at = rng,
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))), 
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Mean Fitted Effects", 
       main="Mean Spatio-temporal Fitted Effects")


#spatial risk
#Marginal Effects
CARmarginals<-MOD$marginals.random$ID[1:254]
CARzeta<-lapply(CARmarginals, function (x) inla.emarginal(exp, x))
risk<-data.frame(CARzeta=unlist(CARzeta), ID=seq(1, 254, 1))
risk<-merge(nTornadospatial, risk, by="ID")
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, risk, by="ID")

rng = c(seq(0, 4, length=5), 9)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
spplot(spdf_img, "CARzeta", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Marginal Effects", 
       main="Spatial Risk \n Occurence of Tornados in Texas")


#exceedence Probability
a=log(1)
stexceed<-lapply(MOD$marginals.random$ID[1:254], 
                 function (X) {
                   1-inla.pmarginal(a, X) 
                 })
stexceed<-unlist(stexceed)
risk<-data.frame(stexceed=stexceed, ID=seq(1, 254, 1))
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, risk, by="ID")
rng = seq(0, 1, length=10)
rnq= rev(RColorBrewer::brewer.pal(9, "RdYlGn"))
spplot(spdf_img, "stexceed", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       main = "More than one tornado per county", 
       sub="Probability")


a=log(2)
stexceed<-lapply(MOD$marginals.random$ID[1:254], 
                 function (X) {
                   1-inla.pmarginal(a, X) 
                 })
stexceed<-unlist(stexceed)
risk<-data.frame(stexceed=stexceed, ID=seq(1, 254, 1))
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, risk, by="ID")
rng = seq(0, 1, length=10)
rnq= rev(RColorBrewer::brewer.pal(9, "RdYlGn"))
spplot(spdf_img, "stexceed", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       main = "More than two tornados per county", 
       sub="Probability")


pop density map

dens<-subset(TornInla, TornInla$Year==2015 )
dens1<-merge(spdf, dens, by="ID")
rng = c(seq(0, 100, length=8),200, 600, 1000, 1200)
rnq= rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
spplot(dens1, "DPop", col = "grey", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=c(0,100,200,600,1000,1200))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       main = "Population Density - Texas (2015)", 
       sub="Individuals per Sq Km")
