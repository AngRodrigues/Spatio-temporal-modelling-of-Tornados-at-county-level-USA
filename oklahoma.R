setwd("D:/OK")

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
TornALL <- readOGR(dsn = "./tmp/torn", layer = "torn", stringsAsFactors = FALSE)
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

TornOK<-subset(TornALL, TornALL$st=="OK")
Torn_ok<-as.data.frame(TornOK)
TornT = as.data.frame(table(Torn_ok$yr, Torn_ok[,11]))
TornT$year = as.numeric(levels(TornT$Var1))
TornT$Fscale = paste("F", TornT$Var2, sep = "")
ggplot(TornT[TornT$Var2 != -9, ], aes(x = year, y = Freq)) + 
  geom_point()+ geom_smooth(span = 0.9, color="forestgreen") + 
  facet_wrap(~Fscale, ncol = 2, scales = "free") + 
  theme_gray()+ 
  ggtitle("Number of Reported Tornados per Year \n (Oklahoma) \n 1970-2015 \n by FScale")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  ylab("Reported Number of Tornadoes")

#Load Boundaries
US.sp <- readOGR(dsn = "./tmp", layer = "cb_2013_us_county_5m", 
                 stringsAsFactors = FALSE)
OK.sp <- US.sp[US.sp$STATEFP == 40, ]
county <- OK.sp$GEOID
county2 <- geometry(spChFIDs(OK.sp, county)) 
counties <- spTransform(county2, CRS.new)
county<-as.numeric(county)

###Insert Pop
Pop <- read.csv("Population.csv", header=T, sep=";", dec=".")
Pop <- Pop[,-2]
Pop$pop2015<- Pop$pop2014
colnames(Pop)[1]<-"FIP"
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
Torn_OK<-subset(TornALL, TornALL$st=="OK" & yr>=1970)
#2. Return Number of tornados, first by state, per year, starting in 1970
ct = over(counties, Torn_OK, returnList = TRUE)
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
spdf$ID<-seq(1:77)

spdf$area = round((rgeos::gArea(counties, byid = TRUE)/10^6), 5) #FOR KM2
spdf$Name = OK.sp$NAME
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


#number of tornados counts over the years
range(spdf$nT)
rng = seq(10, 80, 10)
crq = wesanderson::wes_palette("Moonrise1", 7, "continuous")
spdf@bbox
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-600000,-300000), 
             scale = 150000, fill=c("transparent","black"))
text1 = list("sp.text", c(-600000,-310000), "0")
text2 = list("sp.text", c(-430000,-310000), "300 Km")
text3<-list("sp.text", c(-500000, -370000), "Source: SPC, 2016a", cex=0.6)
text4<-list("sp.text", c( -500000, -400000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-550000, -220000), scale = 100000)
spplot(spdf, "nT", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(
         space = "bottom", labels=list(at=rng)),
       sp.layout=list(scale, text1, text2, text3, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Number of tornados", main="Number of tornados per county in Oklahoma \n 1970-2015")

pc.sp<-setDT(PC.df, keep.rownames = TRUE)[]
colnames(pc.sp)[1]<-"FIP"
spdf<-merge(spdf, pc.sp, by="FIP")
range(spdf$Change.x)

rng = round(seq(-50, 340, length=10), 0)
crq = wesanderson::wes_palette("Zissou", 10, "continuous")
spdf@bbox
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-600000,-300000), 
             scale = 150000, fill=c("transparent","black"))
text1 = list("sp.text", c(-600000,-310000), "0")
text2 = list("sp.text", c(-430000,-310000), "300 Km")
text3<-list("sp.text", c(-500000, -370000), "Source: SPC, 2016a", cex=0.6)
text4<-list("sp.text", c( -500000, -400000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-550000, -220000), scale = 100000)
spplot(spdf, "Change.x", col = "white", at = rng, 
       col.regions = crq,
       colorkey = list(
         space = "bottom", labels=list(at=rng)),
       sp.layout=list(scale, text1, text2, text3, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Percentage of Change", main="Population change per county in Oklahoma \n 1970-2015")




#Prepare data for spatial component
nTornadospatial <- TornInla %>% 
  group_by(County) %>%
  dplyr::summarize(numberTorn = sum(numberTorn))
nTornadospatial$ID<-spdf234$ID

##several controls:
control <- list(
  predictor = list(compute = TRUE),
  results = list(return.marginals.random = TRUE, return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
                 po = TRUE, waic=TRUE, graph=TRUE, gdensity=TRUE, openmp.strategy="huge"), 
  group = list(model="rw2"))

E <- mean(nTornadospatial$numberTorn)/46

#1. Convol model:
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


###################################ADD COVARIATES 
TornInla$Pop<-(Pop.df$pop/100)
TornInla$Lpop<-Pop.df$lpop
a<-data.frame(area=round(spdf$area, 5), County=spdf$FIP)

TornInla<-merge(TornInla, a, by="County")
TornInla <- dplyr::arrange(TornInla, Year)
TornInla$DPop<-""
TornInla$DPop<-TornInla$Pop/TornInla$area

###Insert Roughness Index
wcounty<-unionSpatialPolygons(counties, ID = rep("1", length(row.names(counties))))
Ind<-raster("./data/TPI1.tif")
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


#################LANDCOVER
####INSERT LAND-COVER
###Adding Land-Cover
lc1992<-read.dbf("./data/1992.dbf")
lc2001<-read.dbf("./data/2001.dbf")
lc2006<-read.dbf("./data/2006.dbf")
lc2011<-read.dbf("./data/20112.dbf")


# percentages:
tarea<-data.frame(Area= spdf$area*10^9, GEOID=spdf$FIP)
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

#no covariates covariates
convolformulat3<-numberTorn ~f (ID, model="bym", graph = tornb.inla)+
  f(Year, model="rw1")
convolmodelt3<-inla(formula = convolformulat3, family = "poisson", 
                    quantiles = c(.05, .5, .95),
                    data = TornInla, E=E,
                    control.compute = control$compute,
                    control.predictor = control$predictor)
summary(convolmodelt3)
exp(convolmodelt3$summary.fixed)
-mean(log(convolmodelt3$cpo$cpo))
brier.score(TornInla[["numberTorn"]], convolmodelt3[["summary.fitted.values"]])
goftest::cvm.test(convolmodelt3$cpo$pit, null="punif")



####sp time + covariates
##sp-time correlated
convolformulat1<-numberTorn ~f (ID, model="bym", graph = tornb.inla)+
  f(Year, model="rw1")+
  perc11+perc21+perc31+perc41+perc51+perc91+elevS
convolmodelt1<-inla(formula = convolformulat1, family = "poisson", 
                    quantiles = c(.05, .5, .95),
                    data = TornInla, E=E,
                    control.compute = control$compute,
                    control.predictor = control$predictor)

summary(convolmodelt1)
exp(convolmodelt1$summary.fixed)
-mean(log(convolmodelt1$cpo$cpo))
brier.score(TornInla[["numberTorn"]], convolmodelt1[["summary.fitted.values"]])
goftest::cvm.test(convolmodelt1$cpo$pit, null="punif")

MOD=convolmodelt1
par(mfrow=c(2,4))

plot(MOD$marginals.fixed$elevS, main="SDTPI", xlab="", ylab="", xlim=c(-5,5),  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc21, main="(%) Residential", xlim=c(-50,50),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc51, main="(%) Low-Grass", xlim=c(-50,50), xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc31, main=" (%) Barren", xlim=c(-50,50),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc41, main=" (%) Forest",xlim=c(-50,50),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$perc11, main= "(%) Water", xlim=c(-50,50),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
plot(MOD$marginals.fixed$elevS, main="(%) Wetlands", xlim=c(-5,5),  xlab="", ylab="",  yaxt='n', type="l")
abline(v=0, col="red")
dev.off()


exp(MOD$summary.fixed)


##SPTIME INTERACTIONS

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


#Knorr
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


#sp-time interaction
TornInla$area.year <- seq(1,length(countyTorn))

formTypeI <- numberTorn~ + f(ID, model="bym", graph=tornb.inla)+
  f(Year, model="crw2") + f(Year2, model="iid")+
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


formTypeIa <- numberTorn~ + f(ID, model="bym", graph=tornb.inla)+
  f(Year, model="rw1") + f(Year2, model="iid")+
  f(area.year, model="iid")

mod.intIa <- inla(formTypeIa,family="poisson",data=TornInla,
                 control.predictor=control$predictor,
                 control.compute=control$compute)
summary(mod.intIa)
exp(mod.intIa$summary.fixed)
-mean(log(mod.intIa$cpo$cpo))
brier.score(TornInla[["numberTorn"]], mod.intIa[["summary.fitted.values"]])
goftest::cvm.test(mod.intIa$cpo$pit, null="punif")

formTypeIa <- numberTorn~ + f(ID, model="bym", graph=tornb.inla)+
  f(Year, model="rw1") + f(Year2, model="iid")+
  f(area.year, model="iid")+elevS

mod.intIa <- inla(formTypeIa,family="poisson",data=TornInla,
                  control.predictor=control$predictor,
                  control.compute=control$compute)


MOD<-mod.intI

#spatial risk
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             offset = c(-900000,-1100000), 
             scale = 300000, fill=c("transparent","black"))
text1 = list("sp.text", c(-900000,-1150000), "0")
text2 = list("sp.text", c(-550000,-1150000), "300 Km")
text4<-list("sp.text", c( -730000, -1270000), cex=0.6, "Projection: EPSG 102003")
arrow = list("SpatialPolygonsRescale", layout.north.arrow(), 
             offset = c(-900000, -400000), scale = 200000)

#Marginal Effects
CARmarginals<-MOD$marginals.random$ID[1:77]
CARzeta<-lapply(CARmarginals, function (x) inla.emarginal(exp, x))
risk<-data.frame(CARzeta=unlist(CARzeta), ID=seq(1, 77, 1))
risk<-merge(nTornadospatial, risk, by="ID")
spdf_img<-spdf[,3]
spdf_img<-merge(spdf_img, risk, by="ID")
range(spdf_img$CARzeta)
rng = seq(0, 3, length=4)
rnq = c("#3B9AB2", "#78B7C5", "#EBCC2A", "darkorange1", "#F21A00")
spplot(spdf_img, "CARzeta", col = "white", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       sub = "Marginal Effects", 
       main="Spatial Risk ")


#exceedence Probability
a=log(1)
stexceed<-lapply(MOD$marginals.random$ID[1:77], 
                 function (X) {
                   1-inla.pmarginal(a, X) 
                 })
stexceed<-unlist(stexceed)
risk<-data.frame(stexceed=stexceed, ID=seq(1, 77, 1))
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


#Pop Density Map
dens<-subset(TornInla, TornInla$Year==2015 )
dens1<-merge(spdf, dens, by="ID")
dens1$DPop<-dens1$DPop*100
View(dens1)

rng = c(seq(0, 60, length=5), 200, 500)
rnq= rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
spplot(dens1, "DPop", col = "grey", at = rng, 
       col.regions = rnq,
       colorkey = list(
         space = "bottom", labels=list(
           at=round(rng, 1))),
       sp.layout=list(scale, text1, text2, text4, arrow),
       par.settings = list(axis.line = list(col = NA)),
       main = "Population Density - Oklahoma (2015)", 
       sub="Individuals per Sq Km")
