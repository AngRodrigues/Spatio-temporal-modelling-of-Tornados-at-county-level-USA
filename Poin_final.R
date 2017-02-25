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
  library(nortest)
  library(spacetime)
}
packages()

####DOWNLOAD THE FILE####
download.file("http://www.spc.noaa.gov/gis/svrgis/zipped/tornado.zip",
             "tornado.zip", mode = "wb")
unzip("tornado.zip",exdir="./tmp")
####READ DATA INTO R####
Tornados <- read.dbf("torn/torn.dbf")
Tornados$ref<-""
Tornados$ref<-seq(1,60114, by=1)

####make a simple table to input as point patterns into ARCgis; EPSG: 102003, GET COORDINATES
tornado_events<-as.data.frame(Tornados)
coords<-cbind(tornado_events$slon, tornado_events$slat)
tornados<-SpatialPointsDataFrame(coords, tornado_events, proj4string = CRS("+init=epsg:4326"))
CRS.new <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
tornados <- spTransform(tornados, CRS.new)
tornado_events<-as.data.frame(tornados@coords)
Tornados$x<-""
Tornados$x<-tornado_events$coords.x1
Tornados$y<-tornado_events$coords.x2

#"Total annual counts w/ trend line
##
begin = Tornados$yr[1]
end = as.numeric(Tornados$yr[length(Tornados$yr)])
Count = as.integer(table(Tornados$yr))
AnnualCountALL.df = data.frame(Year = (begin:end), Count = Count)

ggplot(AnnualCountALL.df, aes(x = Year, y = Count)) + geom_line()+
  geom_smooth(method = "gam", color="cadetblue") + ylab("Number of Reported Tornadoes") + 
  theme_gray()+
  ggtitle("Number of Reported Tornados per Year \n (United States of America) \n 1950-2015")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=300, label="Source: SPC 2016a", size=4
  )+
  scale_x_continuous(breaks=c(1950,1960, 1970, 1980, 1990, 2000, 2010))


##plot per Fscale
TornTable = as.data.frame(table(Tornados$yr, Tornados[,11]))
TornTable$year = as.numeric(levels(TornTable$Var1))
TornTable$Fscale = paste("F", TornTable$Var2, sep = "")
ggplot(TornTable[TornTable$Var2 != -9, ], aes(x = year, y = Freq)) + 
  geom_point() + geom_smooth(span = 0.9, color="cadetblue") + 
  facet_wrap(~Fscale, ncol = 2, scales="free_y") + 
  theme_gray()+ 
  ggtitle("Number of Reported Tornados per Year by Fujita Scale\n (United States of America) \n 1950-2015")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  xlab("Year")+ ylab("Reported Number of Tornadoes")+
  scale_x_continuous(breaks=c(1950,1960, 1970, 1980, 1990, 2000, 2010))+
  scale_y_continuous(breaks = scales::pretty_breaks(n=4))


#plot per capital losses
#before 1996, losses are from 0-9; after 1996, they classify it for millions of dollars. 
#give median values for before 1996
Tornados$tloss<-""
Tornados$tloss<-Tornados$loss*1000000
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="0"]<-"0"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="1"]<-"25"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="2"]<-"225"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="3"]<-"2250"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="4"]<-"22500"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="5"]<-"225000"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="6"]<-"2250000"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="7"]<-"22500000"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="8"]<-"225000000"
Tornados$tloss[Tornados$yr %in% c(1950:1995) & Tornados$loss=="9"]<-"5000000000"
Tornados$tloss<-as.numeric(Tornados$tloss)

losses<-aggregate(Tornados$tloss, by=list(Tornados$yr), "sum")
AnnualLossALL.df = data.frame(Year = (begin:end), Count = losses$x/1000000000)

q1<-ggplot(AnnualLossALL.df, aes(x = Year, y = Count)) + 
  geom_point()+geom_smooth(method = "loess", color="gold")+
  ylab("Billion Dollars") + theme_gray()+
  ggtitle("Property Losses")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=-0.5, label="Source: SPC, 2016a", size=2)

AnnualLoss2ALL.df = AnnualLossALL.df[-62,]

q = ggplot(AnnualLoss2ALL.df, aes(x = Year, y = Count)) + geom_point()+geom_smooth(method = "lm")
q + ylab("Billion Dollars") + theme_gray()+
  ggtitle("Total Annual Property Losses in dollars \n (United States of America) \n 1950-2015 \n without 2011")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=-0.5, label="Source: NOAA, SPC, 2016", size=2)


###plot number of fatalities
fatal<-aggregate(Tornados$fat, by=list(Tornados$yr), "sum")
AnnualfatalALL.df = data.frame(Year = (begin:end), Count = fatal$x)

q2<- ggplot(AnnualfatalALL.df, aes(x = Year, y = Count)) + geom_point()+
  geom_smooth(method = "loess", col="skyblue4") + ylab("Number of deaths") + theme_gray()+
  ggtitle("Fatalities")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=-0.5, label="Source: SPC, 2016a", size=2)

###plot number of injuries
inj<-aggregate(Tornados$inj, by=list(Tornados$yr), "sum")
AnnualinjALL.df = data.frame(Year = (begin:end), Count = inj$x)

q3<- ggplot(AnnualfatalALL.df, aes(x = Year, y = Count)) + 
  geom_point()+geom_smooth(method = "loess", col="yellow4")+
  ylab("Number of injuried Individuals") + theme_gray()+
  ggtitle("Injuries")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=-0.5, label="Source: SPC, 2016a", size=2)

multiplot(q1, q2, cols=2)

#plot texas and the other Tornado Alley states
#Estados de Oklahoma, Kansas, Arkansas, Iowa e Missouri 
#Texas, Colorado, Luisiana, Minnesota e Dacota do Sul, 
#Mississippi, Illinois, Indiana, Nebraska, Tennessee e Kentucky Wisconsin.
torn_alley<-subset(Tornados, Tornados$st=="TX" | Tornados$st=="OK" | Tornados$st=="KS" |Tornados$st=="AR"
                   | Tornados$st=="IA" | Tornados$st=="MO"| Tornados$st=="CO"| Tornados$st=="LA"
                   | Tornados$st=="MN"| Tornados$st=="SD"| Tornados$st=="MS"| Tornados$st=="IL"
                   | Tornados$st=="IN"| Tornados$st=="NE"| Tornados$st=="TN"| Tornados$st=="KY"
                   | Tornados$st=="WI")
year<-rep(begin:end, 17)
alleycounts<-as.data.frame(year)
alleycounts$state <- rep(c("TX", "OK", "KS", "AR", "IA", "MO", "CO", "LA", 
                           "MN","SD", "MS", "IL", "IN", "NE", "TN", "KY",
                           "WI"), times=1, each=66)
tx <- subset(Tornados, Tornados$st=="TX")
tx <- as.integer(table(tx$yr))
ok <- subset(Tornados, Tornados$st=="OK")
ok <- as.integer(table(ok$yr))
ks <- subset(Tornados, Tornados$st=="KS")
ks <- as.integer(table(ks$yr))
ar <- subset(Tornados, Tornados$st=="AR")
ar <- as.integer(table(ar$yr))
ia <- subset(Tornados, Tornados$st=="IA")
ia <- as.integer(table(ia$yr))
mo <- subset(Tornados, Tornados$st=="MO")
mo <- as.integer(table(mo$yr))
co <- subset(Tornados, Tornados$st=="CO")
co <- as.integer(table(co$yr))
la <- subset(Tornados, Tornados$st=="LA")
la <- as.integer(table(la$yr))
mn <- subset(Tornados, Tornados$st=="MN")
mn <- as.integer(table(mn$yr))
sd <- subset(Tornados, Tornados$st=="SD")
sd <- as.integer(table(sd$yr))
ms <- subset(Tornados, Tornados$st=="MS")
ms <- as.integer(table(ms$yr))
il <- subset(Tornados, Tornados$st=="IL")
il <- as.integer(table(il$yr))
ind <- subset(Tornados, Tornados$st=="IN")
ind <- as.integer(table(ind$yr))
ne <- subset(Tornados, Tornados$st=="NE")
ne <- as.integer(table(ne$yr))
tn <- subset(Tornados, Tornados$st=="TN")
tn <- as.integer(table(tn$yr))
tn<-c(tn, 0)
ky <- subset(Tornados, Tornados$st=="KY")
ky <- as.integer(table(ky$yr))
ky<-c(ky, 0)
wi <- subset(Tornados, Tornados$st=="WI")
wi <- as.integer(table(wi$yr))
all<-c(tx, ok, ks, ar, ia, mo, co, la, mn, sd, ms, il, ind, ne, tn, ky, wi )
alleycounts$counts<-all
alleycounts1<-alleycounts[1:529,]
p = ggplot(alleycounts1, aes(x = year, y = counts, group=state)) + 
  geom_line(aes(color=state))+
  geom_point(aes(color=state))
p + ylab("Number of Reported Tornadoes") + 
  labs(title="Number of Reported Tornados in Tornado Alley(part 1)")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=300, label="Source: NOAA, SPC, 2016", size=3)
alleycounts2<-alleycounts[529:1122,]
q = ggplot(alleycounts2, aes(x = year, y = counts, group=state)) + 
  geom_line(aes(color=state))+
  geom_point(aes(color=state))
q + ylab("Number of Reported Tornadoes") + 
  labs(title="Number of Reported Tornados in Tornado Alley (part 2)")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  annotate(geom="text", x=2010, y=150, label="Source: NOAA, SPC, 2016", size=3)

#SUBSET TEXAS
torn_texas<-subset(Tornados, Tornados$st=="TX" &yr>=1970 )
torn_texas$date<-as.Date(torn_texas$date, format="%Y-%m-%d")
Tornados$yr<-as.numeric(Tornados$yr)

##plot per Fscale
TornT = as.data.frame(table(torn_texas$yr, torn_texas[,11]))
TornT$year = as.numeric(levels(TornT$Var1))
TornT$Fscale = paste("F", TornT$Var2, sep = "")
ggplot(TornT[TornT$Var2 != -9, ], aes(x = year, y = Freq)) + 
  geom_point()+ geom_smooth(span = 0.9, color="forestgreen") + 
  facet_wrap(~Fscale, ncol = 2, scales = "free") + 
  theme_gray()+ 
  ggtitle("Number of Reported Tornados per Year \n (Texas) \n 1970-2015 \n by FScale")+
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))+
  ylab("Reported Number of Tornadoes")

#ousiders removal
#zone texas - refs 2949, 5732, 9752, 10216, 10641, 57954 are outside
torn_texas<-subset(torn_texas, 
                   torn_texas$ref != 2949 &
                     torn_texas$ref != 5732 &
                     torn_texas$ref != 9752 &
                     torn_texas$ref != 10216 &
                     torn_texas$ref != 10641 &
                     torn_texas$ref != 57954)

#rm(list=setdiff(ls(), "torn_texas"))
###################point processes
#1. Plotting Point Processes
TornALL <- readOGR(dsn = "./tornado/torn", layer = "torn", stringsAsFactors = FALSE)
CRS.new <- CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0+datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  #EPSG:102003
TornALL <- spTransform(TornALL, CRS.new)
TornALL$Ref<-seq(1, 60114)
TornTexas<-subset(TornALL, TornALL$st=="TX" & 
                    TornALL$Ref != 2949 &
                    TornALL$Ref != 5732 &
                    TornALL$Ref != 9752 &
                    TornALL$Ref != 10216 &
                    TornALL$Ref != 10641 &
                    TornALL$Ref != 57954 & yr>=1970)

x<-TornTexas$slon
y<-TornTexas$slat
coords<-cbind(x,y)
we<-SpatialPoints(coords, CRS("+init=epsg:4326"))
we<-spTransform(we, CRS.new)
coords<-as.matrix(we@coords)
num<-as.numeric(torn_texas$mag)
data1<-as.data.frame(num)
sp1<-SpatialPoints(coords=coords, proj4string=CRS.new)
endTime<-as.POSIXct(31-12-2015, format="%d-%m-%Y", origin = "01-01-2016")
TornTexas$date<-as.Date(TornTexas$date, format="%Y-%m-%d")
asd<-STIDF(sp=sp1, time=TornTexas$date, data=data1)
US.sp <- readOGR(dsn = "./tmp", layer = "cb_2013_us_county_5m", 
                 stringsAsFactors = FALSE)
TX.sp <- US.sp[US.sp$STATEFP == 48, ]
county <- paste(TX.sp$STATEFP, TX.sp$COUNTYFP, sep = "")
county2 <- geometry(spChFIDs(TX.sp, county)) 
counties <- spTransform(county2, CRS.new)
a<-RColorBrewer::brewer.pal(5, "Set3")
wcounty<-unionSpatialPolygons(counties, ID = rep("1", length(row.names(counties))))
years<-c(1970:1992)
w<-list("sp.lines", wcounty, col="lavenderblush4", cex=1.5)
scale = list("SpatialPolygonsRescale", layout.scale.bar(), 
             scale = 400000, 
             fill = c("transparent", "black"), offset = c(-850000, -1200000))
text1 = list("sp.text", c(-850000, -1240000), "0", cex=0.5, col="black")
text2 = list("sp.text", c(-500000, -1240000), "300 km", cex=0.5, col="black")
years2<-c(1992:2015)
stplot(asd[3763:6678,], names.attr=years2, number=24, cuts=5, 
       sp.layout=list(w, scale,text1, text2),
       main="Tornados in Texas per Fujita Scale (1992-2015)", cex=0.7, 
       legendEntries = c("F1", "F2", "F3", "F4", "F5"), key.space="right",
       par.settings =list(panel.background=list(col="gray93"), 
                          strip.background=list(col="lightblue4"),
                          add.text=list(col="white"))) 

stplot(asd[3763:6678,], 2, cuts=5, 
       sp.layout=list(w, scale,text1, text2, arrow),
       main="Tornados in Texas per Fujita Scale (1992-2015)", cex=0.5, 
       legendEntries = c("1", "2", "3", "4", "5"), key.space="right",
       par.settings =list(panel.background=list(col="gray93"), 
                          strip.background=list(col="lightblue4"),
                          add.text=list(col="white"))) 

###exploring intensity

torn_texas<-torn_texas[!duplicated(torn_texas$x),]
a<-min(torn_texas$x)
b<-max(torn_texas$x)
c<-min(torn_texas$y)
d<-max(torn_texas$y)

window_texas<-owin(c(a,b), c(c,d))

x <- torn_texas$x
y <- torn_texas$y
#arcpy script in attachment
poly<-read.csv("points_out_zone2.csv", header=T, sep=";", dec = ",")
poly[,1]<-rev(poly[,1])
poly[,2]<-rev(poly[,2])
poly<-as.matrix(poly)
par(mfrow=c(1,2))
plot(poly[,1], poly[,2])
plot(x, y)
options(scipen=4)
tornado_texas_ppp<-ppp(x, y, window=window_texas, check=T)

#1.1. density calculated w/ standard; Diggle and ppl
a <- density.ppp(tornado_texas_ppp, diggle = T)
b <- density.ppp(tornado_texas_ppp, sigma=bw.diggle, adjust=2)
c <- density.ppp(tornado_texas_ppp, sigma=bw.ppl,adjust=2)
d <- density.ppp(tornado_texas_ppp, sigma = 100000)


par(mar=c(2, 2, 2, 2), mfrow=c(2,4), cex=0.8, oma=c(0, 0, 3, 0))
my_palette <- colorRampPalette(c("black", "white"))(n = 299)
plot(a, main=expression(paste(sigma, " = 150 000")), col=my_palette)
plot(b, main=expression(paste(sigma, " = bw.diggle")), col=my_palette) 
plot(c, main=expression(paste(sigma, " = bw.ppl")), col=my_palette)
plot(d, main=expression(paste(sigma, " = 100 000")), col=my_palette)

persp(a, theta=20, phi=20, zlab="density", border=NA, col="grey", shade=0.75, 
      main=expression(paste(sigma, " = 150 000")))
persp(b, theta=20, phi=20, zlab="density", border=NA, col="grey", shade=0.75, 
      main=expression(paste(sigma, " = bw.diggle (30 221)")))
persp(c, theta=20, phi=20, zlab="density", border=NA, col="grey", shade=0.75, 
      main=expression(paste(sigma, " = bw.ppl (17 554)")))
persp(d, theta=20, phi=20, zlab="density", border=NA, col="grey", shade=0.75, 
      main=expression(paste(sigma, " = 100 000")))
mtext("Spatial Density Study for tornado occurence in Texas", outer = T, cex=1.5)

#If the bandwidth sigma is not specified in the call to density.ppp, the default is to take sigma
#equal to one-eighth of the shortest side length of the enclosing rectangle. This is a very rough rule
#of thumb which may be unsatisfactory in many cases (Badelly, 2016)
###As shown in Figure 6.11, a small value of sigma produces an irregular intensity
##surface, while a large value of sigma appears to oversmooth the intensity.

q <- bw.ppl(tornado_texas_ppp)
w <- bw.diggle(tornado_texas_ppp)
par(mfrow=c(1,2))
plot(q, xlim=c(1,15000), main="Likelihood cross validation method")
plot(w, main=expression(paste("Mean square error at bandwidth ", sigma, sep=" ")))

#1.1.2. Standard deviation of intensity
#Estimate of standard error for the kernel estimate of intensity 
#Uniform edge correction, bandwidth 1 metre

f <- density.ppp(tornado_texas_ppp, sigma = 150000, se=T, diggle=T)$SE
q <- density.ppp(tornado_texas_ppp, sigma = 100000, se=T, diggle=T)$SE
g <- density.ppp(tornado_texas_ppp, sigma = 50000, se=T, diggle = T)$SE

my_palette<-heat.colors(10)
Zlist <- list(a=f, b=q, c=g) 
Zrange <- range(unlist(
  lapply(Zlist, function(x){summary(x)$range})))
plot(as.listof(Zlist), zlim=Zrange, ncols=3, 
     main="Standard Error of Intensities",
     sub="(without covariation with RI)", col=my_palette)


#Density w/ covariate Elevation

RI <- raster("Index_Value1.tif")
RI<-as.im.RasterLayer(RI)

w1<-rhohat(tornado_texas_ppp, RI)

RI2 <- raster("final.tif")
RI2 <-as.im.RasterLayer(RI2)

w2<-rhohat(tornado_texas_ppp, RI2)

par(mfrow=c(1,2))
plot(w1, legend=F, ylab=expression(paste(rho, " (TPI)")), 
     xlab="TPI", main="Intensity as a function of TPI")
plot(w2, legend=F, ylab=expression(paste(rho, " (Elevation)")), 
     xlab="Elevation", main="Intensity as a function of elevation")


j<-predict(w)
plot(j, col=my_palette)
d1 <- eval.im(a - j)
d2 <- eval.im(b - j)
d3 <- eval.im(c - j)
d4 <- eval.im(d - j)

Zlist <- list(dens1=d1, dens2=d2, dens3=d3, dens4=d4)
Zrange <- range(unlist(
  lapply(Zlist, function(x){summary(x)$range})))
plot(as.listof(Zlist), zlim=Zrange, ncols=2, 
     main="Comparison with the original densities",
     sub="(without covariation with RI)", col=my_palette)


#tests for the others covariates
#For practical purposes we provisionally assume the process is Poisson. 
#The null hypothesis is
#that the intensity does not depend on Z, while the alternative is that 
#the intensity does depend on Z
#in some unspecified way
#1) create the spatial dataframe, w/ county code and population
#rm(list=setdiff(ls(), c("TX.sp", "window_texas", "torn_texas", "CRS.new")))
population<-read.csv("Population_final2.csv", header=T, sep=";", dec=".")
TX.sp<-TX.sp[,5]
colnames(population)[1] <- "GEOID"
sp_pop<- merge(TX.sp, population, by = "GEOID")
sp_pop$log1970<-log10(population$pop1970)
sp_pop$log1975<-log10(population$pop1975)
sp_pop$log1980<-log10(population$pop1980)
sp_pop$log1985<-log10(population$pop1985)
sp_pop$log1990<-log10(population$pop1990)
sp_pop$log1995<-log10(population$pop1995)
sp_pop$log2000<-log10(population$pop2000)
sp_pop$log2005<-log10(population$pop2005)
sp_pop$log2010<-log10(population$pop2010)
sp_pop$log2015<-log10(population$pop2015)

sp_pop<-spTransform(sp_pop, CRS.new)

#shp2raster function from
#https://github.com/brry/misc/blob/master/shp2raster.R

popr1970 <- shp2raster(shp=sp_pop, column= "log1970", ascname = "l1970", overwrite=TRUE)
popr1975 <- shp2raster(shp=sp_pop, column= "log1975", ascname = "l1975", overwrite=TRUE)
popr1980 <- shp2raster(shp=sp_pop, column= "log1980", ascname = "l1980", overwrite=TRUE)
popr1985 <- shp2raster(shp=sp_pop, column= "log1985", ascname = "l1985", overwrite=TRUE)
popr1990 <- shp2raster(shp=sp_pop, column= "log1990", ascname = "l1990", overwrite=TRUE)
popr1995 <- shp2raster(shp=sp_pop, column= "log1995", ascname = "l1995", overwrite=TRUE)
popr2000 <- shp2raster(shp=sp_pop, column= "log2000", ascname = "l2000", overwrite=TRUE)
popr2005 <- shp2raster(shp=sp_pop, column= "log2005", ascname = "l2005", overwrite=TRUE)
popr2010 <- shp2raster(shp=sp_pop, column= "log2010", ascname = "l2010", overwrite=TRUE)

#1970
t1970<-subset(torn_texas, torn_texas$yr>=1970 &torn_texas$yr<1976)
t1970<-ppp(t1970$x, t1970$y, window=window_texas, check=T)
den1<-rhohat(t1970, as.im.RasterLayer(popr1970))

#1975
t1975<-subset(torn_texas, torn_texas$yr>=1976 &torn_texas$yr<1980)
t1975<-ppp(t1975$x, t1975$y, window=window_texas, check=T)
den2<-rhohat(t1975, as.im.RasterLayer(popr1975))

#1980
t1980<-subset(torn_texas, torn_texas$yr>=1980 &torn_texas$yr<1986)
t1980<-ppp(t1980$x, t1980$y, window=window_texas, check=T)
den3<-rhohat(t1980, as.im.RasterLayer(popr1980))

#1985
t1985<-subset(torn_texas, torn_texas$yr>=1986 &torn_texas$yr<1990)
t1985<-ppp(t1985$x, t1985$y, window=window_texas, check=T)
den4<-rhohat(t1985, as.im.RasterLayer(popr1985))

#1990
t1990<-subset(torn_texas, torn_texas$yr>=1990 &torn_texas$yr<1996)
t1990<-ppp(t1990$x, t1990$y, window=window_texas, check=T)
den5<-rhohat(t1990, as.im.RasterLayer(popr1990))

#1995
t1995<-subset(torn_texas, torn_texas$yr>=1996 &torn_texas$yr<2000)
t1995<-ppp(t1995$x, t1995$y, window=window_texas, check=T)
den6<-rhohat(t1995, as.im.RasterLayer(popr1995))

#2000
t2000<-subset(torn_texas, torn_texas$yr>=2000 &torn_texas$yr<2006)
t2000<-ppp(t2000$x, t2000$y, window=window_texas, check=T)
den7<-rhohat(t2000, as.im.RasterLayer(popr2000))

#2000
t2005<-subset(torn_texas, torn_texas$yr>=2006 &torn_texas$yr<2010)
t2005<-ppp(t2005$x, t2005$y, window=window_texas, check=T)
den8<-rhohat(t2005, as.im.RasterLayer(popr2005))

#2010
t2010<-subset(torn_texas, torn_texas$yr>=2010 &torn_texas$yr<2016)
t2010<-ppp(t2010$x, t2010$y, window=window_texas, check=T)
den9<-rhohat(t2010, as.im.RasterLayer(popr2010))


q1<-(ggplot(den1, aes(x=X, y=rho))+geom_line(aes(X, rho))+
  geom_ribbon(data=den1,aes(ymin=hi,ymax=lo),alpha=0.3))+ 
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="1970-1975")+ theme(axis.title.y = element_text(size = 9),
                                 axis.title.x = element_text(size = 9),
                                 plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q2<-(ggplot(den2, aes(x=X, y=rho))+geom_line(aes(X, rho))+
  geom_ribbon(data=den2,aes(ymin=hi,ymax=lo),alpha=0.3))+
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="1975-1980")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q3<-(ggplot(den3, aes(x=X, y=rho))+geom_line(aes(X, rho))+
  geom_ribbon(data=den3,aes(ymin=hi,ymax=lo),alpha=0.3))+ 
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="1980-1985")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q4<-(ggplot(den4, aes(x=X, y=rho))+geom_line(aes(X, rho))+
  geom_ribbon(data=den4,aes(ymin=hi,ymax=lo),alpha=0.3))+ ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="1985-1990")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))

q5<-(ggplot(den5, aes(x=X, y=rho))+geom_line(aes(X, rho))+
  geom_ribbon(data=den5,aes(ymin=hi,ymax=lo),alpha=0.3))+
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="1190-1995")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q6<-(ggplot(den6, aes(x=X, y=rho))+geom_line(aes(X, rho))+
       geom_ribbon(data=den6,aes(ymin=hi,ymax=lo),alpha=0.3))+
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="1995-2000")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q7<-(ggplot(den7, aes(x=X, y=rho))+geom_line(aes(X, rho))+
       geom_ribbon(data=den7,aes(ymin=hi,ymax=lo),alpha=0.3))+
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="2000-2005")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q8<-(ggplot(den8, aes(x=X, y=rho))+geom_line(aes(X, rho))+
       geom_ribbon(data=den8,aes(ymin=hi,ymax=lo),alpha=0.3))+
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="2005-2010")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))
q9<-(ggplot(den9, aes(x=X, y=rho))+geom_line(aes(X, rho))+
       geom_ribbon(data=den9,aes(ymin=hi,ymax=lo),alpha=0.3))+
  ylab(expression(paste(rho, "(X)"))) + xlab("log(Pop)") +
  labs(title="2010-2015")+
  theme(axis.title.y = element_text(size = 9),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size=9, face="bold", 
                                  margin = margin(10, 0, 10, 0), hjust=0.5))

multiplot(q1, q2, q3, q4, q5, q6, q7, q8, q9, cols=3)
multiplot(q1,q4 ,q7 , q2, q5,q8, q3,q6,q9, cols=3 )
grid.arrange(q1, q2, q3, q4, q5, q6, q7, q8, q9,
             top = title1)
title1=textGrob("Intensity as a function of Population", gp=gpar(fontsize=15, fontface="bold"))


######INHOMOGENEOUS k-Function
######note: maybe do for all years?? to see if it is a poisson ?

kinhomtexas<-envelope(tornado_texas_ppp, Kinhom, nsim=99, 
                      simulate = expression(rpoispp(a)))
kinhomtexas2<-envelope(tornado_texas_ppp, Kinhom, nsim=99, 
                       simulate = expression(rpoispp(d)))
par(mar=c(2,4,2,4))
v<-plot(kinhomtexas, legend=FALSE, main="Inhomogeneous K-function")
legend(-20000, 190022224257, legend=v$meaning, lty=v$lty, col=v$col, cex=0.6, bty="n")

#rm(list=setdiff(ls(), c("torn_texas", "poly")))
plotK()
####space-time
## Space-time inhomogeneous K-function

data<-torn_texas[ ,c(24, 25, 5)]
TX <- as.3dpoints(data[,1]/1000, data[,2]/1000, data[,3])
Poly <- poly/1000

# Estimation of the temporal intensity

Mt <- density(TX[ ,3], n = 1000)
mut <- Mt$y[findInterval(TX[ ,3], Mt$x)] * dim(TX)[1]

# Estimation of the spatial intensity
# Finding the optimal bandwidth for kernel smoothing

h <- mse2d(as.points(TX[,1:2]), Poly, nsmse = 100, range = 4)
h <- h$h[which.min(h$mse)]
Ms <- kernel2d(as.points(TX[ ,1:2]), Poly, h = h, nx = 500, ny = 500)
atx <- findInterval(x = TX[ ,1], vec = Ms$x)
aty <- findInterval(x = TX[ ,2], vec = Ms$y)
mhat <- NULL
for(i in 1:length(atx)) mhat <- c(mhat, Ms$z[atx[i],aty[i]])

# Estimation of the STIK function
#for dx=400km, dt=100days
u <- seq(0,400, leng = 20)
v <- seq(0,100, leng= 20)
stik <- STIKhat(xyt = TX, s.region = Poly, t.region = c(1,24056),
                lambda = mhat*mut/7355, dist = u, times = v, infectious = T)

plotK(stik, L=FALSE,type="persp", theta = 30, phi = 20, legend=TRUE)
plotK(stik, L=TRUE, type="persp",theta=30, phi = 30,legend=T)

#for dx=400 km, and dt=4years
u <- seq(0,400, leng = 20)
v <- seq(0,1460, leng= 20)
stik2 <- STIKhat(xyt = TX, s.region = Poly, t.region = c(1,24056),
                 lambda = mhat*mut/7355, dist = u, times = v, infectious = T)

plotK(stik2, L=FALSE,type="persp", theta = 30, phi = 20, legend=TRUE)
plotK(stik2, L=TRUE, type="persp",theta=30, phi = 30,legend=T)

#for dx=100 km, and dt=4years
u <- seq(0,100, leng = 20)
v <- seq(0,1460, leng= 20)
stik3 <- STIKhat(xyt = TX, s.region = Poly, t.region = c(1,24056),
                 lambda = mhat*mut/7355, dist = u, times = v, infectious = T)
plotK(stik3, L=FALSE,type="persp", theta = 30, phi = 20, legend=TRUE)
plotK(stik3, L=TRUE, type="persp",theta=30, phi = 30,legend=T)

#for dx=100 km, and dt=4years
u <- seq(0,100, leng = 20)
v <- seq(0,1460, leng= 20)
stik3 <- STIKhat(xyt = TX, s.region = Poly, t.region = c(1,24056),
                 lambda = mhat*mut/7355, dist = u, times = v, infectious = T)

plotK(stik4, L=FALSE,type="persp", theta = 20, phi = 20, legend=TRUE)
plotK(stik4, L=TRUE, type="persp",theta=30, phi = 30,legend=T)

#simulation of point pattern for Poisson
#cluster 
simulation<-rpcp(s.region=Poly, t.region=c(1,24056), npoints=60114, lambda=mhat*mut/7355, 
     mc=NULL, nsim=1, cluster="uniform", infectious=TRUE, 
     edge = "without", larger.region=larger.region, tronc=1)

#Poisson inhomogeneous
simulation.poiss<-rpp(lambda=mhat*mut/7355, s.region =Poly, t.region=c(1,24056), npoints=60144, 
    nsim=10, replace=TRUE,
    discrete.time=FALSE, nx=100, ny=100, nt=100, lmax=NULL)

#simulate a log-gaussian cox process
simluation.cox<-rlgcp(s.region=Poly, t.region=c(1,24056), replace=TRUE, 
                      npoints=60114, nsim=10, 
      nx=100, ny=100, nt=100,separable=TRUE,model="exponential",
      param=c(1,1,1,1,1,1), scale=c(1,1),var.grf=1,mean.grf=0,
      lmax=NULL,discrete.time=FALSE,exact=FALSE)
# 