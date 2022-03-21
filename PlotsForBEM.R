# ======================= #
# PLOTS FOR PUBLICATION   #
# ======================= #
# revised December 2021, Carmen L David
#
# Includes:
#    - plots for manuscript after revision
#    - plots for supplementary material
#    - comments on data used in plots

setwd("/.../")
source("BioenergeticIndModel.R")

require(ggplot2) 
require(ggpubr)
require(reshape2)
library(ggOceanMaps)

# ======================= #
# MANUSCRIPT              #
# ======================= #


# ===================== #
# Figure 1 #####
# ===================== #

# Figure 1 Map
all_poly <- list(
  Baff = data.frame(lon = c(-82, -83, -60,-60), 
                    lat = c(70, 78, 78, 70),
                    Region = rep("Baffin Bay", 4)),
  Beaufort = data.frame(lon = c(-135, -135, -115, -115), 
                        lat = c(69, 72, 72, 69),
                        Region = rep("Beaufort Sea", 4)),
  Kitik = data.frame(lon = c(-105, -112,-82.5,-83.5), 
                     lat = c(67, 73, 74, 67),
                     Region = rep("Kitikmeot", 4)),
  Laptev = data.frame(lon = c(160, 160,110, 110), 
                      lat = c(69, 78,78,69),
                      Region = rep("Laptev Sea", 4)))

df_poly <- do.call(rbind, all_poly)
Region = df_poly$Region
coord_labs=data.frame(labels = c("0º", "90º", "180º", "270º", "55º", "65º",  "75º"), 
                      lon=c(0,92,180,-92, 131, 131, 131),  
                      lat=c(57, 57, 57, 57, 56, 66, 76),
                      angle = c(0, 0, 0, 0, 310, 310, 310))
count_labs =data.frame(labels = c("Greenland", "Canada", "USA", "Russia", "Norway"), 
                       lon=c(-40, -120, -145, 115, 25),  
                       lat=c(70, 64, 65, 65, 67))

basemap(55, bathymetry = T,  legends = T) +
  geom_text(data = transform_coord(count_labs[,2:3]), aes(x = lon, y = lat, angle=0), 
            color = "grey80",label= count_labs$labels) +
  geom_text(data = transform_coord(coord_labs[,2:3]), aes(x = lon, y = lat, angle=coord_labs$angle), 
            color = "grey80",label= coord_labs$labels) +
  geom_polygon(data = transform_coord(df_poly[,1:2]), aes(x = lon, y = lat, col = Region), size = 1, fill = NA)

dev.copy2pdf(file="Figure1.pdf", width=6.5, height=5.5)

# =================================== #



# ===================== #
# Figure 2 #####
# ===================== #
# estimate fraction of nauplii
# use an alternative BEM, in which npl is fixed
source("BEM_testing/BioenergeticIndModel.R")

age <- seq(1, 60, 1)

# one larva
larvamm = 8.5 #mm
larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
df1 <- data.frame()
for (x in 1:length(age)){ 
  oneday <- BEM(larvawgt, larvamm, temp=-1, P=0.65, npl = 0.6)
  larvawgt = oneday$larvawgt
  larvamm = oneday$larvamm
  df1 <- rbind(df1, oneday)
}
larvamm = 8.5 #mm
larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
df2 <- data.frame()
for (x in 1:length(age)){ 
  oneday <- BEM(larvawgt, larvamm, temp=-1, P=0.65, npl = 0.8)
  larvawgt = oneday$larvawgt
  larvamm = oneday$larvamm
  df2 <- rbind(df2, oneday)
}
larvamm = 8.5 #mm
larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
df3 <- data.frame()
for (x in 1:length(age)){ 
  oneday <- BEM(larvawgt, larvamm, temp=-1, P=0.65, npl = 1)
  larvawgt = oneday$larvawgt
  larvamm = oneday$larvamm
  df3 <- rbind(df3, oneday)
}

# read the clasic BEM
source("BioenergeticIndModel.R")
larvamm = 8.5 #mm
larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
df4 <- data.frame()
for (x in 1:length(age)){ 
  oneday <- BEM(larvawgt, larvamm, temp=-1, P=0.65)
  larvawgt = oneday$larvawgt
  larvamm = oneday$larvamm
  df4 <- rbind(df4, oneday)
}

df_reg <- data.frame(Age = age,
                     "eq." = df4$larvamm,
                     "f1.0" = df3$larvamm,
                     "f0.8" = df2$larvamm,
                     "f0.6" = df1$larvamm)


df_long <- melt(df_reg, id = "Age")  # convert to long format

main.plot <- ggplot(df_long, aes(y=value, x=Age, linetype=variable)) + geom_line(size=0.6) + 
  labs(y = "Length (mm)", x = "Age (days)", linetype = "") + 
  scale_x_continuous(breaks = seq(0,60,10)) +
  scale_y_continuous(breaks = seq(5,20,2)) +
  scale_linetype(name = "Fraction nauplii", labels = c("eq.","1.0", "0.8", "0.6")) + 
  theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA)) +
  theme(legend.position = c(0.91,0.26)) +
  theme(legend.key.size = unit(2, "lines"))+
  theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11))


# inset #
larvamm <- seq(4.5, 12, length.out = 200)
# equation from Thannassekos 2012 for larvamm<12mm
fn = -0.0042*larvamm^2 + 0.1055*larvamm + 0.1961 

# add fraction nauplii fn=1 for larvae > 25 mm, after Bouchard et al. 2020
larvamm <- c(larvamm, seq(15, 25, length.out = 100))
fn <- c(fn, rep(1, 100))
df <- data.frame(fn=fn, L=larvamm, npl = 1-exp(-0.19786*(larvamm-0.34)))

npl = 1-exp(-0.19786*(larvamm-0.34))

plot(df$L, df$fn, xlab = "Length (mm)", ylab = "Fraction nauplii", pch=20, cex.lab=1.3, cex.axis=1.3)
lines(df$L, df$npl, col="red", lwd=4)

inset.plot <-  ggplot(df, aes(y=fn, x=L)) + geom_point(size=0.6) + 
  labs(y = "Fraction nauplii", x = "Length (mm)", linetype = "") +
  geom_line(data = df, aes(y = npl, x = L, color = "red"), size = 1.2) +
  theme_classic()+ 
  theme(legend.position = c(0.57,0.18), legend.title = element_blank(), legend.text = element_text(size = 11.5))+
  scale_color_hue(labels = expression(paste("f"["nauplii"], " = 1 - e"^" 0.198 · (SL - 0.34)")))


library(cowplot)
ggdraw() +
  draw_plot(main.plot) +
  draw_plot(inset.plot, x = 0.065, y = .53, width = .46, height = .45)

dev.copy2pdf(file="Figure2.pdf", width=8, height=5)


# ========================== #
# Figure 3 ####
# ========================== #
# analysis in: EstimateParam_Validation.R
# includes more plots for Suplementary material
#
# select data from regions and years
# download mean SST in each region and year, save them in df_temp
# fit length-age regressions and add them in df_age
#
require("reshape2")
require(ggplot2)
age <- seq(1, 200, 1) # in grams
#
# Region-years for fitting (Table 1)
# SST data from Table 1
df_temp <- data.frame(Laptev2003 = -0.97,
                      Laptev2007 = -0.64,
                      Beauford2010 = 1.78,
                      Beafort2011 = 1.21,
                      Beauford2014 = 0.30,
                      Kitikmeot2006 = -1.23,
                      Kitikmeot2011 = -0.98,
                      Kitikmeot2015 = -1.39,
                      Baffin2006 = -1.21,
                      Baffin2008 = -1.10)
# length-at-age from otoliths measurements per region-year (Bouchard and Fortier 2011, 2020)
# Individual regressions fitted on data
df_age <- data.frame(age = age,
                     Laptev2003 = 2.840209 + 0.207227*age,
                     Laptev2007 = 8.510625 + 0.189755*age,
                     Beaufort2010 = 5.629976 + 0.171524*age,
                     Beaufort2011 = 2.469326 + 0.198197*age,
                     Beaufort2014 = 8.92386 + 0.17543*age,
                     Kitikmeot2006 = 5.04346 + 0.20379*age,
                     Kitikmeot2011 = 0.955450 + 0.231559*age,
                     Kitikmeot2015 = 2.788851 + 0.199910*age,
                     Baffin2006 = 6.55684 + 0.19465*age,
                     Baffin2008 = 1.776091 + 0.209800*age)


#BaffinBay
ndays = 200
temp1 = apply(df_temp[, grep("Baffin", names(df_temp), value=TRUE)], 1, mean)

larvamm = 4.16
larvawgt = 0.0055 * ((larvamm/10)^3.19) # mm
df = data.frame()#age = 0, W = larvawgt, L = larvamm)
for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp=temp1, P=0.85)
  larvawgt = oneDay$larvawgt
  larvamm = oneDay$larvamm
  df = rbind(df, data.frame(age = i, W = larvawgt, L = larvamm))
}

# initial length
larvamm = 5.5
larvawgt = 0.0055 * ((larvamm/10)^3.19) # mm
df2 = data.frame()#age = 0, W = larvawgt, L = larvamm)
for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp=temp1, P=0.85)
  larvawgt <- oneDay$larvawgt
  larvamm <- oneDay$larvamm
  df2 = rbind(df2, data.frame(age = i, W = larvawgt, L = larvamm))
}

df3 <- data.frame(age = df_age$age,
                  P85=df[,"L"],
                  df_age[,grep('Baffin', names(df_age), value=T)])

df_long <- melt(df3, id = "age")  # convert to long format

df <-  read.table("Datasets/DatabaseLengthAge_Baffin.txt", header=T, sep="\t") 
df <- df[grep("2006|2008", df$Year), ]
df <- df[!is.na(df$Age),]
df_points <- rbind(data.frame(age = df$Age, Year=as.factor(df$Year), variable = "Obs", value = df$LS))

p1 <- ggplot() + geom_point(data = df_points, aes(y=value, x=age, color=Year))+
  ylim(0, 45) + xlim(0,200)+
  geom_line(data=df_long, aes(y=value, x=age, linetype=variable))+
  labs(y = "Length (mm)", x = "Age (days)", linetype = "")+
  annotate(geom="text", x=30, y=40, label="Baffin Bay", cex = 4.5) +
  annotate(geom="text", x=30, y=36, label="T = -1.15ºC") +
  annotate(geom="text", x=30, y=34, label="L_ini = 4.16 mm") + 
  theme_bw() +
  theme(legend.title=element_blank(), axis.title.x = element_blank())+
  guides(color = guide_legend(order = 2), 
         linetype = guide_legend(order = 1))



#Beaufort
ndays = 200
temp1 = apply(df_temp[, grep("Beauf", names(df_temp), value=TRUE)], 1, mean)

larvamm = mean(c(5.629976,2.469326,8.92386))
larvawgt = 0.0055 * ((larvamm/10)^3.19)  # mm
df = data.frame()#age = 0, W = larvawgt, L = larvamm)
for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp=temp1, P=0.6)
  larvawgt <- oneDay$larvawgt
  larvamm <- oneDay$larvamm
  df = rbind(df, data.frame(age = i, W = larvawgt, L = larvamm))
}

df3 <- data.frame(age = df_age$age, 
                  P60=df[,"L"],
                  df_age[,grep('Beauf', names(df_age), value=T)])
df_long <- melt(df3, id = "age")  # convert to long format

df <-  read.table("Datasets/DatabaseLengthAge_Beaufort.txt", header=T, sep="\t") 
df <- df[grep("2010|2011|2014", df$Year), ]
df <- df[!is.na(df$Age),]
df_points <- rbind(data.frame(age = df$Age, Year=as.factor(df$Year), variable = "Obs", value = df$LS))

p2 <- ggplot() + geom_point(data = df_points, aes(y=value, x=age, color=Year))+
  ylim(0, 45) + xlim(0,200)+
  geom_line(data=df_long, aes(y=value, x=age, linetype=variable))+
  labs(y = "Length (mm)", x = "Age (days)", linetype = "")+
  annotate(geom="text", x=30, y=40, label="Beaufort Sea", cex = 4.5) +
  annotate(geom="text", x=30, y=36, label="T = 1.04 ºC") +
  annotate(geom="text", x=30, y=34, label="L_ini = 5.67 mm") +
  theme_bw() +
  theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())+
  guides(color = guide_legend(order = 2), 
         linetype = guide_legend(order = 1))



# Kitikmeot
ndays = 200
temp1=apply(df_temp[, grep("Kitik", names(df_temp), value=TRUE)], 1, mean)

larvamm = mean(c(5.04346,0.955450,2.788851))
larvawgt = 0.0055 * ((larvamm/10)^3.19) # mm
df = data.frame()#age = 0, W = larvawgt, L = larvamm)
for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp=temp1, P=0.82)
  larvawgt <- oneDay$larvawgt
  larvamm <- oneDay$larvamm
  df = rbind(df, data.frame(age = i, W = larvawgt, L = larvamm))
}


df3 <- data.frame(age = df_age$age, 
                  P82=df[,"L"],
                  df_age[,grep('Kitik', names(df_age), value=T)])  

df_long <- melt(df3, id = "age")  # convert to long format

df <-  read.table("Datasets/DatabaseLengthAge_Kitikmeot.txt", header=T, sep="\t") 
df <- df[grep("2006|2011|2015", df$Year), ]
df <- df[!is.na(df$Age),]
df_points <- rbind(data.frame(age = df$Age, Year=as.factor(df$Year), variable = "Obs", value = df$LS))

p3 <- ggplot() + geom_point(data = df_points, aes(y=value, x=age, color=Year))+
  ylim(0, 45) + xlim(0,200)+
  geom_line(data=df_long, aes(y=value, x=age, linetype=variable))+
  labs(y = "Length (mm)", x = "Age (days)", linetype = "")+
  annotate(geom="text", x=30, y=40, label="Kitikmeot", cex = 4.5) +
  annotate(geom="text", x=30, y=36, label="T = -1.20 ºC") +
  annotate(geom="text", x=30, y=34, label="L_ini = 2.93 mm") + 
  theme_bw() +
  theme(legend.title=element_blank()) +
  guides(color = guide_legend(order = 2), 
         linetype = guide_legend(order = 1))





# Laptev sea
ndays = 200
temp1=apply(df_temp[, grep("Laptev", names(df_temp), value=TRUE)], 1, mean)

larvamm = 5.67
larvawgt = 0.0055 * ((larvamm/10)^3.19) # mm
df = data.frame()#age = 0, W = larvawgt, L = larvamm)
for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp=temp1, P=0.71)
  larvawgt <- oneDay$larvawgt
  larvamm <- oneDay$larvamm
  df = rbind(df, data.frame(age = i, W = larvawgt, L = larvamm))
}

df3 <- data.frame(age = df_age$age, 
                  P71=df[,"L"],
                  df_age[,grep('Laptev', names(df_age), value=T)]) 

df_long <- melt(df3, id = "age")  # convert to long format

df <-  read.table("Datasets/DatabaseLengthAge_Laptev.txt", header=T, sep="\t") 
df <- df[grep("2003|2007", df$Year), ]
df <- df[!is.na(df$Age),]
df_points <- rbind(data.frame(age = df$Age, Year=as.factor(df$Year), variable = "Obs", value = df$LS))
  
p4 <- ggplot() + geom_point(data = df_points, aes(y=value, x=age, color=Year))+
  ylim(0, 45) + xlim(0,200)+
  geom_line(data=df_long, aes(y=value, x=age, linetype=variable))+
  labs(y = "Length (mm)", x = "Age (days)", linetype = "")+
  annotate(geom="text", x=30, y=40, label="Laptev Sea", cex = 4.5) +
  annotate(geom="text", x=30, y=36, label="T = -0.81 ºC") +
  annotate(geom="text", x=30, y=34, label="L_ini = 5.67 mm") + 
  theme_bw() +
  theme(legend.title=element_blank(), axis.title.y = element_blank())+
  guides(color = guide_legend(order = 2), 
         linetype = guide_legend(order = 1))




ggpubr::ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2)

dev.copy2pdf(file="Figure3.pdf", width=10, height=7)


# ========================== #
# Figure 4 ####
# ========================== # 

# BEM with multiple hatch dates
source("BioenergeticIndModel.R")
source("extractSST.R")
source("extractDates.R")
require(ncdf4)
require(ggplot2)

#initial larval length for all simulations
larvamm_ini = 6

#Baffin# 
day_temp <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2005-02-15", date_end = "2005-09-15")
DTS <- extractDates(file_sst = "SST_downloads/Baffin.nc", date_start = "2005-02-15", date_end = "2005-06-15")
hatchdates <- c("2005-02-15", "2005-03-15", "2005-04-15", "2005-05-15", "2005-06-15")
hatchID <- match(hatchdates, DTS)

P_val = 0.85
# BEM #
df_hatchers <- data.frame()
for (j in 1:length(hatchID)){ 
  larvamm = larvamm_ini #mm
  larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
  df <- data.frame()
  for (i in hatchID[j]:dim(day_temp)[3]){
    newT = mean(day_temp[,,i], na.rm = T)
    one_larva <- BEM(larvawgt, larvamm, temp = newT, P = P_val)
    larvawgt = one_larva$larvawgt
    larvamm = one_larva$larvamm
    df_day = data.frame(Age = i-hatchID[j]+1, Length = larvamm)
    df <- rbind(df, df_day)
  }
  df$hatchdate <- hatchdates[j]
  df_hatchers <- rbind(df_hatchers, df)
}

df <-  read.table("Datasets/DatabaseLengthAge_Baffin.txt", header=T, sep="\t") 
df <- df[grep("2005", df$Year), ]
df_all <- rbind(data.frame(Age = df$Age, Length = df$LS, Origin = "Observations"),
                data.frame(Age = df_hatchers$Age, Length = df_hatchers$Length, Origin = df_hatchers$hatchdate))


p1 <- ggplot(df_all, aes(y=Length, x=Age, color = Origin))  + geom_point(size=0.7) + xlim(0, 250)+ ylim(0, 60) +#geom_line(aes(linetype=Origin))+
  geom_text(x=40, y=55, label="Baffin Bay (n = 71) \n P = 0.85 \n EF = 0.55", col = "black", size = 3) +
  labs(y = "Length (mm)", x = "Age (days)", linetype = "") + theme_bw() +
  theme(legend.title=element_blank(), axis.title.x = element_blank()) 


#Kitikmeot#
day_temp <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2005-02-15", date_end = "2005-09-15")
DTS <- extractDates(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2005-02-15", date_end = "2005-09-15")
hatchdates <- c("2005-02-15", "2005-03-15", "2005-04-15", "2005-05-15", "2005-06-15")
hatchID <- match(hatchdates, DTS)

P_val = 0.82
# BEM #
df_hatchers <- data.frame()
for (j in 1:length(hatchID)){ 
  larvamm = larvamm_ini #mm
  larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
  df <- data.frame()
  for (i in hatchID[j]:dim(day_temp)[3]){
    newT = mean(day_temp[,,i], na.rm = T)
    one_larva <- BEM(larvawgt, larvamm, temp = newT, P = P_val)
    larvawgt = one_larva$larvawgt
    larvamm = one_larva$larvamm
    df_day = data.frame(Age = i-hatchID[j]+1, Length = larvamm)
    df <- rbind(df, df_day)
  }
  df$hatchdate <- hatchdates[j]
  df_hatchers <- rbind(df_hatchers, df)
}

df <-  read.table("Datasets/DatabaseLengthAge_Kitikmeot.txt", header=T, sep="\t")
df <- df[grep("2005", df$Year), ]
df_all <- rbind(data.frame(Age = df$Age, Length = df$LS, Origin = df$Region),
                data.frame(Age = df_hatchers$Age, Length = df_hatchers$Length, Origin = df_hatchers$hatchdate))

p3 <- ggplot(df_all, aes(y=Length, x=Age, color = Origin))  + geom_point(size=0.7) + xlim(0, 250)+ ylim(0, 60) +#geom_line(aes(linetype=Region))+
  geom_text(x=40, y=55, label="Kitikmeot (n = 32) \n P = 0.82\n EF = 0.94", col = "black", size = 3) +
  labs(y = "Length (mm)", x = "Age (days)", linetype = "") + theme_bw() +
  theme(legend.title=element_blank()) 



#Beaufort#
day_temp <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2005-02-15", date_end = "2005-09-15")
DTS <- extractDates(file_sst = "SST_downloads/Beaufort.nc", date_start = "2005-02-15", date_end = "2005-09-15")
hatchdates <- c("2005-02-15", "2005-03-15", "2005-04-15", "2005-05-15", "2005-06-15")
hatchID <- match(hatchdates, DTS)

P_val = 0.6
# BEM #
df_hatchers <- data.frame()
for (j in 1:length(hatchID)){ 
  larvamm = larvamm_ini #mm
  larvawgt = 0.0055 * ((larvamm/10)^3.19) #grams
  df <- data.frame()
  for (i in hatchID[j]:dim(day_temp)[3]){
    newT = mean(day_temp[,,i], na.rm = T)
    one_larva <- BEM(larvawgt, larvamm, temp = newT, P = P_val)
    larvawgt = one_larva$larvawgt
    larvamm = one_larva$larvamm
    df_day = data.frame(Age = i-hatchID[j]+1, Length = larvamm)
    df <- rbind(df, df_day)
  }
  df$hatchdate <- hatchdates[j]
  df_hatchers <- rbind(df_hatchers, df)
}

df <-  read.table("Datasets/DatabaseLengthAge_Beaufort20056.txt", header=T, sep="\t")
df <- df[grep("2005", df$Year), ]
df <- df[!is.na(df$Age),]
df_all <- rbind(data.frame(Age = df$Age, Length = df$LS, Origin = df$Region),
                data.frame(Age = df_hatchers$Age, Length = df_hatchers$Length, Origin = df_hatchers$hatchdate))

p5 <- ggplot(df_all, aes(y=Length, x=Age, color = Origin))  + geom_point(size=0.7) + xlim(0, 250)+ ylim(0, 60) +#geom_line(aes(linetype=Region))+
  geom_text(x=40, y=55, label="Beaufort Sea (n = 38) \n P = 0.61 \n EF = 0", col = "black", size = 3) +
  labs(y = "Length (mm)", x = "Age (days)", linetype = "") + theme_bw() +
  theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank()) 


#Laptev#
day_temp <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2005-02-15", date_end = "2005-09-15")
DTS <- extractDates(file_sst = "SST_downloads/Laptev.nc", date_start = "2005-02-15", date_end = "2005-09-15")
hatchdates <- c("2005-02-15", "2005-03-15", "2005-04-15", "2005-05-15", "2005-06-15")
hatchID <- match(hatchdates, DTS)
P_val = 0.71
# BEM #
df_hatchers <- data.frame()
for (j in 1:length(hatchID)){ 
  larvamm = larvamm_ini #mm
  larvawgt = 0.0055 * ((larvamm/10)^3.19) #grams
  df <- data.frame()
  for (i in hatchID[j]:dim(day_temp)[3]){
    newT = mean(day_temp[,,i], na.rm = T)
    one_larva <- BEM(larvawgt, larvamm, temp = newT, P = P_val)
    larvawgt = one_larva$larvawgt
    larvamm = one_larva$larvamm
    df_day = data.frame(Age = i-hatchID[j]+1, Length = larvamm)
    df <- rbind(df, df_day)
  }
  df$hatchdate <- hatchdates[j]
  df_hatchers <- rbind(df_hatchers, df)
}


df <-  read.table("Datasets/DatabaseLengthAge_Laptev.txt", header=T, sep="\t") 
df <- df[grep("2005", df$Year), ]
df <- df[!is.na(df$Age),]
df_all <- rbind(data.frame(Age = df$Age, Length = df$LS, Origin = "Laptev"),
                data.frame(Age = df_hatchers$Age, Length = df_hatchers$Length, Origin = df_hatchers$hatchdate))

p7 <- ggplot(df_all, aes(y=Length, x=Age, color = Origin))  + geom_point(size=0.7) + xlim(0, 250)+ ylim(0, 60) + #geom_line(aes(linetype=Region))+
  geom_text(x=40, y=55, label="Laptev Sea (n = 196) \n P = 0.71 \n EF = 0.69", col = "black", size = 3) + theme_bw() +
  labs(y = "Length (mm)", x = "Age (days)", linetype = "") + 
  theme(legend.title=element_blank(), axis.title.y = element_blank()) 




ggpubr::ggarrange(p1,p5,p3,p7, nrow = 2, ncol = 2, widths=c(2,2), common.legend = T)

dev.copy2pdf(file="Figure4dec.pdf", width=8, height=6)



# ========================== #
# Figure 5 ####
# ========================== #
# analysis and dataset produced in: SensitivityAnalysis.R
# here new plot
#
df <- read.table("Datasets/SensitivityAnalysis_Sep2021.txt")
# rm SDA
df <- df[df$Parameter!= "SDA",]

df2 <- dplyr::filter(df, Test=="+20%" | Parameter=="Ref")
df2$group <- as.factor(c(rep("Ref",31), rep("Gr1", 155)))

require(ggplot2)  
p1 <-   ggplot(df2, aes(x=Age, y=L, color = reorder(Parameter, -L))) + 
  geom_line(size = 1.1) + 
  labs(x = "Age (days)", y = "Length (mm)", color="Parameter") +
  #scale_linetype_manual("Effect", labels = c("positive", "negative"), values=c("solid", "dashed")) +
  scale_x_continuous(breaks=c(seq(0,30,5))) +
  scale_y_continuous(breaks=c(seq(6,14,1))) +
  theme_bw()+ 
  theme(text = element_text()) +
  geom_text(x=2.5, y=14, label="A", color = "black", cex = 5)


df3 <- dplyr::filter(df, Age==30)
df3$deltaL <- (df3$L - df3[df3$Parameter=="Ref", "L"]) / df3[df3$Parameter=="Ref", "L"]
df3 <- df3[df3$Parameter != "Ref",]


library(forcats)
p2 <- ggplot(df3, aes(x=reorder(Parameter,deltaL, min), y=deltaL, fill = Test)) + 
  geom_bar(stat="identity", width=0.4) + 
  labs(title = "", x = "Parameter", y = "Ratio of change") + scale_fill_grey() +
  theme(text = element_text(size = 16), legend.title=element_blank()) +
  theme_bw()+
  geom_text(x=1, y=0.08, label="B", color = "black", cex = 5)

ggpubr::ggarrange(p1,p2, ncol = 1, nrow = 2, common.legend = F)

dev.copy2pdf(file="Figure%.pdf", width=5.5, height=7)



# ========================== #
# Figure 6 ####
# ========================== #
# analysis and dataset produced in: SensitivityAnalysis.R
# plot produced in: SensitivityAnalysis.R


# ========================== #
# Figure 7 ####
# ========================== #
#
#  over a range of temperatures
tempr <- seq(-1.5, 10, 0.1) # in grams
# P_val takes values +/-20%
P_val = 0.73
l_ini = 12 #mm
w_ini = 0.0055 * ((l_ini/10)^3.19) 
df <- data.frame()
for (x in 1:length(tempr)){ 
  oneday <- BEM(larvawgt=w_ini, larvamm=l_ini, temp=tempr[x], P=P_val)
  df <- rbind(df, oneday)
}

# experimental growth from Koenker et al.2018
# respiration uses a Q10 expression
source("BEM_testing/BioenergeticInvers.R")
l_ini = 12 #mm
w_ini = 0.0055 * ((l_ini/10)^3.19) 
df4 <- data.frame()
for (x in 1:length(tempr)){ 
  oneday <- BEMLab(larvawgt=w_ini, larvamm=l_ini, temp=tempr[x], P=P_val)
  df4 <- rbind(df4, oneday)
}

df_reg <- data.frame(Temp = tempr,
                     GR = df$GR,
                     C = df$Cmax*0.73*0.8, # real consumption, feeding efficiency included
                     R = df$R,
                     Lab = df4$GR)
df_long <- melt(df_reg, id = "Temp")  # convert to long format

# growth rates with
tempr <- seq(-1.5, 10, 0.1) # in grams
l_ini = 12 #mm
w_ini = 0.0055 * ((l_ini/10)^3.19) 
dfplus <- data.frame()
for (x in 1:length(tempr)){ 
  oneday <- BEM(larvawgt=w_ini, larvamm=l_ini, temp=tempr[x], P=0.87)
  dfplus <- rbind(dfplus, oneday)
}
l_ini = 12 #mm
w_ini = 0.0055 * ((l_ini/10)^3.19) 
dfmin <- data.frame()
for (x in 1:length(tempr)){ 
  oneday <- BEM(larvawgt=w_ini, larvamm=l_ini, temp=tempr[x], P=0.58)
  dfmin <- rbind(dfmin, oneday)
}
df_reg2 <- data.frame(Temp = tempr,
                      GR = df$GR,
                      GRmax = dfplus$GR,
                      GRmin = dfmin$GR)
rm(dfplus,dfmin)

p1 <- ggplot()+
  geom_line(data = df_long, aes(x=Temp, y=value, colour = variable), size = 1.2) +  
  scale_colour_manual(values=c(C="#663399",R="#339999",GR="#CC0033",Lab="#FF9933"))+
  scale_x_continuous(breaks=c(seq(-2,10,2))) +
  labs(y = expression(paste("Metabolic rates (g g"^"-1", "d"^"-1", ")")), x = "Temperature (ºC)")+
  theme_bw() + theme(legend.position = "none") + 
  geom_ribbon(data = df_reg2, aes(x=Temp, ymin = GRmin, ymax = GRmax), fill = "#CC0033", alpha = 0.2)+
  geom_text(data = df_long, aes(x=9.7, y=0.08), label="C", cex = 4, col="#663399") +
  geom_text(data = df_long, aes(x=9.7, y=0.058), label="R", cex = 4, col="#339999") +
  geom_text(data = df_long, aes(x=9.7, y=0.041), label="Lab", cex = 4, col="#FF9933") +
  geom_text(data = df_long, aes(x=9.7, y=0.028), label="GR", cex = 4, col="#CC0033") +
  annotate("text",  x=Inf, y = Inf, label = "A", cex = 5, vjust=1, hjust=1)

# analysis on minP in EstimateParameter.R line 350
p2 <- ggplot(data = minP12, aes(x = T, y = P)) +
  labs(x = "Temperature (ºC)", y = "Min P")+
  scale_y_continuous(breaks=c(seq(0,0.6,0.05))) +
  scale_x_continuous(breaks=c(seq(-2,10,2))) +
  geom_line(size = 1)+
  geom_ribbon(data = minP12, aes(x=T, ymin = minP8$P, ymax = minP27$P), fill = "mediumpurple4", alpha = 0.2)+
  theme_bw() + theme(legend.position="none") +
  geom_text(x=7.5, y=0.15, label="min P = 0.1244 + 0.0324*T \n (T > 0ºC)", color = "black", cex = 3) +
  annotate("text",  x=Inf, y = Inf, label = "B", cex = 5, vjust=1, hjust=1) 

ggarrange(p1,p2, nrow = 2, legend = F)

dev.copy2pdf(file="Figure7.pdf", width=5.5, height=8)

# =============================== #


# ========================== #
# Figure 8 ####
# ========================== #
# create dfs for model, field (Bouchard), labs (Laurel, Koenker, Hop, Kunz)
# plot length-at-age for all dfs

# modified BEM with GR from lab estimates
source("BEM_testing/BioenergeticInvers.R")

# build a function
# Growth of one larva at a given temperature = temp, with an initial length at hatch = larvamm_ini#
larvaGrowth <- function(larvamm_ini, days, temp, P){
  larvamm = larvamm_ini
  larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
  df <- data.frame()
  for (x in 1:days){ 
    oneday <- BEM(larvawgt, larvamm, temp, P, CDn = 1715.8)
    larvamm = oneday$larvamm
    larvawgt = oneday$larvawgt
    df <- rbind(df, oneday)
  } 
  return(df)
}

# apply growth function with different food$T scenarios
df6Cal <- larvaGrowth(larvamm_ini=6, days=30, temp=0, P=0.6)
df3Cal <- larvaGrowth(larvamm_ini=6, days=30, temp=0, P=0.6)


# food levels and type
# Artemia cal = 1715.8*(20.9/27.7)
larvaGrowth <- function(larvamm_ini, days, temp, P){
  larvamm = larvamm_ini
  larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
  df <- data.frame()
  for (x in 1:days){ 
    oneday <- BEM(larvawgt, larvamm, temp, P, CDn = 1294.6)
    larvamm = oneday$larvamm
    larvawgt = oneday$larvawgt
    df <- rbind(df, oneday)
  } 
  return(df)
}
df6Art <- larvaGrowth(larvamm_ini=6, days=30, temp=0, P=0.6)
df3Art <- larvaGrowth(larvamm_ini=6, days=30, temp=0, P=0.3)


# growth over n days #
age <- seq(1, 30, 1) # in grams
larvamm_ini <- 6

# field growth rate
larvamm = larvamm_ini #mm
larvawgt =  0.0055 * (larvamm/10)^3.19
df3 <- data.frame()
for (x in 1:length(age)){ 
  larvamm = larvamm + 0.2 # mm/day growth increment Bouchard and Fortier 2011
  larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
  oneday <- data.frame(larvawgt, larvamm)
  df3 <- rbind(df3, oneday)
} 

# Labs 
# Lab1GR
larvamm = larvamm_ini #mm
larvawgt =  0.0055 * (larvamm/10)^3.19  # W, grams; L, mm
df4 <- data.frame()
for (x in 1:length(age)){ 
  oneday <- BEMLab(larvawgt, larvamm, temp=0, P=0.6)
  larvawgt = oneday$larvawgt
  larvamm = oneday$larvamm
  df4 <- rbind(df4, oneday)
}
# Lab2GR
larvamm = larvamm_ini #mm
larvawgt =  0.0055 * (larvamm/10)^3.19
df5 <- data.frame()
for (x in 1:length(age)){ 
  larvamm = larvamm + 0.068 # mm/day growth increment Graham and Hop 1995
  larvawgt = 0.0055 * (larvamm/10)^3.19 #grams
  oneday <- data.frame(larvawgt, larvamm)
  df5 <- rbind(df5, oneday)
} 

# need to create df_long2
age30=30
df_reg <- data.frame(Age = age30,
                     Laptev2003 = 2.840209 + 0.207227*age30,
                     Laptev2005 = 5.853089 + 0.182449*age30,     # validation
                     Laptev2007 = 8.510625 + 0.189755*age30,
                     Beaufort2005 = 7.77515 + 0.19764*age30,     # validation
                     #Beaufort2006 = 13.45352 + 0.16090*age30,
                     Beaufort2010 = 5.629976 + 0.171524*age30,   
                     Beaufort2011 = 2.469326 + 0.198197*age30,   # validation
                     Beaufort2013 = 4.884622 + 0.185294*age30,
                     Beaufort2014 = 8.92386 + 0.17543*age30,     # validation
                     Kitikmeot2005 = 4.202389 + 0.205647*age30,  # validation
                     Kitikmeot2006 = 5.04346 + 0.20379*age30,
                     Kitikmeot2011 = 0.955450 + 0.231559*age30,  # validation
                     Kitikmeot2014 = 7.21992 + 0.16075*age30,    # validation
                     Kitikmeot2015 = 2.788851 + 0.199910*age30,
                     Baffin2005 = 3.23820 + 0.20863*age30,       # validation
                     Baffin2006 = 6.55684 + 0.19465*age30,
                     Baffin2007 = 7.4770 + 0.19465*age30,
                     Baffin2008 = 1.776091 + 0.209800*age30)

df_long2 <- reshape2::melt(df_reg, id = "Age")
df_long2$variable <- "Field30"

# for larva
df_reg <- data.frame(Age = age,
                     BEM = df6Cal$larvamm,
                     FieldGR = df3$larvamm,
                     Lab1GR = df4$larvamm,
                     Lab2GR = df5$larvamm)
df_long <- reshape2::melt(df_reg, id = "Age")
p1 <- ggplot(df_long, aes(y=value, x=Age, color = variable)) + geom_line() + 
  geom_boxplot(data=df_long2, alpha=0) + 
  scale_x_continuous(breaks=c(seq(0,30,5))) +
  scale_y_continuous(breaks=c(seq(6,14,1)), limits = c(6,14)) + 
  labs(x = "Age (days)", y ="Length (mm)", linetype="Models")+ 
  geom_text(x=2.5, y=13.5, label="A", color = "black", cex = 5) +
  geom_text(x=7, y=12, label=expression(paste("P = 0.6, ", italic("Calanus"))), color = "black") + 
  theme_bw() +
  theme(legend.title=element_blank(), axis.title.x = element_blank())


df_reg <- data.frame(Age = age,
                     BEM = df6Art$larvamm,
                     FieldGR = df3$larvamm,
                     Lab1GR = df4$larvamm,
                     Lab2GR = df5$larvamm)
df_long <- reshape2::melt(df_reg, id = "Age")
p2 <- ggplot(df_long, aes(y=value, x=Age, color = variable)) + geom_line() + 
  geom_boxplot(data=df_long2, alpha=0) + 
  scale_x_continuous(breaks=c(seq(0,30,5))) +
  scale_y_continuous(breaks=c(seq(6,14,1)), limits = c(6,14)) +
  labs(x = "Age (days)", y ="Length (mm)", linetype="Models")+ 
  geom_text(x=2.5, y=13.5, label="B", color = "black", cex = 5) +
  geom_text(x=7, y=12, label=expression(paste("P = 0.6, ", italic("Artemia"))), color = "black") + theme_bw() +
  theme(legend.title=element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())


df_reg <- data.frame(Age = age,
                     BEM = df3Cal$larvamm,
                     FieldGR = df3$larvamm,
                     Lab1GR = df4$larvamm,
                     Lab2GR = df5$larvamm)
df_long <- reshape2::melt(df_reg, id = "Age")
p3 <- ggplot(df_long, aes(y=value, x=Age, color = variable)) + geom_line() + 
  geom_boxplot(data=df_long2, alpha=0) + 
  scale_x_continuous(breaks=c(seq(0,30,5))) +
  scale_y_continuous(breaks=c(seq(6,14,1)), limits = c(6,14)) +
  labs(x = "Age (days)", y ="Length (mm)", linetype="Models")+ 
  geom_text(x=2.5, y=13.5, label="C", color = "black", cex = 5) +
  geom_text(x=7, y=12, label=expression(paste("P = 0.3, ", italic("Calanus"))), color = "black") + theme_bw()


df_reg <- data.frame(Age = age,
                     BEM = df3Art$larvamm,
                     FieldGR = df3$larvamm,
                     Lab1GR = df4$larvamm,
                     Lab2GR = df5$larvamm)
df_long <- reshape2::melt(df_reg, id = "Age")
p4 <- ggplot(df_long, aes(y=value, x=Age, color = variable)) + geom_line() + 
  geom_boxplot(data=df_long2, alpha=0) + 
  scale_x_continuous(breaks=c(seq(0,30,5))) +
  scale_y_continuous(breaks=c(seq(6,14,1)), limits = c(6,14)) +
  labs(x = "Age (days)", y ="Length (mm)", linetype="Models")+ 
  geom_text(x=2.5, y=13.5, label="D", color = "black", cex = 5) +
  geom_text(x=7, y=12, label=expression(paste("P = 0.3, ", italic("Artemia"))), color = "black") + theme_bw() +
  theme(legend.title=element_blank(), axis.title.y = element_blank())


ggpubr::ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2, common.legend = TRUE)

dev.copy2pdf(file="Figure8.pdf", width=8, height=5)


# ======================== #
# SUPPLEMENTAL MATERIAL ####
# ======================== #
# 
source("extractSST.R")
source("extractDates.R")
require(ncdf4)
#
# FIGURES S1-S4 ####
# SST Maps Regions-Years 

# Figure S1 Baffin ####
#
file=nc_open("SST_downloads/Baffin.nc") # 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file, sst)
sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_map[lon<290, lat<72,] <- NA #for Baffin
sst_mean2005 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2006-04-15", date_end = "2006-09-15")
sst_map[lon<290, lat<72,] <- NA #for Baffin
sst_mean2006 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2007-04-15", date_end = "2007-09-15")
sst_map[lon<290, lat<72,] <- NA #for Baffin
sst_mean2007 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2008-04-15", date_end = "2008-09-15")
sst_map[lon<290, lat<72,] <- NA #for Baffin
sst_mean2008 <- apply( sst_map , 1:2 , mean )

new.par=par(mfrow=c(2,2), oma = c(0.2, 0.5, 0.5, 0.1), # two rows of text at the outer left and bottom margin
            mar = c(2, 2.5, 1, 1.7), # space for one row of text at ticks and to separate plots
            mgp = c(1.5, 0.5, 0))    # axis label at 2 rows distance, tick labels at 1 row
plot3D::image2D(z=sst_mean2005, x = 360-lon, y = lat, main="2005", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2),ylim = c(70,79))
plot3D::image2D(z=sst_mean2006, x = 360-lon, y = lat, main="2006", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2),ylim = c(70,79))
plot3D::image2D(z=sst_mean2008, x = 360-lon, y = lat, main="2008", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2),ylim = c(70,79))

dev.copy2pdf(file="FigureS1.pdf", width=10, height=8)


# Fig S2 Beaufort ####
# Beaufort
file=nc_open("SST_downloads/Beaufort.nc") # 
sst <- ncvar_get(file, varid = "sst") 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file, sst)
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_mean2005 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2010-04-15", date_end = "2010-09-15")
sst_mean2010 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2011-04-15", date_end = "2011-09-15")
sst_mean2011 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2014-04-15", date_end = "2014-09-15")
sst_mean2014 <- apply( sst_map , 1:2 , mean )

sst_Beauf <- data.frame(Year = c(2005, 2010, 2011, 2014), 
                        meanSST = c(mean(sst_mean2005, na.rm = T),mean(sst_mean2010, na.rm = T),
                                    mean(sst_mean2011, na.rm = T), mean(sst_mean2014, na.rm = T)),
                        sdSST = c(sd(sst_mean2005, na.rm = T), sd(sst_mean2010, na.rm = T),
                                  sd(sst_mean2011, na.rm = T), sd(sst_mean2014, na.rm = T)))


new.par=par(mfrow=c(2,2), oma = c(0.2, 0.5, 0.5, 0.1), # two rows of text at the outer left and bottom margin
            mar = c(2, 2.5, 1, 1.7), # space for one row of text at ticks and to separate plots
            mgp = c(1.5, 0.5, 0))    # axis label at 2 rows distance, tick labels at 1 row
plot3D::image2D(z=sst_mean2005, x = 360-lon, y = lat, main="2005", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1, 4))
plot3D::image2D(z=sst_mean2010, x = 360-lon, y = lat, main="2010", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1, 4))
plot3D::image2D(z=sst_mean2011, x = 360-lon, y = lat, main="2011", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1, 4))
plot3D::image2D(z=sst_mean2014, x = 360-lon, y = lat, main="2014", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1, 4))

dev.copy2pdf(file="FigureS2.pdf", width=10, height=8)


# Fig S3 Kitik ####
# Kitikmeot #
file=nc_open("SST_downloads/Kitikmeot.nc") # 
sst <- ncvar_get(file, varid = "sst") 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file, sst)
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_mean2005 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2006-04-15", date_end = "2006-09-15")
sst_mean2006 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2011-04-15", date_end = "2011-09-15")
sst_mean2011 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2015-04-15", date_end = "2015-09-15")
sst_mean2015 <- apply( sst_map , 1:2 , mean )


new.par=par(mfrow=c(2,2), oma = c(0.2, 0.5, 0.5, 0.1), # two rows of text at the outer left and bottom margin
            mar = c(2, 2.5, 1, 1.7), # space for one row of text at ticks and to separate plots
            mgp = c(1.5, 0.5, 0))    # axis label at 2 rows distance, tick labels at 1 row
plot3D::image2D(z=sst_mean2005, x = 360-lon, y = lat, main="2005", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2))
plot3D::image2D(z=sst_mean2006, x = 360-lon, y = lat, main="2006",xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2))
plot3D::image2D(z=sst_mean2011, x = 360-lon, y = lat, main="2011", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2))
plot3D::image2D(z=sst_mean2015, x = 360-lon, y = lat, main="2015", xlab="Longitude Wº", ylab="Latitude Nº", zlim=c(-1.8, 2))

dev.copy2pdf(file="FigureS3.pdf", width=10, height=8)




# Figure S4 Laptev ####

file=nc_open("SST_downloads/Laptev.nc") # 
sst <- ncvar_get(file, varid = "sst") 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file, sst)
sst_map <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2003-04-15", date_end = "2003-09-15")
sst_map[, lat>81,] <- NA 
sst_mean2003 <- apply( sst_map , 1:2 , mean )

sst_map <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_map[, lat>81,] <- NA 
sst_mean2005 <- apply( sst_map , 1:2 , mean )

sst_map <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2007-04-15", date_end = "2007-09-15")
sst_map[, lat>81,] <- NA 
sst_mean2007 <- apply( sst_map , 1:2 , mean )

# Laptev, period: Febr 01-Sep 15
new.par=par(mfrow=c(2,2), oma = c(0.2, 0.5, 0.5, 0.1), # two rows of text at the outer left and bottom margin
            mar = c(2, 2.5, 1, 1.7), # space for one row of text at ticks and to separate plots
            mgp = c(1.5, 0.5, 0))    # axis label at 2 rows distance, tick labels at 1 row
plot3D::image2D(z=sst_mean2003, x = lon, y = lat, main="2003", xlab="Longitude Eº", ylab="Latitude Nº", zlim=c(-1.8, 2),ylim = c(70,81))
plot3D::image2D(z=sst_mean2005, x = lon, y = lat, main="2005", xlab="Longitude Eº", ylab="Latitude Nº", zlim=c(-1.8, 2),ylim = c(70,81))
plot3D::image2D(z=sst_mean2007, x = lon, y = lat, main="2007", xlab="Longitude Eº", ylab="Latitude Nº", zlim=c(-1.8, 2),ylim = c(70,81))

dev.copy2pdf(file="FigureS4.pdf", width=10, height=8)




# FIGURES S5:S14 on fitting parameter P versus model validation
# Line 210 in EstimateParam_Validation.R

