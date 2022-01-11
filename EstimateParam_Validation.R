# =========================== #
# Estimate Parameter P        #
# =========================== #
# revised December 2021, Carmen L David
#
# Includes:
#   - Get daily SST
#   - Split observation data for fitting vs validation
#   - Fit parameter P, model validation, estimate model efficiency (EF)
#   - Model validation with daily SST for 2005 data
#   - estimate minimum value of P ~ temperature
#   - Supplemental Material: Figures S5-S14, Table S2
#
#
#
setwd("/../")
source("BioenergeticIndModel.R")
source("extractSST.R")
source("extractDates.R")

require(ncdf4)
require(ggplot2) 
require(ggpubr)
require(reshape2)


# ====================== #
# daily SST           ####
# ====================== #
# 
# get SST_all Regions ####
# Regions selection after Figure 1, SST data downloaded from:
#
# Menemenlis, D, Campin, J, Heimbach, P, Hill, C, Lee, T, Nguyen, A, Schodlok, M, Zhang, H. 2008. 
# ECCO2: High resolution global ocean and sea ice data synthesis. Mercator Ocean Quarterly Newsletter 31: 13–21. 
# Available at http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_daa1_6006_642b.html. Accessed 2021 May 7. 
#
Region_names <- c("Baffin", "Beaufort", "Kitikmeot", "Laptev")
#
#Baffin #
#
l=1 # index for region labels
# read SST file
file=nc_open("SST_downloads/Baffin.nc") # 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file)

sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_map[lon<290, lat<72,] <- NA # cut off cells not included
sst_mean2005 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2006-04-15", date_end = "2006-09-15")
sst_map[lon<290, lat<72,] <- NA # cut off cells not included
sst_mean2006 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Baffin.nc", date_start = "2008-04-15", date_end = "2008-09-15")
sst_map[lon<290, lat<72,] <- NA # cut off cells not included
sst_mean2008 <- apply( sst_map , 1:2 , mean )

sst_all <- data.frame(Region = Region_names[l],
                      Year = c(2005, 2006, 2008), 
                      meanSST = c(mean(sst_mean2005, na.rm = T),mean(sst_mean2006, na.rm = T),mean(sst_mean2008, na.rm = T)),
                      sdSST = c(sd(sst_mean2005, na.rm = T), sd(sst_mean2006, na.rm = T),sd(sst_mean2008, na.rm = T)))
rm(list=ls(pattern="sst_mean"))

#
l=2 # index for region labels
# get SST file
# Beaufort Sea #
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_mean2005 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2010-02-15", date_end = "2010-09-15")
sst_mean2010 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2011-04-15", date_end = "2011-09-15")
sst_mean2011 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Beaufort.nc", date_start = "2014-04-15", date_end = "2014-09-15")
sst_mean2014 <- apply( sst_map , 1:2 , mean )

sst_all <- rbind(sst_all, 
                 data.frame(Region = Region_names[l],
                            Year = c(2005, 2010, 2011, 2014), 
                            meanSST = c(mean(sst_mean2005, na.rm = T),mean(sst_mean2010, na.rm = T),
                                     mean(sst_mean2011, na.rm = T), mean(sst_mean2014, na.rm = T)),
                            sdSST = c(sd(sst_mean2005, na.rm = T), sd(sst_mean2010, na.rm = T),
                                      sd(sst_mean2011, na.rm = T), sd(sst_mean2014, na.rm = T))))
rm(list=ls(pattern="sst_mean"))

#
l=3 # index for region labels
#
file_nc = paste("SST_downloads/Kitikmeot.nc", sep="")
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_mean2005 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2006-04-15", date_end = "2006-09-15")
sst_mean2006 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2011-04-15", date_end = "2011-09-15")
sst_mean2011 <- apply( sst_map , 1:2 , mean )
sst_map <- extractDailySST(file_sst = "SST_downloads/Kitikmeot.nc", date_start = "2015-04-15", date_end = "2015-09-15")
sst_mean2015 <- apply( sst_map , 1:2 , mean )

sst_all <- rbind(sst_all, 
                 data.frame(Region = Region_names[l],
                            Year = c(2005, 2006, 2011, 2015), 
                            meanSST = c(mean(sst_mean2005, na.rm = T),mean(sst_mean2006, na.rm = T),
                                        mean(sst_mean2011, na.rm = T), mean(sst_mean2015, na.rm = T)),
                            sdSST = c(sd(sst_mean2005, na.rm = T), sd(sst_mean2006, na.rm = T),
                                      sd(sst_mean2011, na.rm = T), sd(sst_mean2015, na.rm = T))))
rm(list=ls(pattern="sst_mean"))

# Laptev #
#
l=4 # index for region labels
#
file=nc_open("SST_downloads/Laptev.nc") # 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file)
sst_map <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2003-04-15", date_end = "2003-09-15")
sst_map[, lat>81,] <- NA 
sst_mean2003 <- apply( sst_map , 1:2 , mean )

sst_map <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2005-04-15", date_end = "2005-09-15")
sst_map[, lat>81,] <- NA 
sst_mean2005 <- apply( sst_map , 1:2 , mean )

sst_map <- extractDailySST(file_sst = "SST_downloads/Laptev.nc", date_start = "2007-04-15", date_end = "2007-09-15")
sst_map[, lat>81,] <- NA 
sst_mean2007 <- apply( sst_map , 1:2 , mean )

sst_all <- rbind(sst_all, 
                 data.frame(Region = Region_names[l],
                            Year = c(2003, 2005, 2007), 
                            meanSST = c(mean(sst_mean2003, na.rm = T),mean(sst_mean2005, na.rm = T),
                                        mean(sst_mean2007, na.rm = T)),
                            sdSST = c(sd(sst_mean2003, na.rm = T), sd(sst_mean2005, na.rm = T),
                                      sd(sst_mean2007, na.rm = T))))
rm(list=ls(pattern="sst_mean"))
rm(lon, lat, sst_map)
#
# Get SST_all end #
# =============== #



# ====================== #
# Fit / validation    ####
# ====================== #
#
df_export <- data.frame()
# repeat per region (l = 1:4)
for (l in 1:4){ 
  # extract mean SST in all years
  Year_all <- sst_all[sst_all$Region == Region_names[l],] $Year
  Temp_all <- sst_all[sst_all$Region == Region_names[l],] $meanSST
  #
  # Read data on length-at-age from otolith measurements (see Bouchard and Fortier 2011, 2020)
  filename_obs <- paste("Datasets/DatabaseLengthAge_", Region_names[l], ".txt", sep="")
  df_obs_all <- read.table(filename_obs, header=T, sep="\t", fill = TRUE)
  df_obs_all$Year <- as.factor(df_obs_all$Year)
  df_obs_all <- df_obs_all[!is.na(df_obs_all$Age),] #delete rows with Age=NA

  # repeat per year (k)
  for (k in 1:length(Year_all)){ 
    # select year
    df_obs <- df_obs_all[df_obs_all$Year==Year_all[k],]
    
    # split dataset for parameter fitting or validation:
    # df_fit (60%), df_val (40%)
    sample_size = floor(0.6*nrow(df_obs))
    set.seed(777)

    picked = sample(seq_len(nrow(df_obs)),size = sample_size)
    df_fit =df_obs[picked,]
    df_val =df_obs[-picked,]

    # ============ #
    # FITTING P ####
    # ============ #
    age <- seq(min(df_obs$Age),max(df_obs$Age),1)

    #define intervals pf P to be tested
    P_seq <- seq(0.5, 0.85, 0.01)
    ef_fit <- rep(0,length(P_seq))

    for (j in 1:length((P_seq))){
      # start with min length at hatch in observations data
      larvamm = min(df_obs$LS) # mm
      larvawgt = 0.0055 * ((larvamm/10)^3.19) # grams
      df = data.frame() # age = 0, W = larvawgt, L = larvamm)
      for (i in 1:length(age)){
        # daily growth with bioenergetic model (BEM)
        oneDay <- BEM(larvawgt, larvamm, temp=Temp_all[k], P=P_seq[j])
        larvawgt = oneDay$larvawgt
        larvamm = oneDay$larvamm
        df = rbind(df, data.frame(age = age[i], W = larvawgt, L = larvamm))
      }

      lm2 <- lm(L~age, data = df) 
      # df_fit$Age : observed age in the subset used to fit parametr P
      pred.lm2 <- unlist(lapply(1:length(df_fit$Age), function(x) {df[which(df$age == df_fit$Age[x]),]$L}))

      # df_fit$LS  # observed values in the subset used to fit parametr P
      # EF: MODEL EFFICIENCY calculated using different values of P
      ef_fit[j] <- 1 - sum((df_fit$LS-pred.lm2)^2) / sum((df_fit$LS-mean(df_fit$LS))^2)
    }

    # select max model efficiency for a certain value of OP
    whichP <- which(ef_fit==max(ef_fit))
    P_seq[whichP]


    # ============= #
    # VALIDATION ####
    # ============= #
    # same routine as previously fitting parameter
    age <- seq(min(df_obs$Age),max(df_obs$Age),1)
    larvamm = min(df_obs$LS) # mm
    larvawgt = 0.0055 * ((larvamm/10)^3.19) # grams
    df = data.frame() # age = 0, W = larvawgt, L = larvamm)
    for (i in 1:length(age)){
      # daily growth with bioenergetic model (BEM)
      oneDay <- BEM(larvawgt, larvamm, temp=Temp_all[k], P=P_seq[whichP])
      larvawgt = oneDay$larvawgt
      larvamm = oneDay$larvamm
      df = rbind(df, data.frame(age = age[i], W = larvawgt, L = larvamm))
    }

    lm1 <- lm(LS~Age, data = df_fit) # regression on fitted subset data
    lm2 <- lm(L~age, data = df)  # regression on validation subset data

    # SUPPLEMENTAL MATERIAL
    # FIGURES S5:S14 on fitting parameter P versus model validation 
    title = paste(Region_all[l], Year_all[k], sep = " ")
    plot(LS~Age, data = df_fit, type='p', pch=20, main = title, ylab = "Length (mm)", xlab = "Age (days)")
    points(LS~Age, data = df_val, col=alpha('red', 1), bg = alpha('red', 0.4), pch=22, cex = 0.8)
    abline(lm1, lwd=2)
    lines(df$L~df$age, col='red', lwd=2)
    filename=paste("plots_fitP/", Region_all[l], Year_all[k],".pdf", sep = "")
    dev.copy2pdf(file=filename, width=6, height=5)

    # df_val$Age # observed age, used to fit
    # df_val$LS  # observed values
    pred.LS <- unlist(lapply(1:length(df_val$Age), function(x) {df[which(df$age == df_val$Age[x]),]$L}))
    # EF: MODEL EFFICIENCY
    ef <- 1 - sum((df_val$LS-pred.LS)^2) / sum((df_val$LS-mean(df_val$LS))^2)
    
    # SUPPLEMENTAL MATERIAL
    # Table S2 on model efficiency
    df_export <- rbind(df_export, data.frame(Region = Region_all[l], 
                                          Year = Year_all[k], 
                                          n_fit = nrow(df_fit),
                                          P = P_seq[whichP], 
                                          EF = round(ef, digits = 3),
                                          n_val = nrow(df_val)))

  } # close loop for year
} # close loop for region

# FOR SUPPLEMENTAL MATERIAL
# Table S2 on model efficiency
write.table(df_export, file = "FittingP_ModelEfficiency60-40.txt", sep = "\t")



# ================== #
# Validation 2005 ####
# ================== #
#
# Baffin 
# get SST_all Regions
#
Region_all <- c("Baffin", "Beaufort", "Kitikmeot", "Laptev")
df_export <- data.frame()
l_ini_all = c(4,5,6,7)
P_val_all = c(0.85, 0.61, 0.81, 0.71)

l=4
# get SST file
filename <- paste("SST_downloads/", Region_all[l], ".nc", sep="")
file=nc_open(filename) # 
sst <- ncvar_get(file, varid = "sst") 
lat <- ncvar_get(file, varid = "latitude")
lon <- ncvar_get(file, varid = "longitude")
rm(file, sst)

# Baffin #
sst_map <- extractDailySST(file_sst = filename, date_start = "2005-02-15", date_end = "2005-09-15")
sst_Baf <- apply( sst_map , 3 , mean, na.rm = T)
rm(lat, lon)


# read observation data
filename_obs <- paste("Datasets/DatabaseLengthAge_", Region_all[l], ".txt", sep="")
df_obs <- read.table(filename_obs, header=T, sep="\t", fill = TRUE)
df_obs$Year <- as.factor(df_obs$Year)
df_obs <- df_obs[!is.na(df_obs$Age),] #delete rows with Age=NA
df_obs <- df_obs[df_obs$Year==2005 & df_obs$Age < 213,c("LS","Age")]

# estimate L

for (k in 1:length(l_ini_all)){ 

  ndays = length(sst_Baf)
  l_ini = l_ini_all[k]
  larvawgt = 0.0055 * ((larvamm/10)^3.19) # mm
  P_val = P_val_all[l]

  df_hatch <- data.frame()
  for (i in 1: ndays){ 
      larvamm = l_ini
      larvawgt = 0.0055 * ((larvamm/10)^3.19)
      for (j in i:ndays){
         oneDay <- BEM(larvawgt, larvamm, temp=sst_Baf[j], P=P_val)
         larvawgt = oneDay$larvawgt
         larvamm = oneDay$larvamm
      }
      df_hatch = rbind(df_hatch, data.frame(Age = ndays - i, L = larvamm))
  }

lm1 <- lm(LS~Age, data = df_obs)
#lm2 <- lm(L~Age, data = df_hatch)  
#age <- seq(min(df_obs$Age),max(df_obs$Age),1)

title = paste(Region_all[l], "2005", sep = " ")
plot(LS~Age, data = df_obs, type='p', pch=20, main = title, ylab = "Length (mm)", xlab = "Age (days)")
abline(lm1, lwd=2)
lines(df_hatch$L~df_hatch$Age, col='red', lwd=2)
filename=paste("plots_valid2005/", Region_all[l], l_ini,".pdf", sep = "")
dev.copy2pdf(file=filename, width=6, height=5)


df_obs$pred <- unlist(lapply(1:length(df_obs$Age), function(x) {df_hatch[which(df_hatch$Age == df_obs$Age[x]),]$L}))

ef <- 1 - sum((df_obs$LS-pred.LS)^2) / sum((df_obs$LS-mean(df_obs$LS))^2)

# Select APR-MAY hatchers, EF NOT IMPROVED!!
# df_obs <- df_obs[df_obs$Age>108 & df_obs$Age<168,]
# ef <- 1 - sum((df_obs$LS-df_obs$pred)^2) / sum((df_obs$LS-mean(df_obs$LS))^2)

df_export <- rbind(df_export, data.frame(Region = Region_all[l],
                                         n = nrow(df_obs),
                                         P = P_val, 
                                         L_ini = l_ini,
                                         EF = round(ef, digits = 3)))
}

# Baffin2005 L_ini = 3.23



# ============= #
# Min P f(T) ####
# ============= #
#
#  over a range of temperatures
tempr <- seq(-1.5, 10, 0.1)
# P_val takes values +/-20%
P_val = seq(0, 0.85, 0.001)
l_ini = 12 #mm Alternative calculate for SL = 8, 27 mm (min and max larval length)
w_ini = 0.0055 * ((l_ini/10)^3.19) 

minP <- data.frame()
for (x in 1:length(tempr)){ 
  df <- data.frame()
  for (l in 1:length(P_val)){ 
  Growth <- BEM(larvawgt=w_ini, larvamm=l_ini, temp=tempr[x], P=P_val[l])$GR
  oneday <- data.frame(Temp = tempr[x], P = P_val[l], GR = Growth)
  df <- rbind(df, oneday)
}
  GR_poz <- which(df$GR > 0)
  minP <- rbind(minP, data.frame(T = df[GR_poz[1],]$Temp,
                                P = P_val[GR_poz[1]]))
}


 minP12 <- minP
 # same for min and max larval length
 # minP8 <- minP
 # minP27 <- minP

lmP <- lm(P~T, data = minP[minP$T>=0,],)
summary(lmP)


lmP <- lm(P~T+I(T^2), data = minP)
summary(lmP)

newT <- seq(-1.8,8, 0.1)
minPredict <- predict(lmP,list(T=newT))

#create scatterplot of original data values
plot(minP$T, minP$P, pch=16)
#add predicted lines based on quadratic regression model
lines(newT, minPredict, col='blue')

# Figure 7B #
ggplot(data = minP12, aes(x = T, y = P)) +
  labs(x = "Temperature (ºC)", y = "Min P")+
  scale_y_continuous(breaks=c(seq(0,0.6,0.05))) +
  scale_x_continuous(breaks=c(seq(-2,10,2))) +
  geom_line(size = 1)+
  geom_ribbon(data = minP12, aes(x=T, ymin = minP8$P, ymax = minP27$P), fill = "mediumpurple4", alpha = 0.2)+
  theme_bw() + theme(legend.position="none") +
  geom_text(x=7.5, y=0.15, label="min P = 0.1244 + 0.0324*T \n (T > 0ºC)", color = "black", cex = 3) +
  annotate("text",  x=Inf, y = Inf, label = "B", cex = 5, vjust=1, hjust=1)

dev.copy2pdf(file="Figure.pdf", width=5.5, height=4)

