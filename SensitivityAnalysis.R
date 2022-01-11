# ==================== #
# SENSITIVITY ANALYSES #
# ==================== #
# revised December 2021, Carmen L David
# 
# Includes analyses:
#    - Individual parameter perturbation
#    - Global analysis, generate random parameters matrix with Latin hypercurbe
#
# ============================= #
# Individual Param Perturbation ####
# ============================= #
#
setwd("/../")
source("BioenergeticIndModel.R")
#
# create reference using baseline values for model parameters
ndays = 30
nlarva = 1
temp=0  #(reference temperature)
# set initial larval length/weight
larvamm = 6      # in mm  
larvawgt <- 0.0055 * ((larvamm/10)^3.19)  # in grams
all_larva = data.frame(age = 0,L = larvamm)
for (i in 1:ndays){
  # daily growth
  oneDay <- BEM(larvawgt, larvamm, temp, P = 0.73)
  larvawgt <- oneDay$larvawgt 
  larvamm <- oneDay$larvamm
  all_larva = rbind(all_larva, data.frame(age = i, L = larvamm))
}
df <- all_larva
df$param <- "Ref"
df$test <- "Ref"  


# set variable values +/-20%
pms <- data.frame(temp=0, P=0.73, CDe = 475, CDn = 1715.8, CDf = 1327.3, Qox = 3234.5, SDA = 0.375)
pms[2,] <- 0.8*pms # -20%
pms[3,] <- 1.2*pms[1,] # +20%
pms$temp <- c(0, -1, 1)
pms[4,] <- as.numeric(1, 1, 1, 1, 2, 1, 2)

# individually select a parametr to vary

j=7 # column nb of parameter
larvamm = 6      # initial length in mm  
larvawgt <- 0.0055 * ((larvamm/10)^3.19)  # in grams
all_larva = data.frame(age = 0,L = larvamm)

for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp, SDA = pms[2,j]) # here change the argument to vary, corresponding to j in pms
  larvawgt <- oneDay$larvawgt 
  larvamm <- oneDay$larvamm
  all_larva = rbind(all_larva, data.frame(age = i, L = larvamm))
}
all_larva$param <- names(pms[j])
all_larva$test <- "-20%"
df <- rbind(df, all_larva)

# variate +20%
larvamm = 6      # in mm  
larvawgt <- 0.0055 * ((larvamm/10)^3.19)  # in grams
all_larva = data.frame(age = 0,L = larvamm)
for (i in 1:ndays){
  oneDay <- BEM(larvawgt, larvamm, temp, SDA = pms[3,j]) # here change the argument to vary, corresponding to j in pms
  larvawgt <- oneDay$larvawgt 
  larvamm <- oneDay$larvamm
  all_larva = rbind(all_larva, data.frame(age = i, L = larvamm))
}
all_larva$param <- names(pms[j])
all_larva$test <- "+20%"
df <- rbind(df, all_larva)


write.table(df, "Datasets/SensitivityAnalysis_Sep2021.txt")




# =========================== #
# Global Analsis T & Food  ####
# =========================== #
library(relaimpo)

ndays = 3
temp=seq(-1.8, 10, 0.2)
larvamm_ini = 12      # in mm  
larvawgt_ini <- 0.0055 * ((larvamm_ini/10)^3.19)  # in grams
df = data.frame( W = larvawgt_ini, L = larvamm_ini, GR = 0)

# parameters list for the analysis
n_par <- data.frame(T = 0,  P = 0.73, CDe = 475, CDn = 1715.8, CDf = 1327.3)


# Generate Parameters matrix with Latin hypercurbe
set.seed(1000)
A <- lhs::randomLHS(n=300, k = length(n_par))
B <- matrix(nrow = nrow(A), ncol = ncol(A))
colnames(B) <- names(n_par)
# create distribution for temperatures
B[,1] <- qunif(A[,1], min = min(-1.8), max = max(10))
# create distributions for food parameters
for (i in 2:length(n_par)){ 
  par_seq <- seq(0.8*n_par[1,i], 1.2*n_par[1,i],length=100)
  B[,i] <- rnorm(A[,i], mean = n_par[1,i], sd = sd(par_seq))
}               
B <- rbind(n_par,B)


# Apply growth model on parameters matrix
for (j in 1:nrow(B)) { 
  larvawgt <- larvawgt_ini# in grams
  larvamm =  ((larvawgt/0.0055)^(1/3.19)) * 10
  for (i in 1:ndays){
    one_larva <- BEM(larvawgt, larvamm, temp=B[j,1], P = B[j,2], CDe = B[j,3], 
                     CDn = B[j,4], CDf = B[j,5])
    larvawgt = one_larva$larvawgt
    larvamm = one_larva$larvamm
  }
  df_fin = data.frame(W = larvawgt, L = larvamm, GR = one_larva$GR)
  df <- rbind(df, df_fin)
}

# bind model output (length) to parameters matrix
# correlations
df <- cbind(df[2:nrow(df),], B)
df$S <- (df$L-df$L[1])/df$L[1]



library(ggplot2)
library(ggpubr)
expl_var <- "L"
expl_var_lab <- "Length (mm)"

p4 <- ggscatter(df, x = colnames(df[4]), y = expl_var,  add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",xlab = "Temperature (ºC)", ylab = expl_var_lab) +
  geom_smooth(method='loess', color="red") + ylim(12, 13.50) +
  annotate("text",  x=Inf, y = Inf, label = "A", cex = 5, vjust=1, hjust=1)
p5 <- ggscatter(df, x = colnames(df[5]), y = expl_var,  add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",xlab = colnames(df[5]), ylab = expl_var_lab) +
  geom_smooth(method='loess', color="red") + ylim(12, 13.50)+
  annotate("text",  x=Inf, y = Inf, label = "B", cex = 5, vjust=1, hjust=1)
p6 <- ggscatter(df, x = colnames(df[6]), y = expl_var,  add = "reg.line", conf.int = TRUE, 
                 cor.coef = TRUE, cor.method = "pearson",xlab = colnames(df[6]), ylab = expl_var_lab) +
  geom_smooth(method='loess', color="red") + ylim(12, 13.50)+
  annotate("text",  x=Inf, y = Inf, label = "C", cex = 5, vjust=1, hjust=1)
p7 <- ggscatter(df, x = colnames(df[7]), y = expl_var,  add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",xlab = colnames(df[7]), ylab = expl_var_lab) +
  geom_smooth(method='loess', color="red") + ylim(12, 13.50)+
  annotate("text",  x=Inf, y = Inf, label = "D", cex = 5, vjust=1, hjust=1)
p8 <- ggscatter(df, x = colnames(df[8]), y = expl_var,  add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",xlab = colnames(df[8]), ylab = expl_var_lab) +
  geom_smooth(method='loess', color="red") + ylim(12, 13.50)+
  annotate("text",  x=Inf, y = Inf, label = "E", cex = 5, vjust=1, hjust=1)

ggarrange(p4,p5,p6,p7,p8, ncol = 2, nrow = 3)
dev.copy2pdf(file="Figure8.pdf", width=10, height=9)


# ============================= #
# Estimate relative contribution
# ============================= #
lm1 <- lm(L~T+P+CDe+CDn+CDf, data = df)
summary(lm1)

# Calculate Relative Importance for Each Predictor
# Bootstrap Measures of Relative Importance (1000 samples)
boot <- boot.relimp(lm1, b = 1000, type = c("lmg"), rank = T,diff = T, rela = T)
#boot <- boot.relimp(lm1, b = 1000, type = c("lmg", "last", "first", "pratt"), rank = T,diff = T, rela = T)
booteval.relimp(boot) # print result
par(las = 1)
plot(booteval.relimp(boot,sort=T)) # plot result 
dev.copy2pdf(file="Figure8_inset.pdf", width=5, height=6)

boot.res <- data.frame(parameter = c("CDn", "P", "T", "CDf", "CDe"),
                       value = c(0.3870374,0.3374606,0.1676697,0.1057671,0.0020653)*100,
                       sd = c())
                       #sdplus = c(0.4645,0.4076,0.2297,0.1591, 0.0174),
                       #sdminus = c(0.3106,0.2675, 0.1031, 0.0583, 0.0004))

ggplot(data = boot.res, aes(x=reorder(parameter, -value), y=value), width=.5, ylim = c(0,50)) +
  geom_bar(stat="identity") +  
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2, position=position_dodge(.1))+
  theme_classic() +
  annotate("text",  x=Inf, y = Inf, label = "F", cex = 5, vjust=1, hjust=1)


