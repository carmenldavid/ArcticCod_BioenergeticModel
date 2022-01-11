# =============================== #
# Bioenergetic model 
# =============================== #
# after Thanassekos & Fortier 2012, following Roy et al. 2004
# last update December 2021, Carmen L David

# BASIC EQUATION
# growth (GR) = food consumption (C_real) * assimilation efficiency (A) - respiration (Resp) * (fact + SDA)
# temp limitation on food consumption from Thanassekos & Fortier 2012
# temp limitation for respiration from Grahan & Hop 1997 and Holeton 1978

# BASIC PARAMETERS
# larva length, in mm (larvamm)
# larva weight, in grams (larvawgt)
# water temperature (temp)

# =============================== #
BEM <- function (larvawgt, larvamm, temp, 
                 P = 0.73,        # proportion of food consumption
                 Lyex = 8.5,      # Length (mm) at yolk absorbtion 
                 CDe = 475,       # caloric density of Calanus eggs (cal/g)
                 CDn = 1715.8,    # caloric density of Calanus nauplii (cal/g)
                 CDf = 1327.3,    # caloric density of fish (cal/g)
                 ac = 0.0337,     # coefficient of weight-dependent regression of food consumption
                 bc = -0.2862,    # exponent of weight-dependent regression of food consumption
                 A = 0.8,         # assimilation efficiency
                 Qox = 3234.5,    # oxycaloric coefficient of Calanus prey (cal/g O2)
                 av = 0.0355,     # coefficient (intercept) of weight-dependent respiration equation
                 bv = -0.1699,    # exponent of weight-dependent respiration equation
                 SDA = 0.375)     # specific dynamic action
{
  # ================================== #
  # FOOD CONSUMPTION (C_real) #
  # ================================== #
  # fraction of nauplii in diet increases with larva length
  # equation estimated after Thanasselos & Fortier 2012
  # details in ParameterFitting.R
  npl = 1-exp(-0.19786*(larvamm-0.34))
  
 # Energy consumption when yolk reserves present, f(weight, Temp)
  C_yolk <- ac*(larvawgt^bc) * ifelse(temp<0.7, 0.943 + 0.155*temp - 0.103*temp^2, 1)
  # Maximum daily food consumption (C_max, g/g/day), f(weight, Food, Temp)
  C_max <- ac*(larvawgt^bc) * (npl*CDn + (1-npl)*CDe)/CDf * ifelse(temp<0.7, 0.943 + 0.155*temp - 0.103*temp^2, 1) 
  # Feeding type depends if larva length (larvamm in mm) > Lyex
  # Yolk present or exogenous feeding
  C_real = ifelse (larvamm < Lyex, C_yolk , P*C_max)
  
  # ================================== #
  # METABOLIC LOSSES (Resp) #
  # ================================== #
  # Respiration formulation after Thanassekos and Fortier 2012
  # temp-dependent equation from Hop-Graham1995
  # weight-dependent equation from Holeton 1974
  # grams_O2 converted into grams; metabolic rate becomes g/g/day
  # Fact = 1.5, increase in metabolism due to activity
  Resp <- av * (larvawgt^bv) * (0.0111 * temp + 0.042) * Qox/CDf * (1.5 + SDA)

  # ================================== #
  # GROWTH RATE (GR) #
  # ================================== #
  # GR is a ratio in g/g/day
  GR <- C_real * A - Resp 
  
  # larva is growing with daily increment in weight
  larvawgt <- larvawgt*(1+GR)
  # new daily length estimated from conversion equation of Geoffrey et al. 2016
  larvamm_new <- ((larvawgt/0.0055)^(1/3.19)) * 10
  # 
  GM <- larvamm_new-larvamm # daily growth increment in mm/day
  larvamm = max(larvamm_new, larvamm)

  ibm_param <- data.frame(larvawgt, larvamm, GR, GM, C_max, Resp)
  names(ibm_param) <- c("larvawgt", "larvamm", "GR", "GM", "Cmax", "R")
  return(ibm_param)
}

