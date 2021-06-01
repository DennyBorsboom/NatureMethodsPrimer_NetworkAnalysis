# Reproducible figure for Nature Methods primer paper, Borsboom et al. 2021.
# This examples contains a *subset* of variables collected and modeled in our covid19 paper
# This paper, with full data / paper / code can be found at: https://psyarxiv.com/36xkp
# Eiko Fried, March 14 2021



# -------------------------------------------------------------------------
# --------------- 1. Loading packages & Data ------------------------------
# -------------------------------------------------------------------------

library(ggplot2)
library(dplyr)
library(mlVAR)
library(qgraph)
library(bootnet)
library(reshape)
library(viridis)


# I use the below path structure to run the code
# To reproduce, adjust path structure to your system
setwd("~/Dropbox/Research/My Projects/2020 - Borsboom Fried et al, Nature Networks primer/analysis")

figs <- "./figures/" # figure directory
datapath <- "./data/" # data directory


# -------------------------------------------------------------------------
# --------------- 2. Load data --------------------------------------------
# -------------------------------------------------------------------------

# load data
load(paste0(datapath, "clean_network.RData"))
Data5b <- Data2



# -------------------------------------------------------------------------
# --------------- 3. Detrend data -----------------------------------------
# -------------------------------------------------------------------------

# some of the variables have trends (i.e. linear or other changes over time)
# we remove these trends before estimation network structures
# we work with all 16 variables here that we collected
# afterwards, we select those we want to estimate networks on

# Alpha to detrend:
alpha <- 0.05

# Variables to investigate:
vars <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11","Q12","Q13","Q14","Q15","Q16","Q17","Q18")

# Labels:
varLabs <- c("Relax","Irritable","Worry","Nervous","Future","Anhedonia",
             "Tired","Hungry","Alone","Angry","Social_offline","Social_online","Music",
             "Procrastinate","Outdoors","C19_occupied","C19_worry","Home")

names(Data5b)[names(Data5b) %in% vars] <- varLabs

# Remove items:
Data5b <- Data5b %>% select(-Hungry,-Angry,-Music,-Procrastinate,-Tired,-Outdoors,-Home,-C19_occupied,-C19_worry)
varLabs <- varLabs[!varLabs %in% c("Hungry","Angry","Music","Procrastinate","Tired","Outdoors","Home","C19_occupied","C19_worry")]

# Data frame with empty values for fitted effects (all):
fitted_all <- expand.grid(
  beep = seq(min(Data5b$beep),max(Data5b$beep)),
  day = seq(min(Data5b$day),max(Data5b$day))
)

# Data frame with empty values for day trends:
fitted_day <- data.frame(
  day = seq(min(Data5b$day),max(Data5b$day))
)

# Data frame with empty values for beeps:
fitted_beep <- data.frame(
  beep = seq(min(Data5b$beep),max(Data5b$beep))
)

# Data frame to store p-values:
p_values <- data.frame(
  var = c("day", "beep")
)

# Also empty data frame list for test statistics:
testStatistics <- list()
coefficients <- list()
stdcoefficients <- list()

# Make the beep variable factor in dataset:
Data5b$beepFactor <- factor(Data5b$beep, levels = 0:3, labels = c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00"))
fitted_all$beepFactor <- factor(fitted_all$beep, levels = 0:3, labels = c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00"))
fitted_beep$beepFactor <- factor(fitted_beep$beep, levels = 0:3, labels = c("09:00 - 12:00","12:00 - 15:00","15:00 - 18:00","18:00 - 21:00"))

# Make day variable for dates:
Data5b$date <- as.Date("2020-03-15") + Data5b$day
fitted_all$date <- as.Date("2020-03-15") + fitted_all$day
fitted_day$date <- as.Date("2020-03-15") + fitted_day$day

# Add the midpoints as time variable:
Data5b$midTime <- as.character(factor(Data5b$beep, levels = 0:3, labels = c("10:30","13:30","16:30","19:30")))
Data5b$midTime <- as.POSIXct(paste(Data5b$date,Data5b$midTime), format = "%Y-%m-%d %H:%M", tz = "Europe/Amsterdam")

fitted_all$midTime <- as.character(factor(fitted_all$beep, levels = 0:3, labels = c("10:30","13:30","16:30","19:30")))
fitted_all$midTime <- as.POSIXct(paste(fitted_all$date,fitted_all$midTime), format = "%Y-%m-%d %H:%M", tz = "Europe/Amsterdam")

# Data frame to store detrended data:
data_detrended <- Data5b

# Fix curves:
for (v in seq_along(varLabs)){
  formula <- as.formula(paste0(varLabs[v], " ~ 1 + day + factor(beep)"))
  lmRes <- lm(formula, data = Data5b)
  
  # Fixed effects:
  fixed <- coef(lmRes)
  
  # make zero if not significant at alpha:
  p_values[[varLabs[v]]] <- anova(lmRes)[["Pr(>F)"]][1:2]
  if (p_values[[varLabs[v]]][1] > alpha){
    fixed[2] <- 0
  }
  if (p_values[[varLabs[v]]][2] > alpha){
    fixed[3:5] <- 0
  }
  
  # Add to DFs:
  fitted_all[,varLabs[v]] <- fixed[1] + fixed[2] * fitted_all[["day"]]  +  fixed[3] * (fitted_all[["beep"]] == 1)  + 
  fixed[4] * (fitted_all[["beep"]] == 2) + fixed[5] *  (fitted_all[["beep"]] == 3)
  
  fitted_day[,varLabs[v]] <- fixed[1] + fixed[2] * fitted_day[["day"]]
  
  fitted_beep[,varLabs[v]] <- fixed[1] + fixed[2] * median(fitted_day[["day"]]) +  fixed[3] * (fitted_beep[["beep"]] == 1)  + 
    fixed[4] * (fitted_beep[["beep"]] == 2) + fixed[5] *  (fitted_beep[["beep"]] == 3)
  
  # Detrend data:
  data_detrended[,varLabs[v]] <- Data5b[,varLabs[v]] - (fixed[1] + fixed[2] * Data5b[["day"]]  +  fixed[3] * (Data5b[["beep"]] == 1)  + 
    fixed[4] * (Data5b[["beep"]] == 2) + fixed[5] *  (Data5b[["beep"]] == 3))
  
  ids <- rownames(anova(lmRes))
  testStatistics[[v]] <- cbind(data.frame(var = varLabs[v], effect = ids), anova(lmRes))
  
  coefficients[[v]] <- data.frame(
    var = varLabs[v],
    type = names(coef(lmRes)),
    coef = coef(lmRes),
    std = coef(lm.beta(lmRes))
  )
}



# -------------------------------------------------------------------------
# --------------- 4. Here we estimate network models ----------------------
# -------------------------------------------------------------------------

# Estimate network using multilevel VAR model
res <- mlVAR(data_detrended,
             vars=varLabs,
             idvar="id",
             dayvar="day",
             beepvar="beep",
             lags = 1,
             temporal = "orthogonal",
             contemporaneous = "orthogonal",
             nCores = 8)

# this is how we can save the object after the estimation; careful, over 100mb
# save(res, file=paste0(datapath, "network_orthogonal.RData"))

# you can later load it via
# load(paste0(datapath, "network_orthogonal.RData"))



# Plot :
names <- c("Relax","Irritable","Worry","Nervous","Future", "Anhedonia","Alone","Social-offline", "Social-online")

gr <- list('Stress'=c(1:6), 'Social'=c(7:9))

# Get networks:
cont <- getNet(res, "contemporaneous", layout = "spring", nonsig = "hide", rule = "and")
bet  <- getNet(res, "between", nonsig = "hide", rule = "and")
temp <- getNet(res, "temporal", nonsig = "hide")

L <- averageLayout(cont, temp)

pdf(paste0(figs, "figure.pdf"), width=6, height=2.5)
layout(matrix(c(1,1,2,2,2), nc=5, byrow = TRUE)) # 40% vs 60% widths
n1 <- qgraph(cont, layout = L,
       title="Contemporaneous network", theme='colorblind', negDashed=FALSE,
       groups=gr, legend=FALSE, nodeNames = names, labels=c(1:9),
       vsize=12,color=viridis_pal()(4)[3:4])
n2 <- qgraph(temp, layout = L,
       title="Temporal network", theme='colorblind', negDashed=FALSE, diag=FALSE,
       groups=gr, legend.cex=0.5, legend=TRUE, nodeNames = names, labels=c(1:9),
       vsize=10,color=viridis_pal()(4)[3:4], asize=6, curve=0.75, curveAll=T)
dev.off()
