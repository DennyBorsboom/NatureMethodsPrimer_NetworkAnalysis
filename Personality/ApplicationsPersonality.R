# R code for the cross-sectional personalty example in      #
# "Network approaches to the analysis of multivariate data  #
# to appear in Nature Methods Primer.                       #
# R code by Giulio Costantini (giulio.costantini@unimib.it) #
# and Marco Perugini                                        #


#### data and codebook ####

dt <- read.csv("PersonalityData.csv")

# The dataset includes a subset of the data presented in Study 3 in
# Costantini, Saraulli & Perugini (2020)

# codebook
# - CON = conscientiousness score
# - ORD, IND, IMC = conscientiousness facet scores, respectively orderliness,
#      industriousness and impulse-control
# -precise, ..., imprudent = conscientiousness adjectives
# -G12, ..., G26: Conscientious and unconscientious goals
#   * G12 = Do something well, avoid mistakes
#   * G10 = Have control
#   * G13 = Accomplish something, observe a commitment
#   * G11 = Personal realization
#   * G16 = Avoid or manage things you do not care about
#   * G25 = Think, reflect
#   * G08 = Comply with rules
#   * G17 = Be safe
#   * G26 = Do not think


#### load packages and functions ####

if(!require(pacman)) install.packages("pacman")
pacman::p_load("dplyr", "reshape2", "qgraph", "stringr", "bootnet", "reshape2",
               "ggplot2", "ggthemes", "ggpubr", "psych", "corpcor",  "parallel")

#### load custom functions ####
source("customFunctions.R")

#### define some useful variables ####

# define a color palette for nodes to use often in code below
col1 <- "#FDE624" # goals
col2 <- "#218F8C" # trait
col3 <- "#99CC00" # facet orderliness
col4 <- "#339900" # facet industriousness
col5 <- "#33CC33" # facet impulse-control

# plot margins for figures
mar <- c(1, 1, 1, 1)


#### estimate trait-level network ####

# select data
netdt1 <- select(dt,
                 G12:G26,
                 CON)

# estimate network
net1 <- EBICglasso(cor(netdt1), n = nrow(netdt1))

#### visualize trait-level network ####

# simple network visualization, using a default circular layout
qgraph(net1)

# the following lines in this section make a "prettier" visualization presented
# in the paper, that highlights the personality-goal connections. This plot
# (as well the corresponding ones for facets and items) are obtained by
# overlapping two qgraph plots each, one with faded within-construct edges and
# the other with full between-construct edges

# define graphical parameters
clr <- c(rep(col1, 9), col2)
lbl <- c(1:10)

# define a convenient circular layout
layout <- matrix(0, nrow = 10, ncol = 2)
layout[1:9, 1] = sin(seq(0, 2 * pi, length = 10))[-10]
layout[1:9, 2] = cos(seq(0, 2 * pi, length = 10))[-10]
layout <- layout * .9

# generate the initial qgraph plot without visualizing it
qg1 <- qgraph(net1, color = clr, labels = lbl, layout = layout, vsize = 18,
              maximum = .4, fade = FALSE, DoNotPlot = TRUE,
              theme = "colorblind")

# fade colors within constructs to highlight personality-goal connections
edg <- matrix(NA, nrow = nrow(net1), ncol = ncol(net1))
colnames(edg) <- rownames(edg) <- colnames(net1)
edg[upper.tri(edg)][net1[upper.tri(net1)] != 0] <- 
  qg1$graphAttributes$Edges$color
edg[lower.tri(edg)] <- t(edg[upper.tri(edg)])
edg[,] <-  sapply(edg, Fade, bg = "white", alpha = .2)

# prepare plot for within-constructs transparent edges
qg1 <- qgraph(net1,
              color = edg,
              labels = "", vsize = 0, layout = layout, 
              maximum = .4, edge.color = edg, 
              fade = FALSE, DoNotPlot = TRUE,
              title = "(a) Trait-level",
              title.cex = 2,
              mar = mar, rescale = FALSE)

# prepare plot for edges connecting traits and goals
net1b <- net1
net1b[1:9, 1:9] <- 0

qg1a <- qgraph(net1b, layout = qg1$layout, labels = lbl, color = clr,
               vsize = 18, fade = FALSE, plot = FALSE, theme = "colorblind",
               mar = mar, DoNotPlot = TRUE, rescale = FALSE, 
               label.cex = 1.8)


plot(qg1)
plot(qg1a)


#### estimate facet-level network ####

# select data
netdt2 <- select(dt, 
                 G12:G26,
                 ORD:IMC)
# estimate network
net2 <- EBICglasso(cor(netdt2), n = nrow(netdt2))

#### visualize facet-level network ####
# simple visualization
qgraph(net2)

# nicer visualization

# define graphical parameters
clr <- c(rep(col1, 9),
         col3,
         col4,
         col5)
lbl <- c(1:9, 11:13)

# define a convenient circular layout
maincircX <- sin(seq(0, 2 * pi, length = 4))[-4]/2.5
maincircY <- cos(seq(0, 2 * pi, length = 4))[-4]/2.5

layout <- matrix(0, nrow = 12, ncol = 2)
layout[1:9, 1] = sin(seq(0, 2 * pi, length = 10))[-10]
layout[1:9, 2] = cos(seq(0, 2 * pi, length = 10))[-10]
layout[10:12, 1] = maincircX
layout[10:12, 2] = maincircY
layout <- layout * .9

# generate the initial qgraph plot
qg2 <- qgraph(net2, color = clr, labels = lbl, layout = layout, maximum = .4,
              vsize = 18, fade = FALSE, DoNotPlot = TRUE, theme = "colorblind")


# fade colors within constructs
edg <- matrix(NA, nrow = nrow(net2), ncol = ncol(net2))
colnames(edg) <- rownames(edg) <- colnames(net2)
edg[upper.tri(edg)][net2[upper.tri(net2)] != 0] <- qg2$graphAttributes$Edges$color
edg[lower.tri(edg)] <- t(edg[upper.tri(edg)])
edg[,] <-  sapply(edg, Fade, bg = "white", alpha = .2)


# prepare plot for within-constructs transparent edges
qg2 <- qgraph(net2,
              color = edg,
              labels = "", vsize = 0, layout = layout, 
              maximum = .4, edge.color = edg, 
              fade = FALSE, DoNotPlot = TRUE,
              theme = "colorblind",
              title = "(b) Facet-level",
              title.cex = 2,
              mar = mar,
              rescale = FALSE
)

# prepare plot for edges connecting traits and goals
net2b <- net2
net2b[1:9, 1:9] <- 0
net2b[10:nrow(net2b), 10:ncol(net2b)] <- 0

qg2a <- qgraph(net2b, layout = qg2$layout, labels = lbl, color = clr,
               vsize = 18, fade = FALSE, theme = "colorblind",
               plot = FALSE,
               mar = mar, DoNotPlot = TRUE,
               rescale = FALSE, 
               label.cex = 1.8)

plot(qg2)
plot(qg2a)


#### estimate item-level network ####

# select data
netdt3 <- select(dt,
                 G12:G26,
                 precise:imprudent)

# estimate network
net3 <- EBICglasso(cor(netdt3), n = nrow(netdt3), threshold = FALSE)

#### visualize item-level network ####

# simple visualization
qgraph(net3)

# nicer visualization

# define graphical parameters
clr <- c(rep(col1, 9),
         rep(col3, 10),
         rep(col4, 10),
         rep(col5, 10))
lbl <- c(1:9, 14:43)

vsize <- c(rep(18, 9), rep(11, 30))

# define a convenient circular layout
maincircX <- sin(seq(0, 2 * pi, length = 4))[-4]/2.5
maincircY <- cos(seq(0, 2 * pi, length = 4))[-4]/2.5
smallcircX <- sin(seq(0, 2 * pi, length = 11))[-11]/3.8
smallcircY <- cos(seq(0, 2 * pi, length = 11))[-11]/3.8

layout = matrix(0, nrow = ncol(net3), ncol = 2)
layout[1:9, 1] = sin(seq(0, 2 * pi, length = 10))[-10]
layout[1:9, 2] = cos(seq(0, 2 * pi, length = 10))[-10]
layout[10:19, 1] = maincircX[1] + smallcircX
layout[10:19, 2] = maincircY[1] + smallcircY
layout[20:29, 1] = maincircX[2] + smallcircX
layout[20:29, 2] = maincircY[2] + smallcircY
layout[30:39, 1] = maincircX[3] + smallcircX
layout[30:39, 2] = maincircY[3] + smallcircY
layout <- layout * .9

# generate the initial qgraph plot
qg3 <- qgraph(net3, layout = layout, labels = lbl, color = clr,
              maximum = .4, vsize = vsize, cut = 0,
              fade = FALSE, DoNotPlot = TRUE,
              theme = "colorblind")

# fade colors within constructs
edg <- matrix(NA, nrow = nrow(net3), ncol = ncol(net3))
colnames(edg) <- rownames(edg) <- colnames(net3)
edg[upper.tri(edg)][net3[upper.tri(net3)] != 0] <- qg3$graphAttributes$Edges$color
edg[lower.tri(edg)] <- t(edg[upper.tri(edg)])

edg[,] <-  sapply(edg, Fade, bg = "white", alpha = .2)


# prepare plot for within-constructs transparent edges
qg3 <- qgraph(net3,
              color = edg,
              labels = "", vsize = 0, layout = layout, 
              maximum = .4, edge.color = edg,
              fade = FALSE, DoNotPlot = TRUE,
              theme = "colorblind",
              title = "(c) Item-level",
              title.cex = 2,
              mar = mar,
              rescale = FALSE)

# prepare plot for edges connecting traits and goals
net3b <- net3
net3b[1:9, 1:9] <- 0
net3b[10:nrow(net3b), 10:ncol(net3b)] <- 0

qg3a <- qgraph(net3b, layout = qg3$layout, labels = lbl, color = clr,
               vsize = vsize, fade = FALSE, plot = FALSE,
               theme = "colorblind", DoNotPlot = TRUE,
               mar = mar, DoNotPlot = TRUE,
               rescale = FALSE, 
               label.cex = 1.8)


plot(qg3)
plot(qg3a)

#### network analysis ####

# one can compute Strength centrality and visualize it using a simple built-in
# function, centralityPlot, in package qgraph.
# The following three lines visualize centrality estimates in each network
# separately.

centralityPlot(net1, include = "Strength", standardized = FALSE)
centralityPlot(net2, include = "Strength", standardized = FALSE)
centralityPlot(net3, include = "Strength", standardized = FALSE)


# the following code within this section produces the nicer visualization
# that we used in the paper, that only represents centrality for goals.

# compute centrality estimates
str1 <- centrality_auto(net1)$node.centrality[,"Strength"] %>%
  as.vector()
names(str1) <- colnames(net1)

str2 <- centrality_auto(net2)$node.centrality[,"Strength"] %>%
  as.vector()
names(str2) <- colnames(net2)

str3 <- centrality_auto(net3)$node.centrality[,"Strength"] %>%
  as.vector()
names(str3) <- colnames(net3)


# merge centrality estimates in a dataframe
cnt <- merge(str1[1:9], str2[1:9], by = 0, all = TRUE) %>%
  rename(node = Row.names, trait = x, facet = y) %>%
  merge(str3[1:9], by.x = "node", by.y = 0, all = TRUE) %>%
  rename(item = y)

cnt$node <- factor(cnt$node, levels = c(
  "G08", "G10", "G11", "G12", "G13", "G16", "G17", "G25", "G26"),
  labels = 1:9)

cnt_mlt <- melt(cnt, id.vars = "node",
                variable.name = "level",
                value.name = "centrality")
cnt_mlt <- arrange(cnt_mlt, level, node)

# generate the plot using package ggplot2
cntPlot <- ggplot(cnt_mlt,
                  aes(x = centrality,
                      y = node,
                      group = level,
                      shape = level)) +
  geom_point(aes(size = 2)) + 
  geom_path(color = "grey") +
  guides(color = FALSE, size = FALSE) +
  ylab("goal") +
  ggtitle("(d) Centrality of goals") +
  guides(shape = guide_legend(override.aes = list(size = 5))) +
  theme_few() +
  theme(legend.position = c(.85, .9),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        plot.title = element_text(hjust = 0, size = 18)
        ) 
cntPlot


#### testing out-of-sample predictions ####

# It is interesting to inspect how network models estimated at different levels
# (trait, facets, items) perform in the prediction of goals. We examined the
# out-of-sample predictability of each goal (Haslbeck & Waldorp, 2017) using a 
# 10-fold crossvalidation approach (Yarkoni & Westfall, 2017).
# This analysis estimates out-of-sample R-squared estimates for each goal
# in each network model

K <- 10 # number of cross-validation folds
N <- nrow(dt) # number of observations

set.seed(1)
# randomly assign each case to a fold
rnd <- rep(1:K, length.out = N) %>% 
  sample() 

# split the data according to the folds. Each time, 9/10ths of the sample
# serve as training sample (train_smpl) and the remaining 1/10th serves
# as validation sample (valid_smpl)
smpl <- sapply(1:K, function(x) which(rnd == x))
train_smpl <- lapply(smpl, function(x) dt[-x,])
valid_smpl <- lapply(smpl, function(x) dt[x,])

# create a list to store predicted values in of the three networks
# trait = net1
# facet = net2
# item = net3
PRD <- list()
PRD[["trait"]] <- list()
PRD[["facet"]] <- list()
PRD[["item"]] <- list()

# Run cross validation for each dependent variable and for each model
for(i in 1:K) # loop for each cross-validation fold
{
  # trait network
  netdt1T <- select(train_smpl[[i]],
                   G12:G26,
                   CON)
  net1 <- EBICglasso(cor(netdt1T), n = nrow(netdt1T))
  
  prd1 <- predictability(net1, dt = select(valid_smpl[[i]],
                                           G12:G26,
                                           CON))
  
  PRD$trait[[i]] <- prd1$predicted_orig[,1:9]
  
  # facet network
  netdt2T <- select(train_smpl[[i]], 
                   G12:G26,
                   ORD:IMC)
  
  net2 <- EBICglasso(cor(netdt2T), n = nrow(netdt2T))
  
  prd2 <- predictability(net2, dt = select(valid_smpl[[i]], 
                                           G12:G26,
                                           ORD:IMC))
  
  PRD$facet[[i]] <- prd2$predicted_orig[,1:9]
  
  # item network
  netdt3T <- select(train_smpl[[i]],
                   G12:G26,
                   precise:imprudent)
  
  net3 <- EBICglasso(cor(netdt3T), n = nrow(netdt3T))
  
  prd3 <- predictability(net3, dt = select(valid_smpl[[i]],
                                           G12:G26,
                                           precise:imprudent))
  
  PRD$item[[i]] <- prd3$predicted_orig[,1:9]
}

# define prediced values
PRD_trait <- do.call(rbind, PRD$trait)
PRD_facet <- do.call(rbind, PRD$facet)
PRD_item <- do.call(rbind, PRD$item)

# define observed values (in the same order of predicted)
valid_smpl2 <- lapply(valid_smpl, function(x) select(x, G12:G26) %>% scale())
OBS <- do.call(rbind, valid_smpl2)

# calculate errors
ERR_trait <- OBS - PRD_trait
ERR_facet <- OBS - PRD_facet
ERR_item <- OBS - PRD_item

# estimate R-squared
R2_trait <- 1 - colSums(ERR_trait^2)/colSums(OBS^2)
R2_facet <- 1 - colSums(ERR_facet^2)/colSums(OBS^2)
R2_item <- 1 - colSums(ERR_item^2)/colSums(OBS^2)
R2_tot <- rbind(R2_trait, R2_facet, R2_item)

# Table reported in the paper
# Overall, results indicate that items show slightly better performance in terms
# of prediction closely followed by facets. Trait-conscientiousness is never the
# best in terms of performance
R2_tot[,order(colnames(R2_tot))] %>% round(3)


#### stability and robustness ####

# NOTE: since the code that follows involve bootstrapping, executing this part
# might take a relatively long time, particularly on older computers.
# the code takes advantage of multiple core architectures. The value ncores here
# is automatically set ot the number of available cores minus 1. You can also
# set it to a specific value (e.g., by replacing the line of code below with
# "ncores <- 8")

ncores <- detectCores()-1

# nonparametric bootstrap for stability of edges and of edge differences
# only the one for network 1 is reported in the paper

set.seed(1)
NP_Boot1 <-  bootnet(netdt1,
                     default = "EBICglasso",
                     type = "nonparametric",
                     nCores = ncores,
                     computeCentrality = FALSE,
                     statistics = "edge")

NP_Boot2 <-  bootnet(netdt2,
                     default = "EBICglasso",
                     type = "nonparametric",
                     nCores = ncores,
                     computeCentrality = FALSE,
                     statistics = "edge")
NP_Boot3 <-  bootnet(netdt3,
                     default = "EBICglasso",
                     type = "nonparametric",
                     nCores = ncores,
                     computeCentrality = FALSE,
                     statistics = "edge")

# visualize results

# CI around edges
plot(NP_Boot1, plot = "area", order = "sample", legend = FALSE)
# differences between edges
plot(NP_Boot1, plot = "difference", order = "sample",
                onlyNonZero = FALSE, labels = FALSE)

plot(NP_Boot2, plot = "area", order = "sample", legend = FALSE)
plot(NP_Boot2, plot = "difference", order = "sample",
     onlyNonZero = FALSE, labels = FALSE)

plot(NP_Boot3, plot = "area", order = "sample", legend = FALSE)
plot(NP_Boot3, plot = "difference", order = "sample",
     onlyNonZero = FALSE, labels = FALSE)


# Replication Simulator analysis
set.seed(1)
Sim1 <- netSimulator(
  input = net1, 
  dataGenerator = ggmGenerator(),
  nCases = seq(200, 1000, 200),
  nCores = ncores,
  nReps = 100,
  default = "EBICglasso")

Sim2 <- netSimulator(
  input = net2, 
  dataGenerator = ggmGenerator(),
  nCases = seq(200, 1000, 200),
  nCores = ncores,
  nReps = 100,
  default = "EBICglasso")

Sim3 <- netSimulator(
  input = net3, 
  dataGenerator = ggmGenerator(),
  nCases = seq(200, 1000, 200),
  nCores = ncores,
  nReps = 100,
  default = "EBICglasso")

# the results can be easily visualized as three separate figures using the
# default plot method
plot(Sim1)
plot(Sim2)
plot(Sim3)

# the following code combines them in a single figure, like we in the paper
# using ggplot2
Sim1_mlt <- select(Sim1, nCases, sensitivity:correlation) %>%
  melt(id.vars = "nCases") %>%
  mutate(network = "trait")

Sim2_mlt <- select(Sim2, nCases, sensitivity:correlation) %>%
  melt(id.vars = "nCases") %>%
  mutate(network = "facet")

Sim3_mlt <- select(Sim3, nCases, sensitivity:correlation) %>%
  melt(id.vars = "nCases") %>%
  mutate(network = "item")

Sim_mlt <- rbind(Sim1_mlt, Sim2_mlt, Sim3_mlt)
Sim_mlt$network <- factor(Sim_mlt$network, levels = c("item", "facet", "trait"))

# colorblind palette
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73")


repSim_plt <- ggplot(Sim_mlt, aes(x = factor(nCases), y = value, fill = network)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  xlab("sample size") + 
  ylab("") +
  scale_colour_manual(values = cbbPalette) +
  scale_fill_manual(values = cbbPalette) +
  theme_bw() + 
  theme(legend.text = element_text(size = 20), 
        legend.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size = 18, angle = 90),
        axis.title.x =  element_text(size = 20),
        legend.key.size = unit(1, "in"),
        axis.text.y = element_text(size = 15)
        )


repSim_plt



# case-dropping bootstrap for centrality stability in the three networks
# dropping from 5% to 95% of cases in steps of 5%, with 1000 bootstrap resamples
# for each step.

set.seed(1)
CD_Boot1 <-  bootnet(netdt1,
                     default = "EBICglasso",
                     type = "case",
                     nCores = ncores,
                     computeCentrality = TRUE,
                     statistics = c("strength"),
                     nBoots = 19000,
                     caseMin = .05,
                     caseMax = .95,
                     caseN = 19
)

CD_Boot2 <-  bootnet(netdt2,
                     default = "EBICglasso",
                     type = "case",
                     nCores = ncores,
                     computeCentrality = TRUE,
                     statistics = c("strength"),
                     nBoots = 19000,
                     caseMin = .05,
                     caseMax = .95,
                     caseN = 19
)

CD_Boot3 <-  bootnet(netdt3,
                     default = "EBICglasso",
                     type = "case",
                     nCores = ncores,
                     computeCentrality = TRUE,
                     statistics = c("strength"),
                     nBoots = 19000,
                     caseMin = .05,
                     caseMax = .95,
                     caseN = 19
)

# estimate Correlation stability coefficients for strength centrality
corStability(CD_Boot1)
corStability(CD_Boot2)
corStability(CD_Boot3)

# visualize the results of case-dropping bootstrap
plot(CD_Boot1)
plot(CD_Boot2)
plot(CD_Boot3)


#### References mentioned in comments ####

# Costantini, G., Saraulli, D., & Perugini, M. (2020). Uncovering the motivational core of traits: The case of conscientiousness. European Journal of Personality. https://doi.org/10.1002/per.2237
# Haslbeck, J. M. B., & Waldorp, L. (2017). How well do network models predict future observations? On the importance of predictability in network models. Behavioral Research Methods, 1-9. https://doi.org/10.3758/s13428-017-0910-x
# Yarkoni, T., & Westfall, J. (2017). Choosing Prediction Over Explanation in Psychology: Lessons From Machine Learning. Perspectives on Psychological Science, 12(6), 1100-1122. https://doi.org/10.1177/1745691617693393

