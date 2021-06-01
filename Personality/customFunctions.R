#############################################################
# R functions developed for the paper                       #
# "Network approaches to the analysis of multivariate data  #
# to appear in Nature Methods Primer.                       #
# R code by Giulio Costantini (giulio.costantini@unimib.it) #
# and Marco Perugini                                        #
#############################################################




# Fade is a function that fade edges, that I have taken and very slightly edited
# from package qgraph
Fade <- function(col,alpha,bg)
{
  
  # col = color to fade
  # bg = color to fade to
  # alpha = inverse transparency, 1 = fully visible, 0 = fully transparent.
  
  if(is.na(col)) return(NA)
  
  if (missing(bg)) bg <- par("bg")
  if (length(bg)!=1) stop("'bg' must be of length 1")
  if (length(alpha)==1) alpha <- rep(alpha,length(col))
  if (length(col)==1) col <- rep(col,length(alpha))
  if (length(col)!=length(alpha)) stop("Length of 'col' not equal to length of 'alpha'")
  
  n <- length(col)
  
  rgbCols <- col2rgb(col)/255
  rgbBG <-  col2rgb(bg)/255
  
  colAlpha <- col2rgb(col,alpha=TRUE)[4,]/255
  
  Mix <- rgbCols*rep(alpha,each=3) + rgbBG%*%t(1-alpha)
  
  return(rgb(Mix[1,],Mix[2,],Mix[3,],colAlpha))
}

# function pcor2beta gives you data from a partial correlation/network
#-- input --
# pcor = a partial correlation matrix / network
# -- output --
# a matrix of betas, each column corresponds to a dependent variable
# so that you can get predicted values by a matrix multiplication
# in the form betas %*% data

pcor2beta <- function(pcor)
{
  require(psych)
  require(corpcor)
  diag(pcor) <- 1
  p <- ncol(pcor)
  betas <- matrix(0, ncol = p, nrow = p)
  
  # force the matrix to be symmetric.
  # sometimes numeric differences in the order of 1e-16 cause errors
  cm <- pcor2cor(pcor)
  cm[upper.tri(cm)] <- t(cm)[upper.tri(cm)]
  
  # this code was adapted from psych::matReg
  for(i in 1:p)
  {
    y = i
    x = seq(p)[-i]
    betas[-i, i] <- solve(cm[x, x], cm[x, y])
  }
    
  betas[abs(betas) < 1e-13] <- 0
  betas
}


# function R2 gives you two different types of R2 and of predicted values
#-- input --
# - betas = a matrix of standardized regression coefficients,
#    taken from a network
# - dt = a new data matrix
# - refit: logical, regulates whether R2_refit and predicted_refit are computed.
#    see below

# -- output --
# two types of R2 and predicted values can be returned:
# - R2_orig and predicted_orig use beta weights included in the beta matrix
# - R2_refit and predicted_refit use only the sparsity pattern of the network
#   and then refit linear regression using the same pattern of zeroes in the
#   network, but re-estimating regression coefficients with linear regression.

R2 <- function(betas, dt, refit = FALSE)
{
  # standardize data
  dt <- data.frame(scale(dt))
  
  out <- list()
  p <- ncol(betas)
  
  # predict values
  predicted <- as.matrix(dt) %*% betas
  
  # calculate R squared using the formula in Haslbeck & Waldorp (2018, p. 856)
  R2 <- 1 - apply(predicted - dt, 2, var)
  out$predicted_orig <- predicted
  out$R2_orig <- R2
  
  if(refit)
  {
    # refit regression models using the sparsity indicated in the matrix of
    # betas
    betas_refit <- matrix(0, ncol = p, nrow = p)
    for(i in 1:p)
    {
      if(any(betas[,i] != 0))
      {
        fit <- lm(dt[,i] ~ as.matrix(dt[,betas[,i] != 0]))
        betas_refit[betas[,i] != 0, i] <- fit$coefficients[-1]
      } else betas_refit[,i] <- 0
    }
    
    predicted <- as.matrix(dt) %*% betas_refit
    # R squared
    R2<- 1-apply(predicted - dt, 2,var)
    out$predicted_refit <- predicted
    out$R2_refit <- R2
  }
  out
}

# a wrapper that gives you directly the predictability from the network and
# new data
predictability <- function(net, dt, refit = FALSE)
{
  betas <- pcor2beta(net)
  predict <- R2(betas, dt = dt, refit = refit)
  predict
}


# references
# Haslbeck, J. M. B., & Waldorp, L. J. (2018). How well do network models predict observations? On the importance of predictability in network models. Behavior Research Methods, 50(2), 853-861. https://doi.org/10.3758/s13428-017-0910-x