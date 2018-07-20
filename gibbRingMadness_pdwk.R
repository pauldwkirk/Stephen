#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------

# Author: Stephen Coleman
# Subject: Gibbs sampling from a univariate Gaussian

# My first attempt at writing a Gibbs sampler in R. Based on Yee Whye Teh's
# MAtLab code from https://www.stats.ox.ac.uk/~teh/npbayes.html, although his
# approach is OO whereas I'm more likely to be Functional.

# ------------------------------------------------------------------------------

# === Libraries ================================================================

library(ramcmc) # install.packages("ramcmc", dep = T)
library(glue) # install.packages("glue", ep = T)

# === Functions ================================================================

gibbs_sampler <- function(fm, num_iter = 5) {
  
  # Number of samples (add_count creates a row "n", hence "N")
  # N <- nrow(data)
  
  k            <- fm$k
  N            <- fm$n
  alpha        <- fm$a
  q            <- fm$q
  data         <- fm$data
  class_labels <- fm$class_label
  class_table  <- fm$class_pop
  
  for (i in 1:num_iter) {
    # the number of desired iterations in our Gibbs sampler
    print(paste("Iteration number:", i))
    for (j in 1:N) {
      # for each entry
      print(paste("Individual number:", j))
      
      
      # Remove the current entry from the data and compare models
      # including and excluding it
      
      # Current entry and class
      curr_entry <- data[j] # univariate
      curr_class <- class_labels[j]
      
      if(j == 1)
      {
        print("debug")
      }
      # Update the data to reflect removing the current point
      q[[curr_class]] <- del_item(q[[curr_class]], curr_entry) 
      class_table[curr_class] <- class_table[curr_class] - 1
      
      # Calculate the probabilities for belonging to each class
      # Log, exp to handle awkward numbers
      prob <- log(class_table + alpha / k)
      
      if(any(is.nan(prob))){
        print("Calculated probs: ")
        print(prob)
      }
      
      print(glue("Curr entry value: {value}", value = curr_entry))
      

      # q[k] is a Gaussian, q is a mixture of gaussians
      for (l in 1:k) {
        print(glue("Curr variance value: {value}", value = q[[l]]$C))
        prob[l] <- prob[l] + log_predictive(q[[l]], curr_entry)
      }
      
      print("hi")
      
      prob <- exp(prob - max(prob))
      prob <- prob / sum(prob)
      
      # Rnadom number used to assign class in this iteration
      u <- runif(1)
      
      # Assign to class based on cumulative probabilities
      pred <- 1 + sum(u > cumsum(prob)) #
  
      print(table(class_labels))
      print("Numbers in each class:")
      print(q[[1]]$n)
      print(q[[2]]$n)
      
          
      # add the current entry back into model (component q[k])
      class_labels[j] <- pred
      class_table[pred] <- class_table[pred] + 1
      q[[pred]] <- add_item(q[[pred]], data[j])
      
      print(table(class_labels))
      print("Numbers in each class:")
      print(q[[1]]$n)
      print(q[[2]]$n)

      
      print(class_table)
      
    }
  }
  
  fm$q <- q
  fm$class_label <- class_labels
  fm$class_pop <- class_table
  
  return(fm)
}

fm_init <- function(k, a, q0, data, class_label) {
  # initialize finite mixture model, with
  # k mixture components,
  # a concentration parameter,
  # q0 an empty component with h prior,
  # data
  # class_label initial cluster assignments (between 1 and k).
  
  fm             <- list()
  fm$k           <- k
  fm$n           <- length(data)
  fm$a           <- a
  fm$q           <- list()
  fm$data        <- data
  fm$class_label <- class_label
  fm$class_pop   <- rep(0, k)
  
  pop <- rep(0, k)
  
  for (class in 1:k) {
    fm$q[[class]] <- q0
  }
  
  for (i in 1:fm$n) {
    class               <- class_label[i]
    fm$q[[class]]       <- add_item(fm$q[[class]], data[i])
    fm$class_pop[class] <- fm$class_pop[class] + 1
  }
  return(fm)
}

Gaussian <- function(d, var_coef, df, mean_prior, cluster_cov) {
  # Returns a list of the various components required of a Gaussian distn
  q           <- list()
  precision   <- 1 / var_coef
  S           <- cluster_cov * df
  
  q$d         <- d # Dimensionality of distn
  q$n         <- 0 # Number of members in distn (initialised to 0)
  q$precision <- precision
  q$df        <- df # degress of freedom
  q$C         <- chol(S + precision * mean_prior * mean_prior) # Cholesky factorisation of covariance
  q$X         <- precision * mean_prior # Some score thing? Not sure how this combines with a new point
  q$Z0        <- evidence(q$d, q$n, q$precision, q$df, q$C, q$X) # could use q alone
  
  return(q)
}

rand <- function(q) {
  # I'm not sure this is ever used
  C <- ramcmc::chol_downdate(q$C, q$X) # chol(q$C - (q$X/q$precision) * t(q$X/q$precision))
  C <- mldivide(C, diag(q$d), pinv = TRUE) # left division of matrix C by diag(dd)
}

# chol_update <- function(C, X, sign = "+") {
#   # Update a cholesky factorisation for C using X
#   if (!(sign %in% c("+", "-"))) stop("sign must be one of '+' or '-' ")
#   Y <- as.matrix(X)
#   if (sign == "+") {
#     #new_C <- chol(C + (Y %*% t(Y)))
#     new_C  <- ramcmc::chol_update(C, X)
#   } else {
#     #new_C <- chol(C - (Y %*% t(Y)))
#     new_C  <- ramcmc::chol_downdate(C, X)
#   }
# }

log_predictive <- function(q, new_point) {
  # Compares the evidence of the distribution q with the new_point included and
  # excluded
  if(q$n != 0)
  {
  lp <- (
    evidence(
      q$d,
      q$n + 1,
      q$precision + 1,
      q$df + 1,
      ramcmc::chol_update(q$C, new_point),
      q$X + new_point
    )
    - evidence(
      q$d,
      q$n,
      q$precision,
      q$df,
      q$C,
      q$X
    )
  )
  }
  else
  {
    lp <- (
      evidence(
        q$d,
        q$n + 1,
        q$precision + 1,
        q$df + 1,
        ramcmc::chol_update(q$C, new_point),
        q$X + new_point
      )
    )
  }
}

evidence <- function(d, n, precision, df, C, X) {
  # Computes the normalising constant for given inputs
  marginal_likelihood <- (
    -n * d / 2 * log(pi)
    - d / 2 * log(precision)
    - df * sum(log(diag(ramcmc::chol_downdate(C, X / sqrt(precision))))) # not sure about this line
    + sum(lgamma((df - (0:(d - 1)) )/ 2))
  )
}


add_item <- function(q, x) {
  # Add data point x to distribution q and update contained information
  # accordingly
  q$n <- q$n + 1
  q$precision <- q$precision + 1
  q$df <- q$df + 1
  q$C <- ramcmc::chol_update(q$C, x)
  q$X <- q$X + x
  return(q)
}

del_item <- function(q, x){
  q$n <- q$n - 1
  q$precision <- q$precision - 1
  q$df <- q$df - 1
  q$C <- ramcmc::chol_downdate(q$C, x)
  q$X <- q$X - x
  return(q)
}

# === Demo 1 ===================================================================

set.seed(1)
# demo of finite mixture model in 1D
dd <- 1
KK <- 2
NN <- 10
# xx = [-2+.5*randn(1,5) 2+.5*randn(1,5)];
# xx <- c(rnorm(NN / 2, -2, 1), rnorm(NN / 2, 2, 1))
xx <- unlist(read.csv(file = "dataFile.csv", header= F))  
aa <- 1
s0 <- 3
ss <- 1
numiter <- 1000
hh.dd <- dd
hh.ss <- s0^2 / ss^2
hh.vv <- 5
hh.VV <- ss^2 * diag(dd)
# hh.uu <- matrix(0, nrow = dd, ncol = 1)
hh.uu <- rep(0, dd)

# % set range over which to evaluate density
yy <- seq(-15, 15, by = 0.01)

# % initialize records
record.KK <- matrix(0, nrow = numiter, ncol = length(yy))
record.p0 <- matrix(0, nrow = numiter, ncol = length(yy)) # densities
record.pp <- matrix(0, nrow = numiter, ncol = length(yy)) # densities

# initialize finite mixture with data items.
# xx = num2cell(xx)
zz <- sample(1:KK, NN, replace = TRUE)
q0 <- Gaussian(hh.dd, hh.ss, hh.vv, hh.VV, hh.uu)

fm0 <- fm_init(2, aa, q0, xx, zz)


# This is what we want to run:

# --- Gibbs --------------------------------------------------------------------

fm <- gibbs_sampler(fm0, num_iter = 5)

# Expect problems

# % run
# for (iter in 1:numiter) {
#   # fprintf(1, "finite mixture: iter# %d\r", iter)
#
#   # % gibbs iteration
#   fm <- gibbs_sampler(fm0, num_iter = 5)
# }

# Sample of not working. Go through single iteration of Gibbs sampler step by step
fm           <- fm0
k            <- fm$k
N            <- fm$n
alpha        <- fm$a
q            <- fm$q
data         <- fm$data
class_labels <- fm$class_label
class_table  <- fm$class_pop

num_iter <- 1 # first iteration of gibbs_sampler UDF
j <- 1 # i.e. first data point
# for each entry

# Remove the current entry from the data and compare models
# including and excluding it

# Current entry and class
curr_entry <- data[j] # univariate
curr_class <- class_labels[j]

# Update the data to reflect removing the current point
q[[curr_class]] <- del_item(q[[curr_class]], curr_entry)
class_table[curr_class] <- class_table[curr_class] - 1

# Calculate the probabilities for belonging to each class
# Log, exp to handle awkward numbers
prob <- log(class_table + alpha / k)


# q[k] is a Gaussian, q is a mixture of gaussians

# --- Breakpoint ---------------------------------------------------------------
for (l in 1:k) {
  prob[l] <- prob[l] + log_predictive(q[[l]], curr_entry)
}

prob <- exp(prob - max(prob))
prob <- prob / sum(prob)

# tends to break before here, if interested this is just the inners of
# the UDF gibbs_sampler

# Rnadom number used to assign class in this iteration
u <- runif(1)

# Assign to class based on cumulative probabilities
pred <- 1 + sum(u > cumsum(prob)) #

# add the current entry back into model (component q[k])
class_labels[j] <- pred

class_table[pred] <- class_table[pred] + 1
q[[k]] <- add_item(q[[k]], data[j])


q
