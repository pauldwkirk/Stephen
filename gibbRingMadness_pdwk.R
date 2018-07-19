#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------

# Author: Stephen Coleman
# Subject: Gibbs sampling from a univariate Gaussian

# My first attempt at writing a Gibbs sampler in R. Based on Yee Whye Teh's
# MAtLab code from https://www.stats.ox.ac.uk/~teh/npbayes.html, although his
# approach is OO whereas I'm more likely to be Functional.

# ------------------------------------------------------------------------------

# === Functions ================================================================

gibbs_sampler <- function(fm, num_iter = 1000) {
  
  # Number of samples (add_count creates a row "n", hence "N")
  # N <- nrow(data)
  
  k            <- fm$k
  N            <- fm$n
  alpha        <- fm$a
  q            <- fm$q
  data         <- fm$data
  class_labels <- fm$class_label
  class_table  <- fm$class_pop
  
  # class_table <- rep(0, k)
  # for (i in 1:k) {
  #   class_table[i] <- length(class_labels[class_labels == i])
  # }
  # i <- 1
  # j <- 1
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
      
      # Update the data to reflect removing the current point
      #data                    <- data[-j]
      class_table[curr_class] <- class_table[curr_class] - 1
      
      # Calculate the probabilities for belonging to each class
      # Log, exp to handle awkward numbers
      prob <- log(class_table + alpha / k)
      if(any(is.nan(prob))){
        print("Calculated probs: ")
        print(prob)
      }
      # q[k] is a Gaussian, q is a mixture of gaussians
      for (l in 1:k) {
        prob[l] <- prob[l] + logpredictive(q[[l]], curr_entry)
      }
      
      prob <- exp(prob - max(prob))
      prob <- prob / sum(prob)
      
      # Rnadom number used to assign class in this iteration
      u <- runif(1)
      
      # Assign to class based on cumulative probabilities
      pred <- 1 + sum(u > cumsum(prob)) #
      
      # add the current entry back into model (component q[k])
      class_labels <- c(
        class_labels[1:j - 1],
        pred,
        class_labels[j:length(class_labels)]
      )
      
      class_table[pred] <- class_table[pred] + 1
      q[[k]] <- add_item(q[[k]], data[j])
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
  C <- chol_update(q$C, q$X, sign = "-") # chol(q$C - (q$X/q$precision) * t(q$X/q$precision))
  C <- mldivide(C, diag(q$d), pinv = TRUE) # left division of matrix C by diag(dd)
}

chol_update <- function(C, X, sign = "+") {
  # Update a cholesky factorisation for C using X
  if (!(sign %in% c("+", "-"))) stop("sign must be one of '+' or '-' ")
  Y <- as.matrix(X)
  if (sign == "+") {
    #new_C <- chol(C + (Y %*% t(Y)))
    new_C  <- ramcmc::chol_update(C, X)
  } else {
    #new_C <- chol(C - (Y %*% t(Y)))
    new_C  <- ramcmc::chol_downdate(C, X)
  }
}

evidence <- function(d, n, precision, df, C, X) {
  # Computes the normalising constant for given inputs
  marginal_likelihood <- (
    -n * d / 2 * log(pi)
    - d / 2 * log(precision)
    - df * sum(log(diag(chol_update(C, X / sqrt(precision), sign = "-")))) # not sure about this line
    + sum(lgamma((df - (0:(d - 1)) )/ 2))
  )
}

logpredictive <- function(q, new_point) {
  # Compares the evidence of the distribution q with the new_point included and
  # excluded
  lp <- (
    evidence(
      q$d,
      q$n + 1,
      q$precision + 1,
      q$df + 1,
      chol_update(q$C, new_point),
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
add_item <- function(q, x) {
  # Add data point x to distribution q and update contained information
  # accordingly
  q$n <- q$n + 1
  q$precision <- q$precision + 1
  q$df <- q$df + 1
  q$C <- chol_update(q$C, x)
  q$X <- q$X + x
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
data <- data[-j]
class_table[curr_class] <- class_table[curr_class] - 1

# Calculate the probabilities for belonging to each class
# Log, exp to handle awkward numbers
prob <- log(class_table + alpha / k)


# q[k] is a Gaussian, q is a mixture of gaussians

# --- Breakpoint ---------------------------------------------------------------
for (l in 1:k) {
  prob[l] <- prob[l] + logpredictive(q[[l]], curr_entry)
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
class_labels <- c(
  class_labels[1:j - 1],
  pred,
  class_labels[j:length(class_labels)]
)

class_table[pred] <- class_table[pred] + 1
q[[k]] <- add_item(q[[k]], data[j])


# --- Logpredictive break ------------------------------------------------------

# Closer inspection of logpredictive
q <- fm0$q
l <- 1

# sometimes this works - I think it depends on data generated
logpredictive(q[[l]], curr_entry)

q[[l]]$C
q4 <- chol_update(q[[l]]$C, curr_entry)

# issue is normally in this calculation
x1 <- evidence(
  q[[l]]$d,
  q[[l]]$n + 1,
  q[[l]]$precision + 1,
  q[[l]]$df + 1,
  chol_update(q[[l]]$C, curr_entry),
  q[[l]]$X + curr_entry
)

# specifically we cannot guarantee that q$C > XtX
# for the update in Cholesky factorisation
# We require this for positive definite, as
# q$C <- q$C - q$X %*% t(q$X)

# break calculation of x1 into component parts
new_C <- chol_update(q[[l]]$C, curr_entry)
new_X <- q[[l]]$X + curr_entry

# This is normally the problem point I think
probmat <- chol_update(new_C, new_X / (sqrt(q[[l]]$precision) + 1), sign = "-")

# Sometimes it is this though; as for x1 but excluding new point
retro <- chol_update(q[[l]]$C, q[[l]]$X / sqrt(q[[l]]$precision), sign = "-")

x2 <- evidence(
  q[[l]]$d,
  q[[l]]$n,
  q[[l]]$precision,
  q[[l]]$df,
  q[[l]]$C,
  q[[l]]$X
)
