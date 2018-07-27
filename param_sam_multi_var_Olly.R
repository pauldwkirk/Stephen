#!/usr/bin/env Rscript

# === Libraries ================================================================

library(expm) # install.packages("expm", dep = T)

# string manipulation (not needed, merely for print statements)
library(glue) # install.packages("glue", dep = T)

# generally useful - only using ggplot2 currently
library(tidyverse) # install.packages("tidyverse", dep = T)
library(MASS)

# for pretty heatmaps
library(pheatmap) # install.packages("pheatmap", dep = T)

# for multivariate t distn
library(LaplacesDemon) # install.packages("LaplacesDemon", dep = T)

# for colouring pheatmap
library(RColorBrewer)

# an alternative colour schema
library(ghibli) # install.packages("ghibli")

# for inverse Wishart (riwish)
library(MCMCpack) # install.packages("MCMCpack", dep = T)

# for Olly's stuff
require(pRoloc)
require(pRolocdata)
# require(Rtsne) # for t-SNE in plot2D

# for %<>%
library(magrittr)

# === Functions ================================================================

# --- Generic functions --------------------------------------------------------
data_frame_mean <- function(data) {
  # Calculates the mean of each column in a dataframe (for C++ transfer)
  # data: dataframe of numerical variables
  num_cols <- dim(data)[2]
  mean <- rep(0, num_cols)

  for (j in 1:num_cols) {
    mean[j] <- mean(data[, j])
  }
  return(mean)
}

# --- Parameters in nth iteration ----------------------------------------------
mean_n <- function(lambda_0, mu_0, sample_size, sample_mean) {
  # Calculates mu_n in nth iteration for a given cluster for mean posterior
  # lambda_0: prior for dividing variance for Normal
  # mu_0: int; prior of mean for Normal
  # sample_size: number of elements of current cluster
  # sample_mean: mean of current cluster
  mu_n <- ((lambda_0 * mu_0 + sample_size * sample_mean)
  / (lambda_0 + sample_size)
  )
}

S_n <- function(data, sample_mean, sample_size, num_cols) {
  # Calculates S_n for variance posterior

  # Data: data relevant to current cluster
  # Sample_mean: mean of data
  # sample_size: number of elements of data
  # num_cols: the dimensionality of data (and hence of S_n)

  # Convert data to matrix form for matrix multiplication in later steps
  data <- as.matrix(data)
  sample_covariance <- matrix(0, nrow = num_cols, ncol = num_cols)
  if (sample_size > 0) {
    for (index in 1:sample_size) {
      sample_covariance <- (sample_covariance
      + ((data[index, ] - sample_mean)
        %*% t(data[index, ] - sample_mean)
        )
      )
    }
  }
  return(sample_covariance)
}

scale_n <- function(scale_0,
                    mu_0,
                    lambda_0,
                    sample_covariance,
                    sample_size,
                    sample_mean) {

  # Calculates scale in nth iteration for a given cluster for variance posterior
  # scale_0: int; prior for scale of inverse Wishart
  # mu_0: int; prior of mean for Normal
  # lambda_0: prior for dividing variance for Normal
  # sample_covariance: output of S_n for given cluster
  # sample_size: number of elements of data
  # sample_mean: mean of data
  scale_n <- (scale_0
  + sample_covariance
    + ((lambda_0 * sample_size) / (lambda_0 + sample_size))
    * (sample_mean - mu_0) %*% t(sample_mean - mu_0)
  )
}

# --- Posterior distributions --------------------------------------------------

# Note that both mean_posterior and variance_posterior assume they receive data
# from one given cluster

# data <- cluster_data

mean_posterior <- function(mu_0, variance, lambda_0, data) {
  # Samples from a normal distribution for the mean for current variance

  # mu_0: mean prior mean
  # variance: current estimate of variance
  # lambda_0: scaling factor on variance from prior
  # data: dataframe or vector of points
  if (!(is.data.frame(data) | is.matrix(data))) {
    stop("Data is of the wrong type. ")
  }
  # if (inherits(data, "data.frame")) {
  sample_size <- nrow(data)
  num_cols <- ncol(data) ## PDWK - Added definition of num_cols

  if (sample_size) {
    sample_mean <- colMeans(data) # data_frame_mean(data)
  } else {
    sample_mean <- rep(0, num_cols)
  }
  # } else {
  #   sample_size <- length(data)
  #   if (sample_size) {
  #     sample_mean <- mean(data)
  #   } else {
  #     sample_mean <- 0
  #   }
  # }

  lambda_n <- lambda_0 + sample_size
  mu_n <- mean_n(lambda_0, mu_0, sample_size, sample_mean)

  variance_n <- variance / lambda_n

  # Take a single sample from the posterior
  mu <- mvrnorm(n = 1, mu_n, variance_n)
}

# data <- cluster_data

variance_posterior <- function(df_0, scale_0, lambda_0, mu_0, data) {
  # Samples from a inverse Wishart distribution for variance

  # df_0: prior estimate of df for prior distribution of variance
  # scale_0: prior estimate of scale for prior distribution of variance
  # lambda_0:
  # mu_0: prior estimate of mean for prior distribution of mean
  # data: dataframe or vector of points
  if (!(is.data.frame(data) | is.matrix(data))) {
    stop("Data is of the wrong type. ")
  }
  # if (is.data.frame(data) | is.matrix(data)) {
  sample_size <- nrow(data)

  num_cols <- ncol(data)

  if (sample_size > 0) {
    sample_mean <- colMeans(data)
  } else {
    sample_mean <- rep(0, num_cols)
  }
  # }

  # else {
  #   sample_size <- length(data)
  #   num_cols <- 1
  #   if (sample_size > 0) {
  #     sample_mean <- mean(data)
  #   } else {
  #     sample_mean <- 0
  #   }
  # }

  # Convert data to matrix form for matrix multiplication in later steps
  data <- as.matrix(data)

  # Calculate the component parts of the new parameters
  
  # sample_covariance <- var(data) * (sample_size - 1) # alternative
  
  # print("Components")
  # print(data)
  # print(sample_mean)
  # print(sample_size)
  # print(num_cols)
  
  
  sample_covariance <- S_n(data, sample_mean, sample_size, num_cols)
  
  # print("Hit covariance sample")
  
  # alt_cov <- (var(data) * (sample_size - 1))
  # 
  # if (any(is.na(alt_cov))) {
  #   alt_cov <- matrix(0, ncol = num_cols, nrow = num_cols)
  # }
  # 
  # # if(sample_covariance != (var(data) * (sample_size - 1))){
  # if (!isTRUE(all.equal(sample_covariance, alt_cov, tolerance = 0.001, check.attributes = F))) {
  #   print(sample_covariance)
  #   print((var(data) * (nrow(data) - 1)))
  #   stop("Problem in S_n function")
  # }

  # The effective values of the lambda, df and scale given the data
  lambda_n <- lambda_0 + sample_size
  df_n <- df_0 + sample_size

  # print(sample_covariance)
  
  scale_n_value <- scale_n(
    scale_0,
    mu_0,
    lambda_0,
    sample_covariance,
    sample_size,
    sample_mean
  )

  # print("alcu scale n")
  
  # alt_scale <- (scale_0
  # + sample_covariance
  #   + ((lambda_0 * sample_size) / (lambda_0 + sample_size))
  #   * ((sample_mean - mu_0) %*% t(sample_mean - mu_0))
  # )
  # 
  # if (any(is.na(alt_scale))) {
  #   alt_scale <- scale_0
  # }
  # 
  # if (!isTRUE(all.equal(scale_n_value, alt_scale, tolerance = 0.001, check.attributes = F))) {
  #   print(scale_n_value)
  #   print(alt_scale)
  #   stop("problem in scale_n")
  # }

  # Draw the inverse of variance as a single sample from a Wishart distribution
  # We invert the scale as the first step for moving to an inverse Wishart
  # inverse_variance <- rWishart(1, df_n, solve(scale_n_value)) # solve() inverts
  #
  # # For some reason this produces a 3D object, we want 2D I think
  # inverse_variance <- matrix(
  #   inverse_variance,
  #   dim(inverse_variance)[1],
  #   dim(inverse_variance)[2]
  # )
  #
  # # Solve for the covariance matrix
  # variance <- solve(inverse_variance)

  # Alternatively using MCMCpack
  # print(scale_n_value)
  # print(solve(scale_n_value))
  variance <- riwish(df_n, scale_n_value)
  # return(variance)
}

# This function breaks from the above in it handles all classes simultaneously
class_weight_posterior <- function(concentration_0, class_labels, k) {
  # Sample class weights from the posterior

  # concentration_0: prior estimate of concentration vector
  # class_labels: vector of numbers representing the classes of the
  # corresponding points
  # k: the number of classes
  # class_count <- rep(0, k)
  # class_weight <- rep(0, k)

  class_weight <- rep(0, k)

  # Need count of members of each class to update concentration parameter
  for (i in 1:k) {
    if (any(class_labels == i)) {
      class_count <- sum(class_labels == i)
    } else {
      class_count <- 0
    }

    # This is the parameter of the posterior
    concentration <- concentration_0 + class_count

    class_weight[i] <- rgamma(1, concentration)
  }

  class_weight <- class_weight / sum(class_weight)
}

sample_class <- function(point, data, k, class_weights, class_labels,
                         q = NULL,
                         mu = NULL,
                         variance = NULL) {

  # Samples point's class based on current parameter estimates for each class
  # Point: a data point corresponding to an entry in data
  # data: dataframe of numerical variables omitting point
  # k: the number of classes
  # q: list containing the various bits and pieces as Yee Whye demands
  # (default is NULL)
  # mu: vector of current estimates for mean of each cluster
  # variance: vector of current estimates for variance of each cluster

  if (
    !(
      (
        is.null(mu) & is.null(variance) & !(is.null(q))
      )
      | (!(is.null(mu)) & !(is.null(variance)) & is.null(q))
    )
  ) {
    stop("Either mu and variance should be declared or else q. ")
  }

  prob <- rep(0, k)
  for (i in 1:k) {
    # Use logs for nicer numbers
    curr_weight <- log(class_weights[i])
    sample_size <- sum(class_labels == i)

    if (!is.null(mu)) {
      # Exponent in likelihood function
      exponent <- -1 / 2 * ((as.matrix(point - mu[[i]]))
      %*% solve(variance[[i]])
        %*% t(as.matrix(point - mu[[i]]))
      )

      log_likelihood <- -0.5 * log((det(variance[[i]]))) + exponent

      # Weighted log-likelihood for this class
      # print(curr_weight -(sample_size* 0.5) * log((det(variance[[i]]))) + exponent)
      prob[i] <- curr_weight + log_likelihood
    } else {
      prob[i] <- curr_weight + logpredictive(q[[i]], point)
    }
  }

  # Remove logging
  prob <- exp(prob - max(prob))

  # Normalise probabilities
  prob <- prob / sum(prob)

  # Generate a random number for class allocation
  u <- runif(1)

  # Assign to class based on cumulative probabilities
  pred <- 1 + sum(u > cumsum(prob))
}

# --- MCMC analysis ------------------------------------------------------------

# Compare similarity to other points based on assignment over iterations
point_similarity <- function(cluster_record) {
  # How frequently are any two points in the same cluster across recorded iterations

  # cluster_record: (num_samples x num_iterations) matrix of ints representing
  # the cluster each point is assigned to in a given iteration

  # Record the number of points
  sample_size <- nrow(cluster_record)
  num_iter <- ncol(cluster_record)

  # Initialise the similarity matrix
  similarity_mat <- matrix(0, nrow = sample_size, ncol = sample_size)
  diag(similarity_mat) <- rep(1, sample_size)

  # symmetric matrix so iterate over i in [1, n - 1] and j in [i + 1, n]
  for (point in 1:(sample_size - 1)) {
    for (comparison_point in (point + 1):sample_size) {

      # Divide by num_iter to normalise
      similarity <- sum(
        cluster_record[point, ] == cluster_record[comparison_point, ]
      ) / num_iter

      # Assign value to both points in the matrix (due to symmetry)
      similarity_mat[point, comparison_point] <- similarity
      similarity_mat[comparison_point, point] <- similarity
    }
  }
  return(similarity_mat)
}


entropy <- function(class_weights) {
  # Measure of convergence, entropy from Information theory
  # Class weights: vector of class_weights (numbers from unit interval)
  entropy_components <- class_weights * log(class_weights)
  entropy_components[is.nan(entropy_components)] <- 0
  entropy <- sum(entropy_components)
}

entropy_window <- function(entropy_vec,
                           start = 1,
                           window_length = 25,
                           mean_tolerance = 0.001,
                           sd_tolerance = 0.001) {
  # Find the point at which entropy stabilises in the iterations

  # entropy_vec: vector of numbers corresponding to entropy of each iteration
  # start: integer instructing which iteration to start with (default is 1)
  # window_length: number of iterations to consider when considering convergence
  # (default is 25)
  # mean_tolerance: how close the mean of the two windows must be to be
  # considered converged (default is 0.001)
  # sd_tolerance: how close the sd of the two windows must be to be
  # considered converged (default is 0.001)

  n <- length(entropy_vec)

  search_range <- seq(
    from = start,
    to = n - window_length,
    by = window_length
  )

  for (i in search_range) {

    # Create two windows looking forward from the current iteration and compare
    # their means and standard deviations
    win_1 <- entropy_vec[i:(i + window_length - 1)]
    win_2 <- entropy_vec[(i + window_length)
    :min((i + 2 * window_length - 1), n)]

    mean_1 <- mean(win_1)
    mean_2 <- mean(win_2)

    sd_1 <- sd(win_1)
    sd_2 <- sd(win_2)

    # If the differences are less than the predefined tolerances, return this
    # iteration as the point to burn up to
    if ((abs(mean_1 - mean_2) < mean_tolerance)
    + (abs(sd_1 - sd_2) < sd_tolerance)
    ) {
      return(i)
    }
  }
}

# This function needs tidying
postior_sense_check <- function(data, class_labels, k, scale_0, mu_0, lambda_0, df_0,
                                num_points = 1000) {
  # Returns plots comparing the distribution actually being sampled with the
  # derived distribution for both mean and variance
  # Data: the data used in Gibbs sampler
  # k: int; number of clusters
  # scale_0: int; prior for scale of inverse Wishart
  # mu_0: int; prior of mean for Normal
  # lambda_0: prior for dividing variance for Normal
  # df_0: prior for degrees of freedom for inverse Wishart
  # num_points: int; number of points to draw from distributions

  # Declare a bunch of empty variables
  sampled_mean <- list()
  sampled_variance <- list()
  num_cols <- 1
  scale_n_value <- list()
  df_n <- rep(0, k)
  mu_n <- rep(0, k)
  lambda_n <- rep(0, k)
  posterior_mean <- list()
  posterior_variance <- list()

  # plots is the output, a plot of the sampled vs predicted means and variances
  # for each cluster
  plots <- list()
  plots$mean <- list()
  plots$variance <- list()

  # data <- as.data.frame(data[, 1])
  # scale_0 <-   as.matrix(scale_0[1, 1])
  # mu_0 <- mu_0[1]
  #
  # num_points <- 1000

  # Iterate over clusters
  for (j in 1:k) {

    # Declare the current cluster sample variables
    sampled_mean[[j]] <- rep(0, num_points)
    sampled_variance[[j]] <- matrix(0, ncol = 1, nrow = 1)
    posterior_variance[[j]] <- rep(0, num_points)

    # The various parameters and variables required for the sampling
    cluster_data <- data[class_labels == j, ] # dim(cluster_data)

    if ((!is.null(length(cluster_data))) & !(length(cluster_data) == 0)) {
      sample_mean <- mean(cluster_data)
      sample_size <- length(cluster_data)
    } else {
      sample_mean <- 0
      sample_size <- 0
    }

    sample_covariance <- S_n(cluster_data, sample_mean, sample_size, num_cols)

    scale_n_value[[j]] <- scale_n(
      scale_0,
      mu_0,
      lambda_0,
      sample_covariance,
      sample_size,
      sample_mean
    )

    df_n[j] <- df_0 + sample_size
    lambda_n[j] <- lambda_0 + sample_size

    mu_n[j] <- mean_n(lambda_0, mu_0, sample_size, sample_mean)

    # sample from the calculated posterior for the mean
    posterior_mean[[j]] <- rmvt(
      n = num_points,
      df = df_n[[j]] - d + 1,
      mu = c(mu_n[[j]]),
      S = scale_n_value[[j]] / (lambda_n[[j]] * (df_n[[j]] - d + 1))
    )


    # Generate the desired number of values
    for (i in 1:num_points) {

      # The sampled_param are the values actually being sampled by our programme
      sampled_variance[[j]][i] <- variance_posterior(
        df_0,
        scale_0,
        lambda_0,
        mu_0,
        as.matrix(cluster_data)
      )[1, 1]

      sampled_mean[[j]][i] <- mean_posterior(
        mu_0,
        as.matrix(variance[[j]][1, 1]),
        lambda_0,
        as.matrix(cluster_data)
      )

      # The distribution derived based on conjugacy
      # i.e. what we should be sampling from
      inverse_variance <- rWishart(1, df_n[[j]], solve(scale_n_value[[j]]))

      # For some reason this produces a 3D object, we want 2D I think
      inverse_variance <- matrix(
        inverse_variance,
        dim(inverse_variance)[1],
        dim(inverse_variance)[2]
      )

      # Solve for the covariance matrix
      posterior_variance[[j]][i] <- solve(inverse_variance)[1, 1]
    }

    # PLotting data for the mean values
    mu_data <- data.frame(Sample = sampled_mean[[j]], Actual = posterior_mean[[j]])

    # Mean plot for current cluster
    plots$mean[[j]] <- ggplot(data = mu_data) +
      geom_histogram(aes(x = Sample, y = ..count.. / max(..count..)),
        colour = "black",
        fill = "steelblue2"
      ) +
      geom_density(aes(x = Actual, y = ..scaled..),
        colour = "red",
        size = 0.8
      ) +
      labs(
        title = "MEAN: Comparison sampled vs predicted posteriors",
        subtitle = paste("Cluster", j),
        # caption = "",
        x = "Value",
        y = "Density"
      ) +
      NULL

    # The data for plotting the variance distributions
    var_data <- data.frame(
      Sample = sampled_variance[[j]],
      Actual = posterior_variance[[j]]
    )

    plots$variance[[j]] <- ggplot(data = var_data) +
      geom_histogram(aes(x = Sample, y = ..count.. / max(..count..)),
        colour = "black",
        fill = "steelblue2"
      ) +
      geom_density(aes(x = Actual, y = ..scaled..),
        colour = "red",
        size = 0.8
      ) +
      labs(
        title = "VARIANCE: Comparison sampled vs predicted posteriors",
        subtitle = paste("Cluster", j),
        # caption = "",
        x = "Value",
        y = "Density"
      ) +
      NULL
  }
  return(plots)
}

# --- Gibbs sampling -----------------------------------------------------------

empirical_bayes_initialise <- function(data, mu_0, df_0, scale_0, N, k, d) {
  # Creates priors for the mean, degrees of freedom and scale parameters if not 
  # set
  
  # data: data being analysed
  # mu_0: d-vector; if NULL defaults to mean of data
  # df_0: int; if NULL defaults to d + 2
  # scale_0: inveserse covariance matrix; if NULL defaults to a diagonal matrix
  # N: number of entries in data
  # k: number of classes / clusters
  # d: number of columns in data
  parameters <- list()
  if (is.null(mu_0)) {
    mu_0 <- colMeans(data)
  }

  if (is.null(df_0)) {
    df_0 <- d + 2
  }

  if (is.null(scale_0)) {
    # I think the below line does not work
    scale_0 <- diag(colSums((data - colMeans(data))^2) / N) / (k^(1 / d))
    if(any(is.na(scale_0))){
      scale_0 <- diag(d) / (k^(1 / d))
    }
  }
  parameters$mu_0 <- mu_0
  parameters$df_0 <- df_0
  parameters$scale_0 <- scale_0
  
  return(parameters)
}

gibbs_sampling <- function(data, k, class_labels,
                           d = NULL,
                           N = NULL,
                           num_iter = NULL,
                           burn = NULL,
                           mu_0 = NULL,
                           df_0 = NULL,
                           scale_0 = NULL,
                           lambda_0 = 0.01,
                           concentration_0 = 0.1) {
  # Carries out gibbs sampling of data and returns a similarity matrix for points
  
  # data: data being analysed
  # k: number of classes / clusters
  # class_labels: vector of ints representing the initial cluster of the 
  # corresponding point in data
  # num_iter: int; number of iterations to sample over
  # burn: int; number of iterations to record after
  # mu_0: d-vector; if NULL defaults to mean of data
  # df_0: int; if NULL defaults to d + 2
  # scale_0: inveserse covariance matrix; if NULL defaults to a diagonal matrix
  # lambda_0: number; prior of shrinkage for mean distribution
  # concentration_0: prior for dirichlet distribution of class weights
  
  if(is.null(d)){
    d <- ncol(data)
  }
  
  if( is.null(N)){
    N <- nrow(data)
  }

  # Empirical Bayes
  # if (is.null(mu_0)) {
  #   mu_0 <- colMeans(data)
  # }
  # 
  # if (is.null(df_0)) {
  #   df_0 <- d + 2
  # }
  # 
  # if (is.null(scale_0)) {
  #   scale_0 <- diag(colSums((data - colMeans(data))^2) / N) / (k^(1 / d))
  # }

  # Alternatively, export this to a seperate function
  if(is.null(num_iter)){
    num_iter <- (d^2) * 400
  }
  
  if(is.null(burn)){
    burn <- num_iter / 10
  }
  
  parameters_0 <- empirical_bayes_initialise(data, mu_0, df_0, scale_0, N, k, d)

  mu_0 <- parameters_0$mu_0
  df_0 <- parameters_0$df_0
  scale_0 <- parameters_0$scale_0

  concentration_0 <- rep(concentration_0, k)

  variance <- list()
  mu <- list()

  record <- matrix(0, nrow = N, ncol = num_iter - burn)
  entropy_cw <- rep(0, num_iter) # for comparison with across iterations

  # print("Enter loop")
  
  for (i in 1:num_iter) {
    class_weights <- class_weight_posterior(concentration_0, class_labels, k)
    entropy_cw[i] <- entropy(class_weights)
    for (j in 1:k) {
      cluster_data <- as.data.frame(data[class_labels == j, ])

      # print("Reached variance")
      
      # print(scale_0)
      # print(df_0)
      # print(mu_0)
      # print(lambda_0)
      # print(nrow(cluster_data))

      variance[[j]] <- variance_posterior(
        df_0,
        scale_0,
        lambda_0,
        mu_0,
        cluster_data
      )
      
      # print("Calculated")

      mu[[j]] <- mean_posterior(mu_0, variance[[j]], lambda_0, cluster_data)
    }

    for (index in 1:N) {
      point <- data[index, ]

      class_labels[index] <- sample_class(point, data, k, class_weights, class_labels,
        mu = mu,
        variance = variance
      )
    }
    if (i > burn) {
      record[, i - burn] <- t(class_labels)
    }
  }

  sim <- point_similarity(record)
}

# === Demo =====================================================================

plotting <- FALSE
set.seed(5)

# d <- 3
# N <- d * 50
# k <- 2
# num_iter <- (d^2) * 400
# burn <- num_iter / 10
# 
# # data <- mydatatrain
# data <- matrix(nrow = N, ncol = d)
# for (var_index in 1:d) {
#   data[, var_index] <- c(-2 + 0.1 * rnorm(N / 2), 2 + 0.1 * rnorm(N / 2))
# }
# 
# data <- as.data.frame(data)
# 
# mu_0 <- rep(0, d)
# df_0 <- d + 2
# scale_0 <- diag(d) / (k^(1 / d)) # diag(d) # matrix(1)
# alpha_0 <- 1
# lambda_0 <- 0.01
# concentration_0 <- rep(0.1, k)
# 
# variance <- list()
# mu <- list()
# 
# # class_labels <- c(rep(1, N/2), rep(2, N/2))
# class_labels <- sample(c(1, 2), N, replace = T)
# class_labels_0 <- class_labels
# 
# record <- matrix(0, nrow = N, ncol = num_iter - burn)
# entropy_cw <- rep(0, num_iter)

# --- Gibbs sampling -----------------------------------------------------------

# sim <- gibbs_sampling(data, k, class_labels, num_iter = num_iter)

# for (qwe in 1:num_iter) {
#   class_weights <- class_weight_posterior(concentration_0, class_labels, k)
#   entropy_cw[qwe] <- entropy(class_weights)
#   for (j in 1:k) {
#     cluster_data <- as.data.frame(data[class_labels == j, ])
# 
#     variance[[j]] <- variance_posterior(
#       df_0,
#       scale_0,
#       lambda_0,
#       mu_0,
#       cluster_data
#     )
# 
#     alt_var <- var(cluster_data) # + diag(d)
# 
#     # print("Cluster covariance:")
#     # print(alt_var)
# 
#     mu[[j]] <- mean_posterior(mu_0, variance[[j]], lambda_0, cluster_data)
# 
#     # print(mu[[j]])
#     # print("Cluster mean:")
#     # print(colMeans(cluster_data))
#   }
# 
#   for (i in 1:N) {
#     point <- data[i, ]
# 
#     class_labels[i] <- sample_class(point, data, k, class_weights, class_labels,
#       mu = mu,
#       variance = variance
#     )
#     # print(glue("Class label of individual {i} in iteration {j}:\n{class}",
#     #            i = i,
#     #            j = qwe,
#     #            class = class_labels[i]))
#   }
#   if (qwe > burn) {
#     record[, qwe - burn] <- t(class_labels)
#   }
# }
# 
# sim <- point_similarity(record)

# --- Plotting -----------------------------------------------------------------
if (plotting) {
  pheatmap(sim) # similarity
  pheatmap(1 - sim) # dissimilarity

  # The following plots only work in 1D
  if (d == 1) {
    # look at the data
    hist(data[, 1]) # , breaks = N)
    plot_data <- data.frame(X = data[, 1], Index = 1:N, Class = class_labels)
    ggplot(plot_data, aes(x = X, y = Index, colour = Class)) + geom_point()


    entropy_data <- data.frame(Index = 1:num_iter, Entropy = entropy_cw)

    burn <- entropy_window(entropy_cw,
      window_length = min(25, num_iter / 5),
      mean_tolerance = 0.01,
      sd_tolerance = 0.01
    )

    burn <- ifelse(is.null(burn), 1, burn)

    entropy_plot <- ggplot(data = entropy_data, mapping = aes(x = Index, y = Entropy)) +
      geom_point() +
      geom_vline(mapping = aes(xintercept = burn, colour = "Burn"), lty = 2) +
      # geom_smooth(se = F) +
      ggtitle("Entropy over iterations including recommended burn") +
      xlab("Iteration") + ylab("Entropy") +
      scale_color_manual(name = "", values = c(Burn = "red")) +
      NULL

    x <- entropy(class_weights)

    entropy_plot
    burn

    p <- postior_sense_check(as.data.frame(data[, 1]),
      class_labels,
      k,
      as.matrix(scale_0[1, 1]),
      mu_0[1],
      lambda_0,
      df_0,
      num_points = 1000
    )
    p
  }
}
# --- Heatmapping --------------------------------------------------------------
if(F){
# Trying to add row annotation to pheatmap
dissim <- 1 - sim

# Require names to associate data in annotation columns with original data
colnames(dissim) <- paste0("Test", 1:N)
rownames(dissim) <- paste0("Gene", 1:N)

# Example input for annotation_col in pheatmap
annotation_col <- data.frame(CellType = rep(c("CT1", "CT2"), N / 2))
rownames(annotation_col) <- paste("Test", 1:N, sep = "")

label_names <- paste("Exp", 1:length(unique(class_labels)), sep = "")

# Example input for annotation_row in pheatmap
annotation_row <- data.frame(Class_label = factor(class_labels,
  labels = label_names
))

rownames(annotation_row) <- rownames(dissim)

# Include some NAs
NA_rows <- sample(1:N, N / 5)
annotation_row[NA_rows, 1] <- NA

# This adds annotation as a row - i.e. a comment on the columns
pheatmap(dissim,
  annotation_col = annotation_col,
  annotation_row = annotation_row
)


# Colour scheme
col_pal <- RColorBrewer::brewer.pal(9, "Blues")
col_pal_alt <- ghibli_palette("MononokeLight", 21, type = "continuous")

ann_colors <- list(
  # Time = c("white", "firebrick"),
  Class_label = c(Exp1 = "#C7E9C0", Exp2 = "#00441B"),
  CellType = c(CT1 = "#7570B3", CT2 = "#E7298A")
)

# More full heatmap
pheatmap(dissim,
  annotation_col = annotation_col,
  annotation_row = annotation_row,
  annotation_colors = ann_colors,
  main = "Exmample Heatmap",
  cluster_row = T, # adds hierarchical clustering across rows
  cluster_cols = T, # adds hierarchical clustering across cols
  color = col_pal, # col_pal_alt,
  fontsize = 6.5,
  fontsize_row = 6,
  fontsize_col = 6,
  gaps_col = 50 # only has an effect when cluster_cols = F
)
}
# === Olly =====================================================================
# if(FALSE){



# pRoloc::setStockcol(paste0(pRoloc::getStockcol(), 90)) ## see through colours

data("HEK293T2011") # Human Embroyonic Kidney dataset

# head(exprs(HEK293T2011)) # proteins as rows, MS channels as columns
# 
# head(fData(HEK293T2011)) # qualitative protein information
# pData(HEK293T2011) # information about Density-gradient fractions
# getMarkerClasses(HEK293T2011) # we have 12 classes
# 
# head(fData(HEK293T2011)$markers) # labels and unknowns
# markerMSnSet(HEK293T2011) # This creates a new MSnSet with just the markers
# 
# # Visualisation
# plot2D(object = HEK293T2011, fcol = "markers", method = "PCA") # pca
# plot2D(object = HEK293T2011, fcol = "markers", method = "kpca") # kernal pca
# plot2D(object = HEK293T2011, fcol = "markers", method = "t-SNE") # t-SNE

# this function is hidden but is useful for creating a data frame with just the
# expression data and labels. train = FALSE gives unknowns
mydatatrain <- pRoloc:::subsetAsDataFrame(
  object = HEK293T2011,
  fcol = "markers",
  train = FALSE
)

# train is true gives labelled data.
mydatalabels <- pRoloc:::subsetAsDataFrame(
  object = HEK293T2011,
  fcol = "markers",
  train = TRUE
) %>%
  dplyr::sample_n(100)

# mydatalabels <- mydatalabels[, c(sample(1:(ncol(mydatalabels) - 1), 5), ncol(mydatalabels))]

class_labels <- data.frame(Class = mydatalabels$markers)
rownmaes(class_labels) <- rownames(mydatalabels)
num_data <- mydatalabels %>%
  dplyr::select(- markers)

k <- length(unique(class_labels$Class))
N <- nrow(num_data)
d <- ncol(num_data)

class_labels_key <- data.frame(Class = unique(mydatalabels$markers)) #, Class_num = 1:k)
class_labels_key %<>%
  arrange(Class) %>%
  dplyr::mutate(Class_key = as.numeric(Class))

class_labels %<>%
  mutate(Class_ind = as.numeric(mydatalabels$markers))

class_labels_0 <- sample(1:k, N, replace = T)


sim <- gibbs_sampling(num_data, k, class_labels_0,
                      N = N, 
                      d = d,
                      num_iter = 1000
                      )

# The auxiliary dataset of primary interest is the Gene Ontology Cellular
# Compartment namespace. For convenience the dataset has been put in the same
# format as the MS data.
# data("HEK293T2011goCC") # get goCC data
# head(exprs(HEK293T2011goCC)[,1:10]) # let looks at it

# }

pheatmap(sim) # similarity
pheatmap(1 - sim) # dissimilarity

dissim <- 1 - sim

# Require names to associate data in annotation columns with original data
colnames(dissim) <- rownames(num_data) #paste0("Test", 1:N)
rownames(dissim) <- rownames(num_data) #paste0("Gene", 1:N)

# Example input for annotation_col in pheatmap
annotation_col <- data.frame(CellType = rep(c("CT1", "CT2"), N / 2))
rownames(annotation_col) <- paste("Test", 1:N, sep = "")

label_names <- paste("Exp", 1:length(unique(class_labels)), sep = "")

# Example input for annotation_row in pheatmap
annotation_row <- class_labels %>% dplyr::select(Class) # data.frame(Class_label = factor(class_labels,
                                                  # labels = label_names
# ))

rownames(annotation_row) <- rownames(dissim)

# Include some NAs
# NA_rows <- sample(1:N, N / 5)
# annotation_row[NA_rows, 1] <- NA

# This adds annotation as a row - i.e. a comment on the columns
pheatmap(dissim,
         # annotation_col = annotation_col,
         annotation_row = annotation_row
)


# Colour scheme
col_pal <- RColorBrewer::brewer.pal(9, "Blues")
col_pal_alt <- ghibli_palette("MononokeLight", 21, type = "continuous")

# ann_colors <- list(
#   # Time = c("white", "firebrick"),
#   Class = c(Chromatin associated = "#C7E9C0", Cytosol = "#00441B")
# )

# newCols <- ghibli_palette("MononokeLight", length(unique(class_labels$Class)), type = "continuous")
newCols <- colorRampPalette(grDevices::rainbow(length(unique(class_labels$Class))))

# mycolors <- newCols # if using ghibli
mycolors <- newCols(length(unique(class_labels$Class)))
names(mycolors) <- unique(class_labels$Class)
mycolors <- list(Class = mycolors)

# length(unique(class_labels$Class)

# More full heatmap
pheatmap(dissim,
         # annotation_col = annotation_col,
         annotation_row = annotation_row,
         annotation_colors = mycolors,
         main = "Olly 100",
         cluster_row = T, # adds hierarchical clustering across rows
         cluster_cols = T, # adds hierarchical clustering across cols
         color = col_pal, # col_pal_alt,
         fontsize = 10,
         fontsize_row = 6,
         fontsize_col = 6,
         gaps_col = 50 # only has an effect when cluster_cols = F
)
