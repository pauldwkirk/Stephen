#!/usr/bin/env Rscript

# === Libraries ================================================================

library(expm) # install.packages("expm", dep = T)
library(glue) # install.packages("glue", dep = T)
library(tidyverse) # install.packages("tidyverse", dep = T)
library(MASS)
library(pheatmap) # install.packages("pheatmap", dep = T)

# === Functions ================================================================

# Note that both mean_posterior and variance_posterior assume they receive data
# from one given cluster

mean_posterior <- function(mu_0, variance, lambda_0, data) {
  # Samples from a normal distribution for the mean for current variance

  # mu_0: mean prior mean
  # variance: current estimate of variance
  # lambda_0: scaling factor on variance from prior
  # data: dataframe or vector of points
  if (!(inherits(data, "data.frame") + is.numeric(data))) {
    stop("Data is of the wrong type. ")
  }
  if (inherits(data, "data.frame")) {
    sample_size <- nrow(data)

    if (sample_size) {
      sample_mean <- data_frame_mean(data)
    } else {
      sample_mean <- rep(0, num_cols)
    }
  } else {
    sample_size <- length(data)
  }

  # Calculate the various components of the posterior parameters
  if (sample_size) {
    sample_mean <- mean(data)
  } else {
    sample_mean <- 0
  }

  lambda_n <- lambda_0 + sample_size

  mu_n <- ((lambda_0 * mu_0 + sample_size * sample_mean)
  / (lambda_n)
  )

  variance_n <- variance / lambda_n

  # Take a single sample from the posterior
  mu <- mvrnorm(n = 1, mu_n, variance_n)
}

variance_posterior <- function(df_0, scale_0, lambda_0, mu_0, data) {
  # Samples from a inverse Wishart distribution for variance

  # df_0: prior estimate of df for prior distribution of variance
  # scale_0: prior estimate of scale for prior distribution of variance
  # lambda_0:
  # mu_0: prior estimate of mean for prior distribution of mean
  # data: dataframe or vector of points
  if (!(inherits(data, "data.frame") + is.numeric(data))) {
    stop("Data is of the wrong type. ")
  }
  if (inherits(data, "data.frame")) {
    sample_size <- nrow(data)

    num_cols <- ncol(data)

    if (sample_size > 0) {
      sample_mean <- data_frame_mean(data)
    } else {
      sample_mean <- rep(0, num_cols)
    }
  } else {
    sample_size <- length(data)
    num_cols <- 1
    if (sample_size > 0) {
      sample_mean <- mean(data)
    } else {
      sample_mean <- 0
    }
  }

  # Convert data to matrix form for matrix multiplication in later steps
  data <- as.matrix(data)

  # Calculate the component parts of the new parameters

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

  # The effective values of the lambda, df and scale given the data
  lambda_n <- lambda_0 + sample_size

  df_n <- df_0 + sample_size

  scale_n <- (scale_0
  + sample_covariance
    + ((lambda_0 * sample_size) / (lambda_n))
    * (sample_mean - mu_0) %*% t(sample_mean - mu_0)
  )

  # Draw the inverse of variance as a single sample from a Wishart distribution
  # We invert the scale as the first step for moving to an inverse Wishart
  inverse_variance <- rWishart(1, df_n, solve(scale_n)) # solve() inverts

  # For some reason this produces a 3D object, we want 2D I think
  inverse_variance <- matrix(
    inverse_variance,
    dim(inverse_variance)[1],
    dim(inverse_variance)[2]
  )

  # Solve for the covariance matrix
  variance <- solve(inverse_variance)
}

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

sample_class <- function(point, data, k, class_weights,
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

  # This if statement does not work
  # if (
  #   !(
  #     (
  #       is.null(mu) & is.null(variance) & !(is.null(variance))
  #     )
  #     | (!(is.null(mu)) & !(is.null(variance)) & is.null(variance))
  #   )
  # ) {
  #   stop("Either both mu and variance should be declared or neither. ")
  # }

  prob <- rep(0, k)
  for (i in 1:k) {
    # Use logs for nicer numbers
    curr_weight <- log(class_weights[i])

    if (!is.null(mu)) {
      # Exponent in likelihood function
      exponent <- -1 / 2 * (t(as.matrix(point - mu[[i]]))
      %*% solve(variance[[i]])
        %*% (as.matrix(point - mu[[i]]))
      )

      # Weighted log-likelihood for this class
      prob[i] <- curr_weight + log((det(variance[[i]])^(-0.5))) + exponent
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

# === Demo =====================================================================

N <- 100
k <- 2
data <- c(rnorm(N / 2, -2, 1), rnorm(N / 2, 2, 1)) # hist(data)
mu_0 <- 0
df_0 <- 1
scale_0 <- matrix(1)
alpha_0 <- 0.1
lambda_0 <- 1
concentration_0 <- rep(0.1, k)

variance <- list()
mu <- list()
class_labels <- sample(c(1, 2), N, replace = T)
class_labels_0 <- class_labels

num_iter <- 1000
burn <- 0
record <- matrix(0, nrow = N, ncol = num_iter - burn)
entropy_cw <- rep(0, num_iter)

for (qwe in 1:num_iter) {
  class_weights <- class_weight_posterior(concentration_0, class_labels, k)
  entropy_cw[qwe] <- entropy(class_weights)
  for (j in 1:k) {
    cluster_data <- data[class_labels == j]

    variance[[j]] <- variance_posterior(df_0, scale_0, lambda_0, mu_0, cluster_data)
    mu[[j]] <- mean_posterior(mu_0, variance[[j]], lambda_0, cluster_data)
  }

  for (i in 1:N) {
    point <- data[i]

    class_labels[i] <- sample_class(point, data, k, class_weights,
      mu = mu,
      variance = variance
    )
  }
  if (qwe > burn) {
    record[, qwe - burn] <- t(class_labels)
  }
}

sim <- point_similarity(record)

entropy_data <- data.frame(Index = 1:num_iter, Entropy = entropy_cw)

burn <- entropy_window(entropy_cw)

entropy_plot <- ggplot(data = entropy_data, mapping = aes(x = Index, y = Entropy)) +
  geom_point() +
  geom_vline(mapping = aes(xintercept = burn, colour = "Burn"), lty = 2) +
  # geom_smooth(se = F) +
  # geom_line() +
  ggtitle("Entropy over iterations including recommended burn") +
  xlab("Iteration") + ylab("Entropy") +
  scale_color_manual(name = "", values = c(Burn = "red")) +
  NULL

pheatmap(sim) # similarity
pheatmap(1 - sim) # dissimilarity

# look at the data
hist(data)
plot_data <- data.frame(X = data, Index = 1:N, Class = class_labels)
ggplot(plot_data, aes(x = X, y = Index, colour = Class)) + geom_point()
x <- entropy(class_weights)

entropy_plot
burn
