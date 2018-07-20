#!/usr/bin/env Rscript

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
  } else {
    sample_size <- length(data)
  }

  # Calculate the various components of the posterior parameters
  sample_mean <- mean(data)

  lambda_n <- lambda_0 + sample_size

  mu_n <- ((lambda_0 * mu_0 + sample_size * sample_mean)
  / (lambda_n)
  )

  variance_n <- variance / lambda_n

  # Take a single sample from the posterior
  mu <- rnorm(1, mu_n, variance_n)
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
    sample_mean <- data_frame_mean(data)
  } else {
    sample_size <- length(data)
    sample_mean <- mean(data)
  }

  # Convert data to matrix form for matrix multiplication in later steps
  data <- as.matrix(data)

  # Calculate the component parts of the new parameters
  # Pretty sure this is wrong
  sample_covariance <- sum(
    (data - sample_mean)
    %*% t(data - sample_mean)
  )
  
  # print(glue("Sample covariance pre-sum:\n{x}\n", x = (data - sample_mean) %*% t(data - sample_mean)))
  # print((data - sample_mean) %*% t(data - sample_mean))
  # print(sum((data - sample_mean) %*% t(data - sample_mean)))
  # print(rowSums((data - sample_mean) %*% t(data - sample_mean)))
  # print("")
  # print("")
  

  # The effective values of the lambda, df and scale given the data
  lambda_n <- lambda_0 + sample_size

  df_n <- df_0 + sample_size

  scale_n <- (scale_0
  + sample_covariance
    + ((lambda_0 * sample_size) / (lambda_n))
      * (sample_mean - mu_0) %*% t(sample_mean - mu_0)
  )
  
  # print(glue("Sample covariance: \n{x}\n", x = sample_covariance))
  # print(glue("Scale: \n{x}\n\n", x = scale_n))

  awk_part <- ((lambda_0 * sample_size) / (lambda_n)) * (sample_mean - mu_0) %*% t(sample_mean - mu_0)
  
  first_part <- ((lambda_0 * sample_size) / (lambda_n))
  sec_part <- (sample_mean - mu_0) %*% t(sample_mean - mu_0)
  
  # print(glue("Components of scale_n\nTotal: {total}\nFirst comp: {first}\nSecond comp: {second}\n\n",
  #            total = awk_part,
  #            first = first_part,
  #            second = sec_part)
  #       )
  
  # Draw the variance as a single sample from a Wishart distribution
  variance <- rWishart(1, df_n, scale_n)

  variance <- matrix(variance, dim(variance)[1], dim(variance)[2])
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
    # class_count[i] <- sum(class_labels == i)
    #
    # # This is the parameter of the posterior
    # concenctration <- concentration_0 + class_count[i]
    #
    # class_weight[i] <- rgamma(k, concenctration)
    #
    class_count <- sum(class_labels == i)

    # This is the parameter of the posterior
    concentration <- concentration_0 + class_count

    class_weight[i] <- rgamma(1, concentration)
  }

  class_weight <- class_weight / sum(class_weight)

  # # This is the parameter of the posterior
  # concenctration <- concentration_0 + class_count
  #
  # # Draw k (one for each class) samples from the posterior dirichlet
  # class_weight <- rdirichlet(k, concenctration)
}

sample_class <- function(point, data, class_labels, k, class_weights,
                         q = NULL,
                         mu = NULL,
                         variance = NULL) {

  # Samples point's class based on current parameter estimates for each class
  # Point: a data point corresponding to an entry in data
  # data: dataframe of numerical variables omitting point
  # class_labels: vector of numbers representing the classes of the
  # corresponding points
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
  print(pred)
  return(pred)
}

# === Demo =====================================================================

N <- 10
k <- 2
data <- c(rnorm(N / 2, -4, 1), rnorm(N / 2, 4, 1))
class_labels <- sample(c(1, 2), N, replace = T)
mu_0 <- 0
# variance0 <- matrix(5)
df_0 <- 1
scale_0 <- matrix(3)
alpha_0 <- 0.1
lambda_0 <- 1
concentration_0 <- rep(0.1, k)

variance <- list()
mu <- list()

num_iter <- 10

for (qwe in 1:num_iter) {
  print(glue("\n\nIteration: {iter}", iter = qwe))
  class_weights <- class_weight_posterior(concentration_0, class_labels, k)
  
  print(glue("Class weights: \n{x}, {y}\n", x = class_weights[1], y = class_weights[2]))
  
  for (j in 1:k) {
    cluster_data <- data[class_labels == j]
    
    variance[[j]] <- variance_posterior(df_0, scale_0, lambda_0, mu_0, cluster_data)
    mu[[j]] <- mean_posterior(mu_0, variance[[j]], lambda_0, cluster_data)
  }
  
  print(glue("Class variance: \n{x} \n{y}\n", x = variance[1], y = variance[2]))
  print(glue("Class means: \n{x}\n{y}\n\n", x = mu[1], y = mu[2]))
  
  for (i in 1:N) {
    # print(glue("Individual: {x}", x = i))
    point <- data[i]
    # print(glue("Value: {x}", x = point))
    class_labels[i] <- sample_class(point, data, class_labels, k, class_weights,
                                    mu = mu, 
                                    variance = variance
                                    )
  }
}

prob <- rep(0, k)
for (i in 1:k) {
  # Use logs for nicer numbers
  curr_weight <- log(class_weights[i])
  
  # Exponent in likelihood function
  exponent <- -1 / 2 * (t(as.matrix(point - mu[[i]])) 
                        %*% solve(variance[[i]]) 
                        %*% (as.matrix(point - mu[[i]])))
  
  # Weighted log-likelihood for this class
  prob[i] <- curr_weight - 0.5 * log((det(variance[[i]]))) + exponent
}
# Remove logging
prob <- exp(prob - max(prob))

# Normalise probabilities
prob <- prob / sum(prob)

# Generate a random number for class allocation
u <- runif(1)

# Assign to class based on cumulative probabilities
pred <- 1 + sum(u > cumsum(prob))
print(pred)
