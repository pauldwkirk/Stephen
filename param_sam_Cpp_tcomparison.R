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
## try http:// if https:// URLs are not supported
# source("https://bioconductor.org/biocLite.R")
# biocLite("pRoloc")
# biocLite("pRolocdata")
require(pRoloc)
require(pRolocdata)

# for t-SNE in plot2D
require(Rtsne) # install.packages("Rtsne", dep = T)

# for %<>%
library(magrittr)

# to access C++ files
library(Rcpp)

sourceCpp("sampleRcpp_tcomparison.cpp")

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


# --- MCMC analysis ------------------------------------------------------------

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
    scale_0 <- diag(colSums((data - mean(data))^2) / N) / (k^(1 / d))
    if (any(is.na(scale_0))) {
      scale_0 <- diag(d) / (k^(1 / d))
    }
  }
  parameters$mu_0 <- mu_0
  parameters$df_0 <- df_0
  parameters$scale_0 <- scale_0

  return(parameters)
}

gibbs_sampling <- function(data, k, class_labels, fix_vec,
                           d = NULL,
                           N = NULL,
                           num_iter = NULL,
                           burn = 0,
                           mu_0 = NULL,
                           df_0 = NULL,
                           scale_0 = NULL,
                           lambda_0 = 0.01,
                           concentration_0 = 0.1,
                           thinning = 25,
                           outlier = FALSE,
                           t_df = 4.0) {
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
  # thinning: int; record results whenever the iteration number is a multiple of
  # this after burn in

  if (is.null(d)) {
    d <- ncol(data)
  }

  if (is.null(N)) {
    N <- nrow(data)
  }


  if (is.null(num_iter)) {
    num_iter <- min((d^2) * 1000 / sqrt(N), 10000)
  }

  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }

  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  if (thinning > (num_iter - burn)) {
    if (thinning > (num_iter - burn) & thinning < 5 * (num_iter - burn)) {
      stop("Thinning factor exceeds iterations feasibly recorded. Stopping.")
    } else if (thinning > 5 * (num_iter - burn) & thinning < 10 * (num_iter - burn)) {
      stop("Thinning factor relatively large to effective iterations. Stopping algorithm.")
    } else {
      warning(paste0(
        "Thinning factor relatively large to effective iterations.",
        "\nSome samples recorded. Continuing but please check input"
      ))
    }
  }

  data <- as.matrix(data)

  # Empirical Bayes
  parameters_0 <- empirical_bayes_initialise(data, mu_0, df_0, scale_0, N, k, d)

  mu_0 <- parameters_0$mu_0
  df_0 <- parameters_0$df_0
  scale_0 <- parameters_0$scale_0

  if (is.null(concentration_0)) {
    concentration_0 <- rep(0.1, k)
  } else if (length(concentration_0) < k) {
    print(paste0(
      "Creating vector of ", k + outlier, " repetitions of ", concentration_0,
      " for concentration prior."
    ))
    concentration_0 <- rep(concentration_0, k + outlier)
  }

  sim <- point_comparison(
    num_iter,
    concentration_0,
    scale_0,
    class_labels,
    fix_vec,
    mu_0,
    lambda_0,
    data,
    df_0,
    k,
    burn,
    thinning,
    outlier,
    t_df
  )
}

# === Wrapper ==================================================================

mcmc_out <- function(MS_object,
                     class_labels_0 = NULL,
                     mu_0 = NULL,
                     df_0 = NULL,
                     scale_0 = NULL,
                     lambda_0 = 0.01,
                     concentration_0 = 0.1,
                     train = NULL,
                     num_iter = NULL,
                     burn = NULL,
                     thinning = 25,
                     heat_plot = TRUE,
                     main = "heatmap_for_similarity",
                     cluster_row = T,
                     cluster_cols = T,
                     fontsize = 10,
                     fontsize_row = 6,
                     fontsize_col = 6,
                     gaps_col = 50,
                     entropy_plot = TRUE,
                     window_length = min(25, num_iter / 5),
                     mean_tolerance = 0.0005,
                     sd_tolerance = 0.0005,
                     outlier = FALSE,
                     t_df = 4.0,
                     prediction_threshold = 0.6) {
  # Returns mean, variance and similarity posteriors from Gibbs sampling with
  # option of pheatmap

  # MS_object: a dataset in the format used by pRolocdata
  # class_labels_0: optional prior for classes of MS_object if NULL defaults to
  # a randomly generated set
  # mu_0: d-vector; prior for mean. If NULL defaults to mean of data
  # df_0: int; prior for degrees of freedom for inverse Wishart describing the
  # variance. If NULL, defaults to d + 2
  # scale_0: inveserse covariance matrix; prior of scale in inverse Wishart
  # describing the variance.  if NULL defaults to a diagonal matrix
  # lambda_0: number; prior of shrinkage for mean distribution
  # concentration_0: prior for dirichlet distribution of class weights
  # train: instruction to include all data (NULL), labelled data (TRUE) or
  # unlabelled data (FALSE). Default is NULL.
  # num_iter: int; number of iterations to sample over
  # burn: int; number of iterations to record after
  # thinning: int; record results whenever the iteration number is a multiple of
  # this after burn in
  # heat_plot: bool; instructs saving and printing of heatmap of similarity
  # matrix from Gibbs sampling. Default is TRUE.
  # main: string; title for heatmap, default is "heatmap_for_similarity"
  # cluster_row: bool; instructs pheatmap to cluster rows using a tree
  # cluster_cols: bool; instructs pheatmap to cluster columns using a tree
  # fontsize: int; size of font in pheatmap
  # fontsize_row: int; fontsize in rows
  # fontsize_col: int; fontsize in columns
  # gaps_col: int; not really present

  # Data with labels
  mydata_labels <- pRoloc:::subsetAsDataFrame(
    object = MS_object,
    fcol = "markers",
    train = TRUE
  )

  fixed <- rep(TRUE, nrow(mydata_labels))

  mydata_no_labels <- pRoloc:::subsetAsDataFrame(
    object = MS_object,
    fcol = "markers",
    train = FALSE
  )

  not_fixed <- rep(FALSE, nrow(mydata_no_labels))

  nk <- tabulate(fData(markerMSnSet(HEK293T2011))[, "markers"])

  mydata_no_labels$markers <- NA

  if (is.null(train)) {
    mydata <- bind_rows(mydata_labels, mydata_no_labels)
    fix_vec <- c(fixed, not_fixed)
  } else if (isTRUE(train)) {
    mydata <- mydata_labels
    fix_vec <- fixed
  } else {
    mydata <- mydata_no_labels
    fix_vec <- not_fixed
  }

  class_labels <- data.frame(Class = mydata$markers)

  classes_present <- unique(fData(markerMSnSet(HEK293T2011))[, "markers"])

  rownames(class_labels) <- rownames(mydata)

  # Numerical data of interest for clustering
  num_data <- mydata %>%
    dplyr::select(-markers)

  # Parameters
  k <- length(classes_present)
  N <- nrow(num_data)
  d <- ncol(num_data)

  # Key to transforming from int to class
  class_labels_key <- data.frame(Class = classes_present) # , Class_num = 1:k)
  class_labels_key %<>%
    arrange(Class) %>%
    dplyr::mutate(Class_key = as.numeric(Class))

  class_labels %<>%
    mutate(Class_ind = as.numeric(mydata$markers))

  # Generate class labels
  if (is.null(class_labels_0)) {
    class_weights <- nk / sum(nk)
    if (is.null(train)) {
      fixed_labels <- as.numeric(fData(markerMSnSet(HEK293T2011))[, "markers"])
      class_labels_0 <- c(fixed_labels, sample(1:k, nrow(mydata_no_labels),
        replace = T,
        prob = class_weights
      ))
    } else if (isTRUE(train)) {
      class_labels_0 <- as.numeric(fData(markerMSnSet(HEK293T2011))[, "markers"])
    } else {
      class_labels_0 <- sample(1:k, N, replace = T, prob = class_weights)
    }
  }

  if (is.null(num_iter)) {
    num_iter <- min((d^2) * 1000 / sqrt(N), 10000)
  }

  if (is.null(burn)) {
    burn <- floor(num_iter / 10)
  }

  if (burn > num_iter) {
    stop("Burn in exceeds total iterations. None will be recorded.\nStopping.")
  }

  if (thinning > (num_iter - burn)) {
    if (thinning > (num_iter - burn) & thinning < 5 * (num_iter - burn)) {
      stop("Thinning factor exceeds iterations feasibly recorded. Stopping.")
    } else if (thinning > 5 * (num_iter - burn) & thinning < 10 * (num_iter - burn)) {
      stop("Thinning factor relatively large to effective iterations. Stopping algorithm.")
    } else {
      warning(paste0(
        "Thinning factor relatively large to effective iterations.",
        "\nSome samples recorded. Continuing but please check input"
      ))
    }
  }

  gibbs <- gibbs_sampling(num_data, k, class_labels_0, fix_vec,
    d = d,
    N = N,
    num_iter = num_iter,
    burn = burn,
    mu_0 = mu_0,
    df_0 = df_0,
    scale_0 = scale_0,
    lambda_0 = lambda_0,
    concentration_0 = concentration_0,
    thinning = thinning,
    outlier = outlier,
    t_df = t_df
  )

  print("Gibbs sampling complete")
  
  
  # Create a dataframe for the predicted class
  class_allocation_table <- with(
    stack(data.frame(t(gibbs$class_record))),
    table(ind, values)
  )
  
  eff_iter <- ceiling((num_iter - burn) / thinning)
  
  # Create a column Class_key containing an integer in 1:k representing the most
  # common class allocation, and a Count column with the proportion of times the 
  # entry was allocated to said class
  predicted_classes <- data.frame(
    Class_key =
      as.numeric(colnames(class_allocation_table)
                 [apply(
                   class_allocation_table,
                   1,
                   which.max
                 )]),
    Count = apply(class_allocation_table, 1, max) / eff_iter
  )
  
  # Change the prediction to NA for any entry with a proportion below the input
  # threshold
  predicted_classes[predicted_classes$Count < prediction_threshold, ] = NA
  
  predicted_classes$Class <- class_labels_key$Class[match(
    predicted_classes$Class_key,
    class_labels_key$Class_key)]
  
  gibbs$predicted_class <- predicted_classes
  

  if (heat_plot) {

    # dissimilarity matrix
    dissim <- 1 - gibbs$similarity

    # Require names to associate data in annotation columns with original data
    colnames(dissim) <- rownames(num_data)
    rownames(dissim) <- rownames(num_data)

    # Example input for annotation_row in pheatmap
    annotation_row <- class_labels %>% dplyr::select(Class)
    rownames(annotation_row) <- rownames(dissim)

    # Colour scheme for heatmap
    col_pal <- RColorBrewer::brewer.pal(9, "Blues")

    # Annotation colours
    newCols <- colorRampPalette(grDevices::rainbow(length(classes_present)))
    mycolors <- newCols(length(classes_present))
    names(mycolors) <- classes_present
    mycolors <- list(Class = mycolors)

    # Heatmap
    if (is.null(train) | isTRUE(train)) {
      heat_map <- pheatmap(dissim,
        annotation_row = annotation_row,
        annotation_colors = mycolors,
        main = main,
        cluster_row = cluster_row,
        cluster_cols = cluster_cols,
        color = col_pal,
        fontsize = fontsize,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        gaps_col = gaps_col
      )
    } else {
      heat_map <- pheatmap(dissim,
        main = main,
        cluster_row = cluster_row,
        cluster_cols = cluster_cols,
        color = col_pal,
        fontsize = fontsize,
        fontsize_row = fontsize_row,
        fontsize_col = fontsize_col,
        gaps_col = gaps_col
      )
    }
  }
  if (entropy_plot) {
    entropy_data <- data.frame(Index = 1:num_iter, Entropy = gibbs$entropy)

    rec_burn <- entropy_window(gibbs$entropy,
      window_length = window_length,
      mean_tolerance = mean_tolerance,
      sd_tolerance = sd_tolerance
    )

    # Check if instantly ok
    rec_burn <- ifelse(is.null(rec_burn), 1, rec_burn)


    entropy_scatter <- ggplot(data = entropy_data, mapping = aes(x = Index, y = Entropy)) +
      geom_point() +
      geom_vline(mapping = aes(xintercept = rec_burn, colour = "Reccomended"), lty = 2) +
      geom_vline(mapping = aes(xintercept = burn, colour = "Implemented"), lty = 4) +
      ggtitle("Entropy over iterations including recommended and implemented burn") +
      xlab("Iteration") + ylab("Entropy") +
      scale_color_manual(name = "Burn", values = c(
        Reccomended = "red",
        Implemented = "blue"
      ))
  }
  if (heat_plot & entropy_plot) {
    return(list(
      gibbs = gibbs,
      heat_map = heat_map,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn
    ))
  }
  if (heat_plot) {
    return(list(
      gibbs = gibbs,
      heatmap = heat_map
    ))
  }
  if (entropy_plot) {
    return(list(
      gibbs = gibbs,
      entropy_plot = entropy_scatter,
      rec_burn = rec_burn
    ))
  }

  return(list(gibbs = gibbs))
}

# === Olly =====================================================================

set.seed(5)

# MS object
data("HEK293T2011") # Human Embroyonic Kidney dataset

t1 <- Sys.time()

stuff <- mcmc_out(HEK293T2011,
  num_iter = 100,
  burn = 10,
  thinning = 10,
  outlier = TRUE,
  heat_plot = FALSE
  # main = "Gene clustering by organelle"
)

t2 <- Sys.time()

t2 - t1 # how long does it take

# To plot the entropy over iterations
stuff$entropy_plot
str(stuff$gibbs$class_prob)

str(stuff$gibbs$class_record)

stuff$gibbs$predicted_class

y <- stuff$gibbs$class_record
z <- with(stack(data.frame(t(y))), table(ind, values))
predicted_classes <- data.frame(Class_key = as.numeric(colnames(z)[apply(z, 1, which.max)]),
                                Count = apply(z, 1, max))

predicted_classes$Class <- df2$B[match(df1$Class_key, df2$Class_key)]

predicted_classes %<>%
  dplyr::mutate(Class = 0)

summary(stuff$gibbs$predicted_class)
