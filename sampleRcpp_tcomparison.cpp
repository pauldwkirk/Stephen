# include <RcppArmadillo.h>
# include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// Compares how similar two points are with regards to their clustering across 
// all iterations.
// Works for unsupervised methods (i.e. allows label flipping)
double point_similarity(arma::uword point, 
                        arma::uword comparison_point,
                        arma::umat cluster_record,
                        arma::uword num_iter) {
  double out = 0.0;

  for (arma::uword i = 0; i < num_iter; i++){
    if(cluster_record(point, i) == cluster_record(comparison_point, i)){
      out++;
    }
    
  }
  out = out / num_iter;
  return out;
}

// Constructs a similarity matrix comparing all points clustering across the 
// iterations
arma::mat similarity_mat(arma::umat cluster_record){
  arma::uword sample_size = cluster_record.n_rows;
  arma::uword num_iter = cluster_record.n_cols;
  arma::mat out(sample_size, sample_size);
  
  for (arma::uword point = 0; point < sample_size; point++){ // if not doing diagonal, restrict to sample size - 1
    for (arma::uword comparison_point = point; // + 1; 
         comparison_point < sample_size;
         comparison_point++){
      out(point, comparison_point) = point_similarity(point, 
          comparison_point,
          cluster_record,
          num_iter);
      out(comparison_point, point) = out(point, comparison_point);
    }
  }
  return out;
}

// Returns a variable involved in updating the scale, it is similar to sample 
// covariance
arma::mat S_n(arma::mat data,
              arma::vec sample_mean,
              arma::uword sample_size,
              arma::uword num_cols){
  arma::mat sample_covariance(num_cols, num_cols);
  sample_covariance.zeros();
  if(sample_size > 0){
    for(arma::uword i = 0; i < sample_size; i++){
      
      arma::vec row_i = trans(data.row(i));
      
      sample_covariance = (sample_covariance 
                             + ((row_i - sample_mean) 
                                  * arma::trans(row_i - sample_mean)
                             )
      );
      
    }
  }
  return sample_covariance;
}

// Returns a vector of the mean of a cluster of size n
arma::vec mean_n(double lambda_0,
                 arma::vec mu_0,
                 arma::uword sample_size,
                 arma::uword num_cols,
                 arma::vec sample_mean){
  arma::vec mu_n(num_cols);
  mu_n = ((lambda_0 * mu_0 + sample_size * sample_mean)
            / (lambda_0 + sample_size));
  return mu_n;
}

// Returns the matrix for the scale hyperparameter for n observations
arma::mat scale_n(arma::mat scale_0,
                  arma::vec mu_0,
                  double lambda_0,
                  arma::mat sample_covariance,
                  arma::uword sample_size,
                  arma::vec sample_mean){
  arma::mat scale_out;
  
  scale_out = (scale_0
                 + sample_covariance
                 + ((lambda_0 * sample_size) / (lambda_0 + sample_size))
                 * (sample_mean - mu_0) * arma::trans(sample_mean - mu_0)
  ); 
  return scale_out;
}

// sample a multivariate normal
arma::mat mvrnormArma(arma::uword n,
                      arma::vec mu,
                      arma::mat sigma) {
  arma::uword ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() +
    Y * arma::chol(sigma);
}

// sample from the mean posterior
arma::vec mean_posterior(arma::vec mu_0,
                         arma::mat variance,
                         double lambda_0,
                         arma::mat data){
  arma::uword ncols = data.n_cols;
  arma::vec sample_mean(ncols); sample_mean.zeros();
  arma::vec mu_out(ncols);
  arma::vec mu_n(ncols);
  arma::mat variance_n(ncols, ncols);
  
  arma::uword sample_size = data.n_rows;
  if (sample_size > 0){
    sample_mean = trans(arma::mean(data, 0));
  }
  
  // std::cout << "\nPast initial if\n";
  
  double lambda_n = lambda_0 + sample_size;
  
  mu_n = mean_n(lambda_0, mu_0, sample_size, ncols, sample_mean);
  
  // std::cout << "\nPast initial mu_n\n";
  // std::cout << "\n" << mu_n << "\n";
  
  variance_n = variance / lambda_n;
  
  // std::cout << "\n" << variance_n << "\n";
  
  arma::mat x = mvrnormArma(1, mu_n, variance_n);
  
  mu_out = arma::conv_to<arma::vec>::from(x.row(0));
  
  return mu_out;
  
}

// sample the concentration parameter for a dirichlet distribution after n 
// observations
arma::vec concentration_n(arma::vec concentration_0,
                          arma::uvec class_labels,
                          arma::uword k){
  
  arma::uword n = class_labels.n_elem;
  
  // std::cout << n ;
  
  arma::uword class_count;
  arma::vec concentration(k);
  
  for (arma::uword i = 1; i < k + 1; i++) {
    class_count = 0;
    
    for (arma::uword j = 0; j < n; j++ ) {
      if (class_labels(j) == i) {
        class_count++;
      }
    }
    
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

// update the cluster weights based on current class_labels
arma::vec class_weight_posterior(arma::vec concentration_0,
                                 arma::uvec class_labels,
                                 arma::uword k){
  arma::vec class_weight = arma::zeros<arma::vec>(k);

  // std::cout << concentration_0;
  
  arma::vec concentration(k);
  concentration = concentration_n(concentration_0,
                                  class_labels,
                                  k);
  
  
  for (arma::uword i = 1; i < k + 1; i++) {
    
    class_weight(i - 1) = Rf_rgamma(arma::as_scalar(concentration(i - 1)), 1);

  }

  // std::cout << class_weight;
  
  double total_class_weight = sum(class_weight);
  class_weight = class_weight / total_class_weight;
  return class_weight;
}


// Calculate the entropy for the current cluster weights
// [[Rcpp::export]]
double entropy(arma::vec class_weights){
  arma::uword n = class_weights.n_elem;
  arma::vec entropy_components(n);
  // std::cout << "\nDeclared\n";
  
  
  for(arma::uword i = 0; i < n; i++){
    entropy_components(i) = - class_weights(i) * log(class_weights(i));
    if (entropy_components.has_nan()){
      entropy_components(i) = 0;
    }
  }
  // std::cout << "Inter";
  double entropy_out = sum(entropy_components);
  return entropy_out;
}

// Sample the variance for after n observations
arma::mat variance_posterior(int df_0,
                             arma::mat scale_0,
                             double lambda_0,
                             arma::vec mu_0,
                             arma::mat data){
  
  arma::uword sample_size = data.n_rows, num_cols = data.n_cols;

  // std::cout << "\nCluster data:\n" << data;
  
  arma::vec sample_mean(num_cols); sample_mean.zeros();
  
  int df_n = df_0 + sample_size;
  
  // std::cout << "Reached another dec\n";
  
  arma::mat scale_n_value(num_cols, num_cols);
  arma::mat sample_covariance(num_cols, num_cols);
  arma::mat variance_out(num_cols, num_cols);

  if (sample_size > 0){
    
    sample_mean = arma::trans(arma::mean(data, 0));
    
  } else{
    sample_mean.fill(0.0);
  }
  
  // std::cout << "\nSample covariance reached\n";
  
  sample_covariance = S_n(data, sample_mean, sample_size, num_cols);

  // std::cout << "Scale_n reached\n";
  
  arma::mat samp_cov(num_cols, num_cols);
  samp_cov = (sample_size - 1) * arma::cov(data);
  
  // std::cout  << samp_cov << "\n\n";
  
  scale_n_value = scale_n(
    scale_0,
    mu_0,
    lambda_0,
    sample_covariance,
    sample_size,
    sample_mean
  );
  
  variance_out = arma::iwishrnd(scale_n_value, df_n);
  
  return variance_out;
  
} 

arma::vec sample_class(arma::vec point,
                       arma::mat data,
                       arma::uword k,
                       arma::vec class_weights,
                       arma::uvec class_labels,
                       arma::cube mu,
                       arma::cube variance,
                       bool outlier = false,
                       arma::vec global_mean = arma::zeros<arma::vec>(1),
                       arma::mat global_variance = arma::zeros<arma::mat>(1, 1),
                       double t_df = 4.0){

  double curr_weight;
  double exponent;
  double log_likelihood;
  
  
  double log_det;
  arma::vec prob_vec(k);

  arma::uvec count_probs;
  
  arma::uword d = data.n_cols;;
  
  for(arma::uword i = 1; i < k + 1; i++){
    curr_weight = log(class_weights(i - 1));

    if(outlier && i == k){

      exponent = arma::as_scalar(arma::trans(point - global_mean) 
                                   * arma::inv(global_variance)
                                   * (point - global_mean));
                                   
      log_det = arma::log_det(global_variance).real();
      
      log_likelihood = (lgamma((t_df + d)/2.0) - lgamma(t_df/2.0)
                          + d/2.0 * log(t_df * M_PI) + log_det - ((t_df + d)/2.0) *
                            log(1 + (1/t_df) * exponent));
                        
    }
    else {
      exponent = arma::as_scalar(arma::trans(point - mu.slice(i - 1)) 
                                   * arma::inv(variance.slice(i - 1))
                                   * (point - mu.slice(i - 1)));
                                   
      log_det = arma::log_det(variance.slice(i - 1)).real();
      log_likelihood = -0.5 *(log_det + exponent + d * log(2 * M_PI));
    }
    prob_vec(i - 1) = curr_weight + log_likelihood;
  } 
  prob_vec = exp(prob_vec - max(prob_vec));
  prob_vec = prob_vec / sum(prob_vec);
  
  return prob_vec;
}

// Predicts the cluster for a point given the vector of probabilities associated
// with each one
arma::uword select_class(arma::vec prob_vec){
  
  // std::cout << prob_vec << "\n\n";
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  pred = 1 + sum(u > cumsum(prob_vec));
  return pred;
}

// Returns a 2-D field of cubes for the current iterations means and variances 
// across each cluster (hence a field)
arma::field<arma::cube> mean_variance_sampling(arma::mat data,
                                               arma::uvec cluster_labels,
                                               arma::uword k,
                                               int df_0,
                                               arma::uword num_cols,
                                               arma::mat scale_0,
                                               double lambda_0,
                                               arma::vec mu_0){
  arma::mat cluster_data;
  arma::mat variance(num_cols, num_cols);
  arma::vec mu(num_cols);
  
  arma::field<arma::cube> mean_variance_field(2);
  
  arma::cube var_entry = arma::zeros<arma::cube>(num_cols, num_cols, k);
  arma::cube mu_entry = arma::zeros<arma::cube>(num_cols, 1, k);
  
  mean_variance_field(0) = var_entry;
  mean_variance_field(1) = mu_entry;
  
  for (arma::uword j = 1; j < k + 1; j++) {
    // std::cout << "\nj for loop";
    cluster_data = data.rows(find(cluster_labels == j ));
    
    // std::cout << mean_variance_field(1).slice(j - 1) << "\n";
    
    mean_variance_field(0).slice(j - 1) = variance_posterior(
      df_0,
      scale_0,
      lambda_0,
      mu_0,
      cluster_data
    );
    
    // std::cout << "\nVariance sampled\n";
    
    // std::cout << mean_variance_field(1).slice(j - 1) << "\n";
    
    mean_variance_field(1).slice(j - 1) = mean_posterior(mu_0, 
                                                         mean_variance_field(0).slice(j - 1), 
                                                         lambda_0,
                                                         cluster_data);
    
    // std::cout << "\nAccessed cubes";
  }
  return mean_variance_field;
}

// Unused function - considered replacing a for loop in point_comparison, but 
// need multiple non-similar outputs
// Could use an RcppList or something, but think it's best to leave this out
// int update_clusters(arma::uword N,
//                             arma::uword i,
//                             arma::mat data,
//                             arma::uword k,
//                             arma::vec class_weights,
//                             arma::uvec cluster_labels,
//                             arma::cube mu,
//                             arma::cube variance,
//                             bool outlier,
//                             arma::vec global_mean,
//                             arma::mat global_variance,
//                             arma::uword t_df,
//                             arma::uvec fix_vec,
//                             arma::uword num_cols,
//                             arma::uword burn,
//                             arma::uword thinning,
//                             arma::uword eff_count){
//   
//   arma::vec point(num_cols);
//   arma::vec curr_cluster_probs(k);
//   
//   arma::cube class_probs(eff_count, k, N);
//   
//   for (arma::uword jj = 0; jj < N; jj++){
//     // if the current point is not fixed, sample its class
//     point = arma::trans(data.row(jj));
//     
//     curr_cluster_probs = sample_class(point, 
//                                     data,
//                                     k, 
//                                     class_weights, 
//                                     cluster_labels,
//                                     mu,
//                                     variance,
//                                     outlier,
//                                     global_mean,
//                                     global_variance,
//                                     t_df
//     );
//     
//     // std::cout << curr_class_probs << "\n\n";
//     
//     if (i >= burn && (i - burn) % thinning == 0) {
//       // std::cout << "record accessed" << "\n";
//       class_probs.slice(jj).row((i - burn) / thinning) = arma::trans(curr_cluster_probs);
//       
//     }
//     if(! fix_vec[jj]){
//       cluster_labels(jj) = select_class(curr_cluster_probs);
//       // std::cout << "New label\n" << class_labels(jj) << "\n";
//     }
//   }
//   return 0;
// }


// The actual clustering/sampling
// [[Rcpp::export]]
Rcpp::List point_comparison(arma::uword num_iter,
                            arma::vec concentration_0,
                            arma::mat scale_0,
                            arma::uvec class_labels,
                            std::vector<bool> fix_vec,
                            arma::vec mu_0,
                            double lambda_0,
                            arma::mat data,
                            int df_0,
                            arma::uword k,
                            arma::uword burn,
                            arma::uword thinning,
                            bool outlier = false,
                            double t_df = 4.0,
                            bool record_posteriors = false){
  
  // std::cout << "In function";
  arma::uword N;
  arma::uword num_cols;
  N = data.n_rows;
  num_cols = data.n_cols;

  // for use in the outlier distribution
  arma::mat global_variance(num_cols, num_cols);
  
  // std::cout << arma::cov(data) << "\n";
  
  global_variance = 0.5 * arma::cov(data); // Olly's rec
  
  arma::vec global_mean(num_cols);

  global_mean = arma::trans(arma::mean(data, 0));

  arma::vec entropy_cw(num_iter);

  arma::uword eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::uword record_ind;
  
  arma::umat record(N, eff_count);
  record.zeros();
  
  // std::cout << "Record out \n";
  
  arma::mat sim(N, N);
  arma::mat cluster_data;
  
  // Add +1 to k to allow outlier class
  if(outlier){
    k++;
  }
  
  arma::vec class_weights(k);
  
  // These are the lists recording the posterior mean and 
  // variance for each class for each recorded iteration
  ListMatrix variance(eff_count, k);
  ListMatrix mu(eff_count, k);

  arma::vec point;
  
  // std::cout << "Faux output sentence\n";
  
  arma::cube class_probs(eff_count, k, N);
  arma::vec curr_class_probs(k);
  
  // std::cout << "Output sentence\n";
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::field<arma::cube> loc_mu_variance(2);
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    class_weights = class_weight_posterior(concentration_0, class_labels, k);
    
    // std::cout << class_weights << "\n\n";
    // std::cout << "\nENTROPY";
    
    entropy_cw(i) = entropy(class_weights);
   
    // std::cout << "\nBegin sampling parameters\n";
    
    loc_mu_variance = mean_variance_sampling(data,
                           class_labels,
                           k,
                           df_0,
                           num_cols,
                           scale_0,
                           lambda_0,
                           mu_0);
    
    // std::cout << "\nAccessed cubes\n";
    
    for (arma::uword jj = 0; jj < N; jj++){
       // if the current point is not fixed, sample its class
        point = arma::trans(data.row(jj));
        
        curr_class_probs = sample_class(point, 
                                        data,
                                        k, 
                                        class_weights, 
                                        class_labels,
                                        loc_mu_variance(1),
                                        loc_mu_variance(0),
                                        outlier,
                                        global_mean,
                                        global_variance,
                                        t_df
        );
        
        // std::cout << curr_class_probs << "\n\n";
        
        if (i >= burn && (i - burn) % thinning == 0) {
          // std::cout << "record accessed" << "\n";
          record_ind = (i - burn) / thinning;
          class_probs.slice(jj).row(record_ind) = arma::trans(curr_class_probs);
          
        }
        if(! fix_vec[jj]){
        class_labels(jj) = select_class(curr_class_probs);
        // std::cout << "New label\n" << class_labels(jj) << "\n";
      }
    }
    // std::cout << "Labels\n" << class_labels << "\n";
    // std::cout << "Generic message\n" << "\n";
    
    if (i >= burn && (i - burn) % thinning == 0) {
      // std::cout << i << "\n";

      record_ind = (i - burn) / thinning;
      record.col(record_ind) = class_labels;
      // std::cout << "record accessed" << "\n";
      
      if(record_posteriors){
        for(arma::uword j = 0; j < k; j ++){
          
          // std::cout << "Recording params" << j << "\n";
            
          mu(record_ind, j) = loc_mu_variance(1).slice(j);
          variance(record_ind, j) = loc_mu_variance(0).slice(j);
        }
      }

    }
  }

  // std::cout << "Issue is here";
  sim = similarity_mat(record);
  // std::cout << "Y is here";
  // return sim;
  
  if(record_posteriors){
  
    return List::create(Named("similarity") = sim,
                        Named("class_record") = record,
                        Named("mean_posterior") = mu,
                        Named("variance_posterior") = variance,
                        Named("entropy") = entropy_cw,
                        Named("class_prob") = class_probs);
  }
  
  return List::create(Named("similarity") = sim,
                      Named("class_record") = record,
                      Named("entropy") = entropy_cw,
                      Named("class_prob") = class_probs);
}
