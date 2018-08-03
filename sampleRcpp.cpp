# include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//#include <Rcpp.h>
// using namespace Rcpp;
// using namespace arma;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
//  timesTwo(42)
// */

// [[Rcpp::export]]
double Cpoint_similarity(int point, 
                       int comparison_point,
                       arma::mat cluster_record,
                       int num_iter) {
  double out = 0.0;
  // int ncol = cluster_record.ncol();
  for (int i = 0; i < num_iter; i++){
    if(cluster_record(point, i) == cluster_record(comparison_point, i)){
      out++;
    }
    
  }
  out = out / num_iter;
  return out;
}

// [[Rcpp::export]]
arma::mat similarity_mat(arma::mat cluster_record){
  int sample_size = cluster_record.n_rows;
  int num_iter = cluster_record.n_cols;
  arma::mat out(sample_size, sample_size);
  
  for (int point = 0; point < sample_size; point++){ // if not doing diagonal, restrict to sample size - 1
    for (int comparison_point = point; // + 1; 
         comparison_point < sample_size;
         comparison_point++){
      out(point, comparison_point) = Cpoint_similarity(point, 
                                                      comparison_point,
                                                      cluster_record,
                                                      num_iter);
      out(comparison_point, point) = out(point, comparison_point);
    }
  }
  return out;
}


// [[Rcpp::export]]
arma::mat S_n(arma::mat data,
              arma::vec sample_mean,
              int sample_size,
              int num_cols){
  arma::mat sample_covariance(num_cols, num_cols);
  if(sample_size > 0){
    for(int i = 0; i < sample_size; i++){
      arma::vec col_i = data.col(i);
      sample_covariance = (sample_covariance 
                             + ((col_i - sample_mean) 
                                  * arma::trans(col_i - sample_mean)
                                  )
                             );
      
    }
  }
  return sample_covariance;
}

// [[Rcpp::export]]
arma::vec mean_n(double lambda_0,
                     arma::vec mu_0,
                     int sample_size,
                     arma::vec sample_mean){
  arma::vec mu_n;
  mu_n = ((lambda_0 * mu_0 + sample_size * sample_mean)
            / (lambda_0 + sample_size));
  return mu_n;
}

// [[Rcpp::export]]
arma::mat scale_n(arma::mat scale_0,
                      arma::vec mu_0,
                      double lambda_0,
                      arma::mat sample_covariance,
                      int sample_size,
                      arma::vec sample_mean){
  arma::mat scale_out;
  scale_out = (scale_0
                + sample_covariance
                + ((lambda_0 * sample_size) / (lambda_0 + sample_size))
                * (sample_mean - mu_0) * arma::trans(sample_mean - mu_0)
  ); 
  return scale_out;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n,
                      arma::vec mu,
                      arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() +
    Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::vec mean_posterior(arma::vec mu_0,
                             arma::mat variance,
                             double lambda_0,
                             arma::mat data){
  int ncols = data.n_cols;
  arma::vec sample_mean(ncols);
  
  int sample_size = data.n_rows;
  if (sample_size > 0){
    sample_mean = arma::mean(data, 0);
  }
  
  double lambda_n = lambda_0 + sample_size;
  arma::vec mu_n;
  mu_n = mean_n(lambda_0, mu_0, sample_size, sample_mean);
  arma::mat variance_n = variance / lambda_n;
  
  arma::vec mu;
  mu = mvrnormArma(1, mu_n, variance_n);
  return mu;
  
}

// [[Rcpp::export]]
arma::vec class_weight_posterior(double concentration_0,
                                 arma::vec class_labels,
                                 int k){
  arma::vec class_weight(k);
  int n = class_labels.n_elem;
  int class_count = 0;
  double concentration = 0.0;
  double total_class_weight = 0.0;
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < n; j++ )
      if (class_labels(j) == i) {
        class_count++;
      }
      concentration = concentration_0 + class_count;
      class_weight(i) = rgamma(1, concentration, 1)(1);
      total_class_weight += class_weight(i);
  }
  for (int i = 0; i < k; i++) {
    class_weight(i) = class_weight(i) / total_class_weight;
  }
  return class_weight;
}

// [[Rcpp::export]]
double entropy(arma::vec class_weights){
  int n = class_weights.n_elem;
  arma::vec entropy_components(n);
  entropy_components = class_weights * log(class_weights);
  if (entropy_components.has_nan()){
    for(int i = 0; i < n; i++){
      entropy_components(i) = 0;
    }
  }
  double entropy_out = sum(entropy_components);
  return entropy_out;
}

// [[Rcpp::export]]
arma::mat variance_posterior(int df_0,
                             arma::mat scale_0,
                             double lambda_0,
                             arma::vec mu_0,
                             arma::mat data){
  int sample_size = data.n_rows, num_cols = data.n_cols;
  arma::vec sample_mean(num_cols);
  
  int df_n = df_0 + sample_size;
  
  arma::mat scale_n_value(num_cols, num_cols);
  arma::mat sample_covariance(num_cols, num_cols);
  arma::mat variance(num_cols, num_cols);
  
  if (sample_size > 0){
    sample_mean = mean(data, 0);
  } else{
    sample_mean.fill(0.0);
  }
  
  sample_covariance = S_n(data, sample_mean, sample_size, num_cols);
  
  scale_n_value = scale_n(
              scale_0,
              mu_0,
              lambda_0,
              sample_covariance,
              sample_size,
              sample_mean
          );
  
  variance = arma::iwishrnd(scale_n_value, df_n);
  
  return variance;
   
} 

// [[Rcpp::export]]
int sample_class(arma::vec point,
                   arma::mat data,
                   int k,
                   arma::vec class_weights,
                   arma::vec class_labels,
                   arma::cube mu,
                   arma::cube variance){
  arma::vec prob(k);
  double curr_weight;
  double exponent;
  double log_likelihood;
  arma::uvec curr_sample;
  double log_determinant;
  arma::vec prob_vec(k);
  double u;
  int pred;
  arma::uvec count_probs;
  // arma::mat variance_loc;
  // arma::vec mu_loc;
  for(int i = 0; i < k; i++){
    curr_weight = log(class_weights(i));
    curr_sample = find(class_weights == i);
    
    // mu_loc = mu(i);
    // variance_loc = variance(i);
    
    exponent = -0.5 * arma::as_scalar((point - mu.slice(i)) 
                                        * arma::inv(variance.slice(i))
                                        * arma::trans(point - mu.slice(i)));
                                        
    log_determinant = arma::log_det(variance.slice(i)).real();
    log_likelihood = -0.5 * log_determinant  + exponent;
    prob_vec(i) = curr_weight + log_likelihood;
  } 
  prob_vec = exp(prob_vec - max(prob_vec));
  prob_vec = prob_vec / sum(prob_vec);
  u = arma::randu<double>( );
  
  count_probs = arma::find(prob_vec > u);
  pred = count_probs.n_elem;
  return pred;
}

 
// Rcpp::List

// [[Rcpp::export]
arma::mat point_comparison(int num_iter,
                            double concentration_0,
                            arma::mat scale_0,
                            arma::vec class_labels,
                            arma::vec mu_0,
                            double lambda_0,
                            arma::mat data,
                            int df_0,
                            int k,
                            int burn,
                            int thinning){
    int N, num_cols;
    N = data.n_rows;
    num_cols = data.n_cols;
    arma::mat record(N, num_iter - burn);
    arma::vec entropy_cw(num_iter);
    arma::mat sim(N, N);
    arma::mat cluster_data;
    arma::vec class_weights(k);
    
    arma::cube variance(num_cols, num_cols, k);
    arma::cube mu(1, num_cols, k);
    
    arma::rowvec point;
    
    for(int i = 0; i < num_iter; i++){
      
      class_weights = class_weight_posterior(concentration_0, class_labels, k);
      
      entropy_cw(i) = entropy(class_weights);
        for (int j = 0; j < k; j++) {
          cluster_data = data.elem( find(class_labels == j ));
          
          variance.slice(j) = variance_posterior(
              df_0,
              scale_0,
              lambda_0,
              mu_0,
              cluster_data
          );
          
          mu.slice(j) = mean_posterior(mu_0, variance.slice(j), lambda_0, cluster_data);
          
            
        }
        
        for (int jj = 0; jj < N; jj++){
          point = data.row(jj);
          
          class_labels(jj) = sample_class(point, data, k, class_weights, class_labels,
                                          mu,
                                          variance
          );
  
        }
        if (i > burn && (i - burn) % thinning == 0) {
          
          record.col(i - burn) = trans(class_labels);
        }
    }
    
    sim = similarity_mat(record);
    return sim;
    // return Rcpp::List::create(Rcpp::named("Similarity") = sim);
  }