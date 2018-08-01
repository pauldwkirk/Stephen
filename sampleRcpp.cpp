#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
//  timesTwo(42)
// */

// [[Rcpp::export]]
float Cpoint_similarity(int point, 
                       int comparison_point,
                       NumericMatrix cluster_record,
                       int num_iter) {
  float out = 0;
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
NumericMatrix similarity_mat(NumericMatrix cluster_record){
  int sample_size = cluster_record.nrow(), num_iter = cluster_record.ncol();
  NumericMatrix out(sample_size, sample_size);
  
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
// 
// // [[Rcpp::export]]
// NumericMatrix point_comparison(int num_iter,
//                                NumericVector concentration_0,
//                                NumericVector class_labels,
//                                NumericMatrix data,
//                                float df_0,
//                                int k,
//                                int burn,
//                                int thinning
// ){
//   int N = data.nrow();
//   NumericMatrix record(N, num_iter - burn);
//   NumericVector entropy_cw(num_iter);
//   NumericMatrix sim;
//     
//   NumericVector class_weights;
// 
//   Rcpp::List variance;
//   Rcpp::List mu;
// 
//   for(int i = 0; i < num_iter; i++){
// 
//     class_weights = class_weight_posterior(concentration_0, class_labels, k);
//     entropy_cw(i) = entropy(class_weights)
//     for (int j = 0; j < k; j++) {
//       cluster_data <- as.data.frame(data[class_labels == j, ])
//       
//       
//       variance[[j]] <- variance_posterior(
//           df_0,
//           scale_0,
//           lambda_0,
//           mu_0,
//           cluster_data
//       )
//       
//       mu[[j]] <- mean_posterior(mu_0, variance[[j]], lambda_0, cluster_data)
//       */
//     
//     }
//     
//     for (int jj = 0; jj < N; jj++){
//       /***R
//       point <- data[jj, ]
//       
//       class_labels[jj] <- sample_class(point, data, k, class_weights, class_labels,
//                                           mu = mu,
//                                           variance = variance
//       )
// 
//     
//       */
//     }
//     if (i > burn && (i - burn) % thinning == 0) {
//       /***R
//       record[, i - burn] <- t(class_labels)
//       */
//       
//     }
//   }
//   
//   sim = similarity_mat(record);
//   return sim;
// }
// 
// 
// 
// 
// 
//                                
// NumericVector class_weight_posterior(float concentration_0,
//                                      NumericVector class_labels,
//                                      int k){
//   NumericVector class_weight(k);
//   int n = class_labels.size();
//   int class_count = 0;
//   float concentration = 0;
//   float total_class_weight = 0;
//   for (int i = 0; i < k; i++) {
//     for (int j = 0; j < n; j++ )
//     if (class_labels(j) == i) {
//       class_count++;
//     }
//     concentration = concentration_0 + class_count;
//     class_weight(i) = rgamma(1, concentration, 1)(1);
//     total_class_weight += class_weight(i);
//   }
//   for (int i = 0; i < k; i++) {
//     class_weight(i) = class_weight(i) / total_class_weight;
//   }
//   return class_weight;
// }
// 
// float entropy(NumericVector class_weights){
//   
// }
// 
// 
// entropy <- function(class_weights) {
// # Measure of convergence, entropy from Information theory
// # Class weights: vector of class_weights (numbers from unit interval)
//   entropy_components <- class_weights * log(class_weights)
//   entropy_components[is.nan(entropy_components)] <- 0
//   entropy <- sum(entropy_components)
// }