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
float point_similarity(int point, 
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
      out(point, comparison_point) = point_similarity(point, 
                                                      comparison_point,
                                                      cluster_record,
                                                      num_iter);
      out(comparison_point, point) = out(point, comparison_point);
    }
  }
  return out;
}
