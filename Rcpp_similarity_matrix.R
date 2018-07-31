#!/usr/bin/env Rscript

library(Rcpp)

# check Rcpp set up correctly
evalCpp("2 + 2")

cppFunction("
bool isOddCpp(int num = 10) {
            bool result = (num % 2 == 1);
            return result;
            }")
isOddCpp(42)

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

# cppFunction('float point_similarity(int point, 
#             int comparison_point,
#             NumericMatrix cluster_record,
#             int num_iter) {
#   float out = 0;
#   // int ncol = cluster_record.ncol();
#   for (int i = 0; i < num_iter; i++){
#     if(cluster_record(point, i) == cluster_record(comparison_point, i)){
#       out++;
#     }
#     
#   }
#   out = out / num_iter;
#   return out;
# }')
# 
# cppFunction('NumericMatrix similarity_mat(NumericMatrix cluster_record){
#   int sample_size = cluster_record.nrow(), num_iter = cluster_record.ncol();
#   NumericMatrix out(sample_size, sample_size);
#   
#   for (int point = 0; point < sample_size - 1; point++){
#     for (int comparison_point = point + 1; comparison_point < sample_size; comparison_point++){
#       out(point, comparison_point) = point_similarity(point, comparison_point, cluster_record, num_iter);
#       out(comparison_point, point) = out(point, comparison_point);
#     }
#   }
#   return out;
# }')

sourceCpp("sampleRcpp.cpp")

cluster_record = matrix(c(1, 1, 0, 1), 2, 2)

similarity_mat(cluster_record)
