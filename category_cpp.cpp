# include <RcppArmadillo.h>
# include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;


double point_similarity(int point, 
                        int comparison_point,
                        arma::Mat<int> cluster_record,
                        int num_iter) {
  double out = 0.0;
  
  for (int i = 0; i < num_iter; i++){
    if(cluster_record(point, i) == cluster_record(comparison_point, i)){
      out++;
    }
    
  }
  out = out / num_iter;
  return out;
}

// [[Rcpp::export]]
arma::mat similarity_mat(arma::Mat<int> cluster_record){
  int sample_size = cluster_record.n_rows;
  int num_iter = cluster_record.n_cols;
  arma::mat out(sample_size, sample_size);
  
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


// update the concentration parameter in the Dirichlet distribution
arma::vec concentration_n(arma::vec concentration_0,
                          arma::Col<int> cluster_labels,
                          int num_cat){
  
  int n = cluster_labels.n_elem;
  
  int class_count;
  arma::vec concentration(num_cat);
  
  for (int i = 1; i < num_cat + 1; i++) {
    class_count = 0;
    
    for (int j = 0; j < n; j++ ) {
      if (cluster_labels(j) == i) {
        class_count++;
      }
    }
    
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

// sample parameters for a dirichlet distribution
// [[Rcpp::export]]
arma::vec cluster_weight_posterior(arma::vec concentration_0,
                                   arma::Col<int> cluster_labels,
                                   int num_clusters){
  arma::vec cluster_weight = arma::zeros<arma::vec>(num_clusters);
  
  arma::vec concentration(num_clusters);
  concentration = concentration_n(concentration_0,
                                  cluster_labels,
                                  num_clusters);
  
  
  for (int i = 1; i < num_clusters + 1; i++) {
    
    cluster_weight(i - 1) = Rf_rgamma(arma::as_scalar(concentration(i - 1)), 1);
    
  }
  
  double total_cluster_weight = sum(cluster_weight);
  cluster_weight = cluster_weight / total_cluster_weight;
  return cluster_weight;
}

// count unique entries in a vector
arma::uword unique_counter(arma::vec v){
  std::sort(v.begin(), v.end());
  arma::uword unique_count = std::unique(v.begin(), v.end()) - v.begin();
  return unique_count;
}

// returns a vector of the number of unqieu values in each column
// [[Rcpp::export]]
arma::uvec cat_counter(arma::mat data){
  int num_cols = data.n_cols;
  arma::uvec num_categories(num_cols);
  for(int i = 0; i < num_cols; i++){
    num_categories(i) = unique_counter(data.col(i));
  }
  return num_categories;
}

// find the number of categories in each covariate and declare the appropriate
// matrix to record the associated probabiltieis for each cluster

arma::field<arma::mat> declare_class_probs_field(arma::uvec cat_count,
                                                 arma::uword num_cols,
                                                 arma::uword num_clusters){
  
  arma::field<arma::mat> class_probabilities(num_cols);
  for(arma::uword i = 0; i < num_cols; i++){
    arma::mat phi_j(num_clusters, cat_count(i));
    phi_j.zeros();
    class_probabilities(i) = phi_j;
  }
  return class_probabilities;
}

// Sample the probabilities for each category across all clusters for each covariate
arma::field<arma::mat> sample_class_probabilities(arma::field<arma::mat> class_probabilities,
                                                  arma::field<arma::vec> phi_prior,
                                                  arma::Col<int> cluster_labels,
                                                  arma::uvec cat_count,
                                                  arma::uword num_clusters,
                                                  arma::uword num_cols
                                                    
){

  for(arma::uword j = 0; j < num_cols; j++){
    for(arma::uword k = 0; k < num_clusters; k++){
      class_probabilities(j).row(k) = arma::trans(cluster_weight_posterior(phi_prior(j),
                                                                           cluster_labels,
                                                                           cat_count(j)
                                                                          )
      );
    }
  }
  return class_probabilities;
}

// Sample the cluster membership of point
arma::vec cluster_probabilities(arma::rowvec point,
                                arma::mat data,
                                arma::field<arma::mat> class_probabilities,
                                arma::uword num_clusters,
                                arma::uword num_cols){
  
  // std::cout << "In function cluster_probabilities\n";
  arma::vec probabilities = arma::zeros<arma::vec>(num_clusters);

  // std::cout << "\n\n" << class_probabilities << "\n\n";
  
  for(arma::uword i = 0; i < num_clusters; i++){
    for(arma::uword j = 0; j < num_cols; j++){
      // std::cout << "Adding to probability score thing\n";
      // std::cout << probabilities(i) << "\n";
      // std::cout << point(j) << "\n";
      // 
      // std::cout << "Phi field entry " << j << point(j) << i << "\n";
      // std::cout << class_probabilities(j) << "\n";
      // 
      // std::cout << class_probabilities(j)(i, point(j)) << "\n";
      // 
      // std::cout << "iteration " << i << "\n";
      
      probabilities(i) = probabilities(i) + std::log(class_probabilities(j)(i, point(j)));
    }
    
  }
  probabilities = exp(probabilities - max(probabilities));
  probabilities = probabilities / sum(probabilities);
  
  return probabilities;
}

int cluster_predictor(arma::vec probabilities){
  double u;
  arma::uword pred;
  u = arma::randu<double>( );
  
  pred = 1 + sum(u > cumsum(probabilities));
  return pred;
}


// [[Rcpp::export]]
arma::mat sampling(arma::mat data,
                   arma::field<arma::vec> phi_prior,
                   arma::Col<int> cluster_labels,
                   arma::vec fix_vec,
                   arma::vec cluster_weight_priors,
                   arma::uword num_clusters,
                   arma::uword num_iter,
                   arma::uword burn,
                   arma::uword thinning){
  
  arma::uword n = data.n_rows;
  arma::uword num_cols = data.n_cols;
  
  arma::uvec cat_count(num_cols);
  cat_count = cat_counter(data);
  
  arma::field<arma::mat> class_probabilities(num_cols);
  
  // std::cout << "Declaring field of matrices for class probs\n";
  
  class_probabilities = declare_class_probs_field(cat_count,
                            num_cols,
                            num_clusters);
  
  // std::cout << "Class probabilities declared\n";
  
  arma::vec cluster_weights(num_clusters);
  
  arma::vec curr_class_probs(num_clusters);
  
  int eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  arma::Mat<int> record(n, eff_count);
  record.zeros();
  
  // std::cout << "Reached for loop\n";
  
  for(arma::uword i = 0; i < num_iter; i++){
    
    cluster_weights = cluster_weight_posterior(cluster_weight_priors,
                                               cluster_labels,
                                               num_clusters);
    
    // std::cout << "Cluster weights calculated\n";
    
    class_probabilities = sample_class_probabilities(class_probabilities,
                                                     phi_prior,
                                                     cluster_labels,
                                                     cat_count,
                                                     num_clusters,
                                                     num_cols
                                 
    );
    
    // std::cout << "Class probs calculated\n";
    
    for(arma::uword j = 0; j < n; j++){
      // sample cluster for each point here
      curr_class_probs = cluster_probabilities(data.row(j),
                                               data,
                                               class_probabilities,
                                               num_clusters,
                                               num_cols);
      
      // std::cout << "Cluster sampled\n";
      
      if(! fix_vec(j)){
        cluster_labels(j) = cluster_predictor(curr_class_probs);
      }
      if (i >= burn && (i - burn) % thinning == 0) {

        record.col((i - burn) / thinning) = cluster_labels;
      }
    }
  }
    
  arma::mat sim(n, n);
  sim = similarity_mat(record);
  return sim;  
}
              
