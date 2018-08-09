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

arma::mat S_n(arma::mat data,
              arma::vec sample_mean,
              int sample_size,
              int num_cols){
  arma::mat sample_covariance(num_cols, num_cols);
  sample_covariance.zeros();
  if(sample_size > 0){
    for(int i = 0; i < sample_size; i++){
      
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

arma::vec mean_n(double lambda_0,
                 arma::vec mu_0,
                 int sample_size,
                 int num_cols,
                 arma::vec sample_mean){
  arma::vec mu_n(num_cols);
  mu_n = ((lambda_0 * mu_0 + sample_size * sample_mean)
            / (lambda_0 + sample_size));
  return mu_n;
}

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

arma::mat mvrnormArma(int n,
                      arma::vec mu,
                      arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() +
    Y * arma::chol(sigma);
}

arma::vec mean_posterior(arma::vec mu_0,
                         arma::mat variance,
                         double lambda_0,
                         arma::mat data){
  int ncols = data.n_cols;
  arma::vec sample_mean(ncols); sample_mean.zeros();
  arma::vec mu_out(ncols);
  arma::vec mu_n(ncols);
  arma::mat variance_n(ncols, ncols);
  
  int sample_size = data.n_rows;
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
  // arma::vec y(ncols);
  // y = arma::conv_to<arma::vec>::from(x.row(0));
  
  // std::cout << x << "\n";
  // std::cout << y;
  
  mu_out = arma::conv_to<arma::vec>::from(x.row(0));
  
  // std::cout << x << "\n";
  // std::cout << x.row(0);
  // 
  // mu_out = trans(mvrnormArma(1, mu_n, variance_n).row(0));
  
  return mu_out;
  
}

arma::vec concentration_n(arma::vec concentration_0,
                          arma::Col<int> class_labels,
                          int k){
  
  int n = class_labels.n_elem;
  
  // std::cout << n ;
  
  int class_count;
  arma::vec concentration(k);
  
  for (int i = 1; i < k + 1; i++) {
    class_count = 0;
    
    for (int j = 0; j < n; j++ ) {
      if (class_labels(j) == i) {
        class_count++;
      }
    }
    
    concentration(i - 1) = arma::as_scalar(concentration_0(i - 1)) + class_count;
  }
  return concentration;
}

arma::vec class_weight_posterior(arma::vec concentration_0,
                                 arma::Col<int> class_labels,
                                 int k){
  arma::vec class_weight = arma::zeros<arma::vec>(k);

  // std::cout << "\n\n" << class_weight << "\n\n";
  // 
  // std::cout << concentration_0;
  
  arma::vec concentration(k);
  concentration = concentration_n(concentration_0,
                                  class_labels,
                                  k);
  
  
  for (int i = 1; i < k + 1; i++) {
    
    class_weight(i - 1) = Rf_rgamma(arma::as_scalar(concentration(i - 1)), 1);

  }
  // for (int i = 0; i < k; i++) {
  //   class_weight(i) = class_weight(i) / total_class_weight;
  // }
  // std::cout << "Finished loop\n";
  
  // std::cout << class_weight;
  
  double total_class_weight = sum(class_weight);
  class_weight = class_weight / total_class_weight;
  return class_weight;
}



// [[Rcpp::export]]
double entropy(arma::vec class_weights){
  int n = class_weights.n_elem;
  arma::vec entropy_components(n);
  // std::cout << "\nDeclared\n";
  
  
  for(int i = 0; i < n; i++){
    entropy_components(i) = - class_weights(i) * log(class_weights(i));
    if (entropy_components.has_nan()){
      entropy_components(i) = 0;
    }
  }
  // std::cout << "Inter";
  double entropy_out = sum(entropy_components);
  return entropy_out;
}

arma::mat variance_posterior(int df_0,
                             arma::mat scale_0,
                             double lambda_0,
                             arma::vec mu_0,
                             arma::mat data){
  
  int sample_size = data.n_rows, num_cols = data.n_cols;
  
  // std::cout <<"\nSample size:\n" << sample_size 
  //   << "\nDimensionality:\n" << num_cols;
  
  // std::cout << "\nCluster data:\n" << data;
  
  arma::vec sample_mean(num_cols); sample_mean.zeros();
  
  int df_n = df_0 + sample_size;
  
  // std::cout << "Reached another dec\n";
  
  arma::mat scale_n_value(num_cols, num_cols);
  arma::mat sample_covariance(num_cols, num_cols);
  arma::mat variance_out(num_cols, num_cols);
  
  // std::cout << "have an issue?\n";
  
  if (sample_size > 0){
    
    sample_mean = arma::trans(arma::mean(data, 0));
    
  } else{
    sample_mean.fill(0.0);
  }
  
  // std::cout << "\nSample covariance reached\n";
  
  sample_covariance = S_n(data, sample_mean, sample_size, num_cols);
  
  // std::cout << sample_covariance << "\n";
  
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
  
  // std::cout << scale_n_value;
  // 
  // std::cout << "\nVariance sampling reached\n";
  // std::cout << arma::iwishrnd(scale_n_value, df_n);
  variance_out = arma::iwishrnd(scale_n_value, df_n);
  
  return variance_out;
  
} 

arma::vec sample_class(arma::vec point,
                       arma::mat data,
                       int k,
                       arma::vec class_weights,
                       arma::Col<int> class_labels,
                       arma::cube mu,
                       arma::cube variance,
                       bool outlier = false,
                       arma::vec global_mean = arma::zeros<arma::vec>(1),
                       arma::mat global_variance = arma::zeros<arma::mat>(1, 1),
                       double t_df = 4.0){
  
  // arma::vec prob(k);
  double curr_weight;
  double exponent;
  double log_likelihood;
  
  // arma::uvec curr_sample;
  
  double log_det;
  arma::vec prob_vec(k);

  arma::uvec count_probs;
  
  int d = data.n_cols;;
  // double df = 0.0;
  // arma::vec concentration(k); concentration.zeros();
  
  // if(outlier){
    // 
    // concentration = concentration_n(concentration_0,
    //                                 class_labels,
    //                                 k);
    // df = concentration(k - 1) - (double)d + 1.0;
    // std::cout << df << "\n";
  // }
  
  
  for(int i = 1; i < k + 1; i++){
    curr_weight = log(class_weights(i - 1));

    if(outlier && i == k){

      exponent = arma::as_scalar(arma::trans(point - global_mean) 
                                   * arma::inv(global_variance)
                                   * (point - global_mean));
                                   
      // std::cout << exponent << "\n";
                                   
      log_det = arma::log_det(global_variance).real();
      
      // std::cout << log_det << "\n";
      
      // log_likelihood = (lgamma((d + t_df)/2.0) - (lgamma(t_df/2.0) +
      //   0.5 * log_det + d/2.0 * std::log(M_PI * t_df)) - (t_df + d/2.0) * 
      //   log(1 + (1/t_df) * exponent));
      
      log_likelihood = (lgamma((t_df + d)/2.0) - lgamma(t_df/2.0) + d/2.0 * log(t_df)
                        + d/2.0 * log(M_PI) + 0.5 * log_det - (t_df + d/2.0) *
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

int select_class(arma::vec prob_vec){
  
  // std::cout << prob_vec << "\n\n";
  double u;
  int pred;
  u = arma::randu<double>( );
  
  // count_probs = arma::find(prob_vec > u);
  // pred = count_probs.n_elem;
  
  pred = 1 + sum(u > cumsum(prob_vec));
  return pred;
}

// [[Rcpp::export]]
Rcpp::List point_comparison(int num_iter,
                            arma::vec concentration_0,
                            arma::mat scale_0,
                            arma::Col<int> class_labels,
                            std::vector<bool> fix_vec,
                            arma::vec mu_0,
                            double lambda_0,
                            arma::mat data,
                            int df_0,
                            int k,
                            int burn,
                            int thinning,
                            bool outlier = false,
                            double t_df = 4.0){
  
  // std::cout << "In function";
  int N;
  int num_cols;
  N = data.n_rows;
  num_cols = data.n_cols;
  
  // std::cout << "Declaration of int";
  
  
  // for use in the outlier distribution
  arma::mat global_variance(num_cols, num_cols);
  
  // std::cout << arma::cov(data) << "\n";
  
  global_variance = 0.5 * arma::cov(data); // Olly's rec
  
  arma::vec global_mean(num_cols);
  
  // std::cout << arma::mean(data, 0) << "\n";
  
  // std::cout << num_cols << "\n";
  
  global_mean = arma::trans(arma::mean(data, 0));
  
  // std::cout << "Declared global var \n";
  
  arma::vec entropy_cw(num_iter);
  
  // arma::vec class_probs(k);

  
  // std::cout << "Declaration of record";
  
  int eff_count = ceil((double)(num_iter - burn) / (double)thinning);
  int record_ind;
  
  // std::cout << eff_count << "\n";
  // std::cout << (double)(num_iter - burn) / (double)thinning << "\n";
  
  arma::Mat<int> record(N, eff_count);
  record.zeros();
  
  // std::cout << "Record out \n";
  
  arma::mat sim(N, N);
  arma::mat cluster_data;
  
  // Add +1 to k to allow outlier class
  if(outlier){
    k++;
  }
  
  arma::vec class_weights(k);
  
  // These are the fields containing cubes recording the posterior mean and 
  // variance for each class for each recorded iteration
  // arma::field<arma::mat> variance(eff_count, k);
  // arma::field<arma::mat> mu(eff_count, k);
  
  ListMatrix variance(eff_count, k);
  ListMatrix mu(eff_count, k);
  
  // These are the local cubes of posterior mean and variance overwritten each
  // iteration
  arma::cube loc_variance(num_cols, num_cols, k);
  arma::cube loc_mu(num_cols, 1, k);
  
  arma::vec point;
  
  // std::cout << "Faux output sentence\n";
  
  arma::field<arma::vec> curr_class_probs(N, 1);
  List class_probs(eff_count);
  
  // std::cout << "Output sentence\n";
  
  for(int i = 0; i < num_iter; i++){
    // std::cout << "Output sentence";
    class_weights = class_weight_posterior(concentration_0, class_labels, k);
    // std::cout << class_weights << "\n\n";
    // std::cout << "\nENTROPY";
    entropy_cw(i) = entropy(class_weights);
    for (int j = 1; j < k + 1; j++) {
      // std::cout << "\nj for loop";
      cluster_data = data.rows(find(class_labels == j ));
      
      // std::cout << scale_0 << "\n";
      
      loc_variance.slice(j - 1) = variance_posterior(
        df_0,
        scale_0,
        lambda_0,
        mu_0,
        cluster_data
      );
      
      // std::cout << "\nVariance sampled\n";
      
      loc_mu.slice(j - 1) = mean_posterior(mu_0, loc_variance.slice(j - 1), lambda_0, cluster_data);
      
      // std::cout << "\nAccessed cubes";
    }
    
    for (int jj = 0; jj < N; jj++){
       // if the current point is not fixed, sample its class
        point = arma::trans(data.row(jj));
        
        // std::cout << "Point initialised\n" << "Iteration:\n" << jj << "\n";
        
        // std::cout << class_labels(jj) << "\nCheck index \n";
        
        curr_class_probs(jj, 0) = sample_class(point, 
                                            data,
                                            k, 
                                            class_weights, 
                                            class_labels,
                                            loc_mu,
                                            loc_variance,
                                            outlier,
                                            global_mean,
                                            global_variance,
                                            t_df
        );
      
      if(! fix_vec[jj]){
        class_labels(jj) = select_class(curr_class_probs(jj, 0));
        // std::cout << "New label\n" << class_labels(jj) << "\n";
      }
    }
    // std::cout << "Labels\n" << class_labels << "\n";
    // std::cout << "Generic message\n" << "\n";
    
    if (i >= burn && (i - burn) % thinning == 0) {
      // std::cout << i << "\n";
      // std::cout << (i - burn) / thinning << "\n";
      record_ind = (i - burn) / thinning;
      record.col(record_ind) = class_labels;
      // std::cout << "record accessed" << "\n";
      class_probs(record_ind) = curr_class_probs;
      
      for(int j = 0; j < k; j ++){
        // variance
        // std::cout << "Recording params" << j << "\n";
        mu(record_ind, j) = loc_mu.slice(j);
        variance(record_ind, j) = loc_variance.slice(j);
      }
      
    }
  }
  // std::cout << "Record\n" << record << "\n";
  // for(int i = 0; i <num_iter - burn; i++){
  //   std::cout << record.col(i) << "\n";
  // }
  
  
  
  // std::cout << "Issue is here";
  sim = similarity_mat(record);
  // std::cout << "Y is here";
  // return sim;
  return List::create(Named("similarity") = sim,
                      Named("mean_posterior") = mu,
                      Named("variance_posterior") = variance,
                      Named("entropy") = entropy_cw,
                      Named("class_prob") = class_probs);
}
