#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset1;  
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> max_invalid_si;
  real <lower = -100> min_invalid_si;
  real <lower = 0> width;
  int M;
  // Vector of SI from min_invalid_si to max_invalid_si, offset by
  // a small amount to avoid boundary issues.
  real si_vec[M];
  int first_valid_nu;
}
parameters{
  real <lower = 0, upper = 1> pinvalid;
  real <lower = 0, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
}
model{
  real recall = 0;  
  real valid;
  real invalid;
  real denominator_valid;
  matrix[M, N] pdf_mat;
  pinvalid ~ beta(4, 10);
  // Do this once when alpha and beta are sampled.
  // Make sure nus are in increasing order. This will work even if
  // some values of nu are repeated
  pdf_mat = pdf_matrix_recall(nu, si_vec, max_shed, offset1, 
                              recall, alpha1, beta1,
                              alpha2, beta2, first_valid_nu, width);
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n] | min_invalid_si, max_invalid_si);
    if ((si[n] > offset1) && (nu[n] > offset1)) {
      valid = valid_beta_lpdf(si[n] | nu[n], max_shed, offset1, 
                              recall, alpha1, beta1,
                              alpha2, beta2, width);
      denominator_valid = sum(col(pdf_mat, n));      
      valid = valid - log(denominator_valid);      
      target += log_mix(pinvalid, invalid, valid);          
    } else {
      target += log(pinvalid) + invalid + log(1 - pinvalid);
    }
  }
}
