#include likelihoods_nf.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> max_invalid_si;
  real <lower = -100> min_invalid_si;
  real <lower = 0> width;
  int M;
  // Vector of SI from min_invalid_si to max_invalid_si, offset by
  // a small amount to avoid boundary issues.
  real si_vec[M];
}
parameters{
  real <lower = 0, upper = 1> pinvalid;
  real <lower = 0> a;
  real <lower = 0> b;
  real <lower = 0, upper = 1> c;
  real <lower = -20, upper = 10> tmax;
  real <lower = 0, upper = 1> pleak;    
}
model{
  real recall = 0;  
  real valid;
  real invalid;
  real denominator_valid;
  matrix[M, N] pdf_mat;
  int first_valid_nu = 1;
  a ~ normal(4, 1);
  b ~ normal(1, 0.5);
  pinvalid ~ beta(4, 10);
  // Do this once when alpha and beta are sampled.
  // Make sure nus are in increasing order. This will work even if
  // some values of nu are repeated
  pdf_mat = leaky_pdf_matrix(nu, si_vec, max_shed, a, b, c, tmax, 
                             recall, alpha2, beta2, pleak, width, first_valid_nu);
  for (n in 1:N) {
      invalid = invalid_lpdf(si[n] | min_invalid_si, max_invalid_si);
      valid = validnf_leaky_lpdf(si[n] | nu[n], max_shed, a, b, c, tmax, 
                                 recall, alpha2, beta2, pleak, width);
      denominator_valid = sum(col(pdf_mat, n));
      valid = valid - log(denominator_valid);
      target += log_mix(pinvalid, invalid, valid);
  }
}
