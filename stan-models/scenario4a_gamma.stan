#include likelihoods_gamma.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> width;
  int M;
  // Vector of SI from min_invalid_si to max_invalid_si, offset by
  // a small amount to avoid boundary issues.
  real si_vec[M];
}
parameters{
  real <lower = 0> alpha1;
  real <lower = 0> beta1;
}
model{
  real recall = 0;  
  real valid;
  real invalid;
  real denominator_valid;
  matrix[M, N] pdf_mat;
  int first_valid_nu = 1;
  pdf_mat = pdf_matrix(nu, si_vec, max_shed, alpha1, beta1, offset, 
                       recall, alpha2, beta2, width, first_valid_nu);
  for (n in 1:N) {
    valid = validgamma_lpdf(si[n] | nu[n], max_shed, alpha1, beta1, offset, 
                            recall, alpha2, beta2, width);
    denominator_valid = sum(col(pdf_mat, n));
    target += valid - log(denominator_valid);
  }
}
