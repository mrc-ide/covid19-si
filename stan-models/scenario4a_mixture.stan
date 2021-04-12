#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real <lower = 0> max_valid_si;
  real <lower = -100> min_valid_si;
  // min_invalid_si can be a lot smaller than the min_valid_si.
  real <lower = 0> max_invalid_si;
  real <lower = -100> min_invalid_si;
  real <lower = 0> width;
  int M;
  // Vector of SI from min_invalid_si to max_invalid_si, offset by
  // a small amount to avoid boundary issues.
  real si_vec[M];
  int first_valid_nu;
  //real t_max;
}
parameters{
  real <lower = 0, upper = 1> pinvalid;
  real <lower = 0> a;
  real <lower = 0> b;
  real <lower = 0, upper = 1> c;
  real <lower = -20, upper = 10> tmax;
}
model{
  real valid;
  real invalid;
  real denominator_valid;
  matrix[M, N] pdf_mat;
  a ~ normal(4.28, 0.74);
  b ~ normal(1.44, 0.12);
  pinvalid ~ beta(1, 8);
  // Do this once when alpha and beta are sampled.
  // Make sure nus are in increasing order. This will work even if
  // some values of nu are repeated
  pdf_mat = pdf_matrix(nu, si_vec, max_shed, a, b, c, tmax, 
                       alpha2, beta2, width, first_valid_nu);
  for (n in 1:N) {
      invalid = invalid_lpdf(si[n] | min_invalid_si, max_invalid_si,
                           alpha_invalid, beta_invalid);
      valid = scenario4a_lpdf(si[n] | nu[n], max_shed, a, b, c, tmax, 
                              alpha2, beta2, width);
      denominator_valid = sum(col(pdf_mat, n));
      valid = valid - log(denominator_valid);
      target += log_mix(pinvalid, invalid, valid);
  }
}
