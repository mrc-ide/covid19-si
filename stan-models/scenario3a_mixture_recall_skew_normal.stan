#include likelihoods_skew_normal.stan
data{
  int N; // number of data points
  real si[N];
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
  real <lower = -10, upper = 10> a;
  real <lower = 0, upper = 5> b;
  real <lower = -10, upper = 10> c;
  real <lower = 0, upper = 5> recall;
}
model{
  real valid;
  real invalid;
  real denominator_valid;
  real denominator;
  matrix[M, N] pdf_mat;
  int first_valid_nu = 1;
  pinvalid ~ beta(4, 10);
  pdf_mat = pdf_matrix(nu, si_vec, max_shed, a, b, c, 
                       recall, alpha2, beta2, width, first_valid_nu);
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n]|min_invalid_si, max_invalid_si);
    valid = validnf_lpdf(si[n]|nu[n], max_shed, a, b, c, 
                         recall, alpha2, beta2, width);
    denominator_valid = sum(col(pdf_mat, n));
    valid = valid - log(denominator_valid);      
    target += log_mix(pinvalid, invalid, valid);    
  }
}
