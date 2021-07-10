#include likelihoods_skew_normal.stan
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> width;
  int M;
  real <lower = 0> max_invalid_si;
  real <lower = -100> min_invalid_si;
  // Vector of SI from min_invalid_si to max_invalid_si, offset by
  // a small amount to avoid boundary issues.
  real si_vec[M];
}
parameters{
  // a is what stanc calls xi
  // b is what stan calls omega
  // c is what stan calls alpha;
  real <lower = -10, upper = 10> a;
  real <lower = 0, upper = 5> b;
  real <lower = -10, upper = 10> c;
  real <lower = 0, upper = 1> pinvalid;
}
model{
  real recall = 0;
  real invalid;
  real valid;
  real denominator_valid;
  real denominator;
  matrix[M, 1] pdf_mat;
  real dummy[1];
  int first_valid_nu = 1;
  pinvalid ~ beta(4, 10);
  // Since this model doesn't need nu, we set nu to be a value larger
  // than max_shed so that the division by F(nu) never takes place.
  dummy[1] = max_shed + 10;
  pdf_mat = pdf_matrix(dummy, si_vec, max_shed, a, b, c, 
                       recall, alpha2, beta2, width, first_valid_nu);
  denominator_valid = sum(col(pdf_mat, 1));
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n]|min_invalid_si, max_invalid_si);
    valid = validnf_lpdf(si[n] |dummy[1], max_shed, a, b, c, 
                         recall, alpha2, beta2, width);
    valid = valid - log(denominator_valid);      
    target += log_mix(pinvalid, invalid, valid);
  }
}
