#include likelihoods_gamma.stan
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real offset;
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
  real <lower = 0, upper = 100> alpha1;
  real <lower = 0, upper = 100> beta1;
  real <lower = 0, upper = 5> recall;
  
}
model{
  real valid;
  real denominator_valid;
  real denominator;
  matrix[M, 1] pdf_mat;
  real dummy[1];
  int first_valid_nu = 1;
  // Since this model doesn't need nu, we set nu to be a value larger
  // than max_shed so that the division by F(nu) never takes place.
  dummy[1] = max_shed + 10;
  pdf_mat = pdf_matrix(dummy, si_vec, max_shed, alpha1, beta1, offset,
                       recall, alpha2, beta2, width, first_valid_nu);
  denominator_valid = sum(col(pdf_mat, 1));
  for (n in 1:N) {
    valid = validgamma_lpdf(si[n] |dummy[1], max_shed, alpha1, beta1, offset, 
                         recall, alpha2, beta2, width);
    target += valid - log(denominator_valid);      
  }
}
