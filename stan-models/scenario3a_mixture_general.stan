#include likelihoods.stan
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real offset1;
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
}
parameters{
  // simplex[2] theta;
  real <lower = 0, upper = 1> pinvalid;
  real <lower = alpha_invalid, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
}
model{
  real valid;
  real invalid;
  real denominator_valid;
  real denominator;
  matrix[M, 1] pdf_mat;
  real dummy[1];
  pinvalid ~ beta(4, 10);
  // Since this model doesn't need nu, we set nu to be a value larger
  // than max_shed so that the division by F(nu) never takes place.
  dummy[1] = max_shed + 10;
  pdf_mat = pdf_matrix(dummy, si_vec, max_shed, offset1, 
                       alpha1, beta1, alpha2, beta2, width, first_valid_nu);
  denominator_valid = sum(col(pdf_mat, 1));
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n] | min_invalid_si, max_invalid_si,
                           alpha_invalid, beta_invalid);
    if (si[n] > offset1) {
      valid = scenario3a_lpdf(si[n] | max_shed, offset1, alpha1, beta1,
                              alpha2, beta2, width);
      valid = valid - log(denominator_valid);      
      target += log_mix(pinvalid, invalid, valid);    
    } else {
      target += invalid + log(pinvalid);
    }
  }
}
