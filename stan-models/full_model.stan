#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset1;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> width;
  int M;
  real y_vec[M];
  //real recall; 
  //real alpha1;
  //real beta1;
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;  
  real <lower = 0> max_invalid_si;
  real <lower = -50> min_invalid_si;
  real <lower = 0> max_valid_si;
  
}
parameters{
  real <lower = alpha_invalid, upper = 50> alpha1; // infectious profile parameter
  real <lower = 0, upper = 50> beta1;  // infectious profile parameter
  real <lower = 0, upper = 1> pinvalid;
  real <lower = 0, upper = 5> recall;
}
model{
  real valid;
  real invalid;
  real denominator_valid;
  real denominator_invalid;
  real denominator;
  for (outer in 1:N) {
    invalid = invalid_lpdf(si[outer] | max_invalid_si, min_invalid_si,
                         alpha_invalid, beta_invalid);
    denominator_invalid = (max_valid_si - offset1) / (max_invalid_si - min_invalid_si);

    if ((si[outer] > offset1) && (nu[outer] > offset1)) {
      valid = full_model_lpdf(si[outer]| nu[outer], max_shed, offset1,
                              recall, alpha1, beta1, alpha2, beta2,
                              width);
      denominator_valid = normalising_constant(y_vec, nu[outer], max_shed, 
                                       offset1, recall, alpha1, beta1, 
                                       alpha2, beta2, width);
      denominator = log(pinvalid * denominator_invalid +
                      (1 - pinvalid) * denominator_valid);
      target += log_mix(pinvalid, invalid, valid) - denominator;
      
      } else {
       target += invalid - log(denominator_invalid);
      }
   }
}
