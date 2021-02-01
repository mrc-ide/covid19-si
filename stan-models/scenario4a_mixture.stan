#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset1;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real <lower = 0> max_si;
  real <lower = -100> min_si;
  real <lower = 0> width;
  int M;
  real y_vec[M];  
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
  real denominator;
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n] | min_si, max_si, alpha_invalid, beta_invalid);
    if ((si[n] > offset1) && (nu[n] > offset1)) {
      denominator = pinvalid + 
        (1 - pinvalid) * s4_normalising_constant(y_vec, nu[n], max_shed,
                                                 offset1, alpha1, beta1,
                                                 alpha2, beta2, min_si,
                                                 max_si, width);
      valid = scenario4a_lpdf(si[n] | nu[n], max_shed, offset1, alpha1,
                              beta1, alpha2, beta2, width);
      target += log_mix(pinvalid, invalid, valid) - log(denominator);
    } else {
      target += invalid; //+ log(1 - pinvalid);
    }
  }
}
