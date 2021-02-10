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
  real <lower = 0> max_valid_si;
  real <lower = -100> min_valid_si;
  // min_invalid_si can be a lot smaller than the min_valid_si.
  real <lower = 0> max_invalid_si;
  real <lower = -100> min_invalid_si;
  real <lower = 0> width;
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
  real denominator_invalid;
  real denominator;
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n] | min_invalid_si, max_invalid_si,
                           alpha_invalid, beta_invalid);
    denominator_invalid = (max_valid_si - min_valid_si) /
      (max_invalid_si - min_invalid_si);
    if ((si[n] > offset1) && (nu[n] > offset1)) {
      valid = scenario4a_lpdf(si[n] | nu[n], max_shed, offset1, alpha1,
                              beta1, alpha2, beta2, width);
      denominator_valid = s4_normalising_constant(nu[n], max_shed,
                                                  offset1, alpha1,
                                                  beta1, alpha2, beta2,
                                                  max_valid_si, width);
      denominator = log(pinvalid * denominator_invalid +
                        (1 - pinvalid) * denominator_valid);
      target += log_mix(pinvalid, invalid, valid) - denominator;
    } else {
      target += invalid - log(denominator_invalid);
    }
  }
}
