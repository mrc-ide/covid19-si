#include likelihoods.stan
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real offset;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real <lower = 0> max_si;
  real <lower = -100> min_si;
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
  for (n in 1:N) {
    //print("alpha1 = ", alpha1);
    //print("beta1 = ", beta1);    
    //print("valid pdf = ", valid);
    invalid = invalid_lpdf(si[n] | max_si, min_si, alpha_invalid, beta_invalid);
    if (si[n] > -offset) {
      valid = scenario3a_v2_lpdf(si[n] | max_shed, offset, alpha1, beta1, alpha2, beta2, width);
      target += log_mix(pinvalid, invalid, valid);    
    } else {
      target += invalid;
    }
  }
}
