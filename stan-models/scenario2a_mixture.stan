#include likelihoods.stan
data{
  int N; // number of data points
  real si[N];
  real nu[N];
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> max_si;  
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real width;
}
parameters{
  //simplex[2] theta;
  real <lower = 0, upper = 1> pinvalid;
  real <lower = alpha_invalid, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
}
model{
  for (n in 1:N) {
    target += log_mix(pinvalid,
                      beta_lpdf(si[n]/max_si | alpha_invalid, beta_invalid),
                      scenario2a_lpdf(si[n] | nu[n], max_shed, alpha1,
                                              beta1, alpha2, beta2,
                                              width)
                      );
  }
}
