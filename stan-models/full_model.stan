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
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;  
  
}
parameters{
  real <lower = alpha_invalid, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
  real <lower = 0, upper = 1> pinvalid;
  real <lower = 0, upper = 5> recall;
}
model{
  real si_inner;
  real denominator;
  real valid;
  real invalid;
  for (outer in 1:N) {
    invalid = invalid_lpdf(si[outer] | max_si, min_si,
                           alpha_invalid, beta_invalid);

    if (si[outer] && nu[outer] > offset1) {
      // For each nu and for sampled recall, we have to normalise over
      // all possible SIs. Hence this loop.
      denominator = 0;
      si_inner = 0;
      while (si_inner <= max_si) {
        denominator = denominator +
          exp(-recall * fabs(si_inner - nu[outer]));
        // No real reason to increment this by width as well,
        // can increase it for speed
        si_inner = si_inner + width; 
      }
      denominator = log(denominator);
      valid = full_model_lpdf(si[outer]| nu[outer], max_shed, offset1,
                              recall, alpha1, beta1, alpha2, beta2,
                              width) - denominator;
      target += log_mix(pinvalid, invalid, valid);
     } else {
       target += invalid;
     }
   }
}
