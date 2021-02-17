#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset1;
  real max_si;
  real min_si;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> width;
  //int M;
  //real y_vec[M];
  //real recall; 
  //real alpha1;
  //real beta1;
  //real <lower = 0> alpha_invalid;
  //real <lower = 0> beta_invalid;  
  
}
parameters{
  real <lower = 0, upper = 50> alpha1; // infectious profile parameter
  real <lower = 0, upper = 50> beta1;  // infectious profile parameter
  //real <lower = 0, upper = 1> pinvalid;
  real <lower = 0, upper = 5> recall;
}
model{
  //real valid;
  //real invalid;
  //real denominator;
  for (outer in 1:N) {
    //invalid = invalid_lpdf(si[outer] | max_si, min_si,
    //                     alpha_invalid, beta_invalid);

    //if ((si[outer] > offset1) && (nu[outer] > offset1)) {
    //denominator = normalising_constant(y_vec, nu[outer], max_shed, 
                                       //offset1, recall, alpha1, beta1, 
                                       //alpha2, beta2, width, max_si,
                                       //min_si);
    target += full_model_lpdf(si[outer]| nu[outer], max_shed, offset1,
                              recall, alpha1, beta1, alpha2, beta2,
                              width, max_si, min_si); //- denominator;
      //target += log_mix(pinvalid, invalid, valid);
      //} else {
      // target += invalid;
      //}
   }
}
