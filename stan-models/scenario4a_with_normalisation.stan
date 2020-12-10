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
}
parameters{
  real <lower = 0, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
}
model{
  //print("alpha1 ", alpha1);
  //print("beta1 ", beta1);
  real denominator;
  for (n in 1:N) {
    denominator = s4_normalising_constant(y_vec, nu[n], max_shed, 
                                          offset1, alpha1, beta1, 
                                          alpha2, beta2, width);
    target += scenario4a_lpdf(si[n]| nu[n], max_shed, offset1, alpha1,
                              beta1, alpha2, beta2, width) - denominator;;
      //print("target = ", target());

  }
}
