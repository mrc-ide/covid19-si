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
}
parameters{
  real <lower = 0, upper = 30> alpha1; // infectious profile parameter
  real <lower = 0, upper = 30> beta1;  // infectious profile parameter
}
model{
  //print("alpha1 ", alpha1);
  //print("beta1 ", beta1);
  for (n in 1:N) {
    if (si[n] > 0) {
      target += scenario4a_lpdf(si[n]| nu[n], max_shed, offset1, alpha1,
                                beta1, alpha2, beta2, width);
      //print("target = ", target());

    } else {
      target += 0;
    }  
  }
}
