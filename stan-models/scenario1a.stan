#include likelihoods.stan
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> width;
}
parameters{
  real <lower = 0, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter

}
model{
  //alpha1 ~ uniform(1, 10);
  //beta1 ~ uniform(1, 10);
  for (n in 1:N) {
    if (si[n] > 0) {
      target += scenario1a_lpdf(si[n]| max_shed, alpha1, beta1, alpha2, beta2, width);
    } else {
      target += 0;
    }  
  }
}
