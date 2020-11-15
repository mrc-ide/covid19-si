#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval
  real nu[N]; // Time of isolation
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> width;
}
parameters{
  real <lower = 0, upper = 10> alpha1; // infectious profile parameter
  real <lower = 0, upper = 10> beta1; // // infectious profile parameter
} 
model{
  //alpha1 ~ uniform(1, 2000);
  //beta1 ~ uniform(1, 2000);
  for (n in 1:N) {
    si[n] ~ scenario2a(nu[n], max_shed, alpha1, beta1, alpha2, beta2, width);
  }
}
