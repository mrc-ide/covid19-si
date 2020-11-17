#include likelihoods.stan
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real offset1;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter  
}
parameters{
  real <lower = 1, upper = 30> alpha1; // infectious profile parameter
  real <lower = 1, upper = 30> beta1;  // infectious profile parameter
}
model{
  //alpha1 ~ uniform(1, 10);
  //beta1 ~ uniform(1, 10);
  for (n in 1:N) {
    si[n] ~ scenario3a(max_shed, offset1, alpha1, beta1, alpha2, beta2);
  }
}
