#include likelihoods.stan
data{
  int N; // number of data points
  real si[N];
  real max_shed;
  real width;
}
parameters{
  real <lower = 0, upper = 30> alpha1; // infectious profile parameter
  real <lower = 0, upper = 30> beta1; // // infectious profile parameter
  real <lower = 0, upper = 30> alpha2; // incubation period parameter
  real <lower = 0, upper = 30> beta2; // incubation period parameter
}
model{
  for (n in 1:N) {
    si[n] ~ scenario1b(max_shed, alpha1, beta1, alpha2, beta2, width);
  }

}
