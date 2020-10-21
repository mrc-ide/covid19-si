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
  real <lower = 0, upper = 10> beta1; // infectious profile parameter
  real <lower = 0, upper = 2> recall;
} 
model{
  // This is supposed to normalise the probability;
  // For a given SI, we want the probability to be normalised over all
  // possible nu.
  //
  real numerator;
  real denominator;
   for (outer in 1:N) {
     // For each SI and for sampled recall, we have to normalise over
     // all possible nus. Hence this loop.
     denominator = 0;
     for (inner in 1:N) {
       denominator = denominator +
         exp(-recall * fabs(si[outer] - nu[inner]));
     }
     denominator = log(denominator);
     target += scenario2a_with_recall_lpdf(si[outer]| nu[outer], max_shed, alpha1, beta1, recall, alpha2, beta2, width) - denominator;
  }
}

