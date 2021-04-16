#include likelihoods.stan
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset1; //infectiousness prior to symptoms
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real <lower = 0> max_valid_si;
  real <lower = -100> min_valid_si;
  real <lower = 0> max_invalid_si;
  real <lower = -100> min_invalid_si;
  real <lower = 0> width;
  int M;
  real si_vec[M];
  int first_valid_nu;
}
parameters{
  real <lower = 0, upper = 50> alpha1; // infectious profile parameter
  real <lower = 0, upper = 50> beta1;  // infectious profile parameter
  real <lower = 0, upper = 1> pinvalid; // proportion of SIs invalid
  real <lower = 0, upper = 5> recall; // strength of recall bias
}
model{
  real valid;
  real invalid;
  real denominator_valid;
  matrix[M, N] pdf_mat;
  pinvalid ~ beta(1, 8); // prior on pinvalid
  
  // Do this once when alpha and beta are sampled.
  // Make sure nus are in increasing order. This will work even if
  // some values of nu are repeated
  pdf_mat = pdf_matrix_recall(nu, si_vec, max_shed, offset1, alpha1, beta1,
                       alpha2, beta2, width, first_valid_nu, recall);
  for (n in 1:N) {
    invalid = invalid_lpdf(si[n] | min_invalid_si, max_invalid_si,
                           alpha_invalid, beta_invalid);

    if ((si[n] > offset1) && (nu[n] > offset1)) {
    valid = full_model_lpdf(si[n]| nu[n], max_shed, offset1,
                              recall, alpha1, beta1, alpha2, beta2,
                              width);    
    denominator_valid = sum(col(pdf_mat, n));
    valid = valid - log(denominator_valid);

      target += log_mix(pinvalid, invalid, valid);
      } else {
      target += log(pinvalid) + invalid;
      }
   }
}
