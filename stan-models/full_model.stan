functions {
    // x : SI
  // nu: Delay from symptom onset to isolation
  // max_shed: Maximum possible time of infection of secondary case
  // offset1: Minimum possible time of infection of secondary case
  // Should be negative of time is before symptom onset
  // recall: coefficient of recall bias component
  // alpha 1, beta1 : shape1 and shape2 parameters of infectious profile
  // alpha 2, beta2 : shape and rate parameters of incubation profile
  // width: width of rectangle for aproximation of integral. Should be
  // smaller than the samllest SI.
  // min_si cannot be achieved without an incubation period of 0
  // so when x = min_si, this function will return -Inf
  real full_model_lpdf(real x, real nu, real max_shed, real offset1,
                       real recall, real alpha1, real beta1,
                       real alpha2, real beta2, real width,
                       real max_si, real min_si) {
    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real max_shed_shifted;
    real nu_shifted;
    real denominator;
    
    if (x > max_shed) ulim = max_shed;
    else ulim = x;
    if(ulim > nu) ulim = nu;
    out = 0;
    // the whole infectious profile is shifted
    // right if offset < 0  and left if offset > 0
    max_shed_shifted = max_shed - offset1;
    nu_shifted = nu - offset1;
    //print("nu shifted = ", nu_shifted);
    //print("max shed shifted = ", max_shed_shifted);
    // Mapping s which varies from -offset to ulim to
    // interval (0, 1). We want to integrate from offset to ulim
    s = offset1 + width;
    // Now map it into (0, 1)
    //s = (s - offset1) / max_shed_shifted;
    
    while (s < ulim) {
      inf_density = beta_lpdf((s - offset1)/max_shed_shifted |alpha1, beta1);
      inc_density = gamma_lpdf(x - s|alpha2, beta2); 
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    if(nu < max_shed) {
      out = out - beta_lcdf(nu_shifted / max_shed_shifted|alpha1, beta1);
    }
    //out = out - recall * fabs(x - nu);
    //denominator = normalising_constant(nu, max_si, min_si, recall);
    //print("denominator = ", denominator);
    //denominator = approx_normalising_constant(x, nu, max_si, min_si,
    //                                          recall, width);
    //out = out - denominator;
    //print("out = ", out);
    return out;
  }

  real partial_sum(real[] y_slice, int start, int end,
                   real nu, real max_shed, real offset1,
                   real recall, real alpha1, real beta1,
                   real alpha2, real beta2, real width,
                   real max_si, real min_si) {
    real out = 0;
    for (idx in start:end) {
      out += full_model_lpdf(y_slice[idx] | nu, max_shed, offset1, 
                             recall, alpha1, beta1, alpha2, beta2, 
                             width, max_si, min_si);
    }
    return out;
  }
  
}
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
  int M;
  real y_vec[M];
  //real <lower = 0> alpha_invalid;
  //real <lower = 0> beta_invalid;  
  
}
parameters{
  real <lower = 0, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
  //real <lower = 0, upper = 1> pinvalid;
  real <lower = 0, upper = 5> recall;
}
model{
  //real valid;
  //real invalid;
  real denominator;
  int grainsize = 1;
  for (outer in 1:N) {
    //invalid = invalid_lpdf(si[outer] | max_si, min_si,
    //                     alpha_invalid, beta_invalid);

    //if ((si[outer] > offset1) && (nu[outer] > offset1)) {
    denominator = normalising_constant(y_vec, nu[outer], max_shed, 
                                       offset1, recall, alpha1, beta1, 
                                       alpha2, beta2, width, max_si,
                                       min_si);
    target += full_model_lpdf(si[outer]| nu[outer], max_shed, offset1,
                              recall, alpha1, beta1, alpha2, beta2,
                              width, max_si, min_si) - denominator;
      //target += log_mix(pinvalid, invalid, valid);
      //} else {
      // target += invalid;
      //}
   }
}
