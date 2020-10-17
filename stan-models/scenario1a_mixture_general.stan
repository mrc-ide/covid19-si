functions{
  real invalid_lpdf(real x, real max_si, real min_si, real alpha_invalid, real beta_invalid) {
    real out;
    real y;
    y = (x + fabs(min_si))/ (max_si - min_si);
    out = beta_lpdf(y| alpha_invalid, beta_invalid);
    return(out);
  }
  
  real scenario1a_lpdf(real x,
                       real max_shed,
                       real alpha1,      
                       real beta1,   
                       real alpha2,    
                       real beta2,
                       real width) {   

    // out = 0;
    // for i in 0:x
    // out = out + exp(log(f(i)) + log(g(x - i)))
    // return log(out)
    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    // width should be smaller than the smallest SI otherwise
    // x - width will be -ve and gamma_lpdf will be 0, messing up
    // everything.
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // s is the time of infection
    inf_density = beta_lpdf(width/max_shed|alpha1, beta1);
    inc_density = gamma_lpdf(x - width|alpha2, beta2);
    out =  exp(inf_density + inc_density);
    s = 2 * width;
    while(s < ulim) {
      inf_density = beta_lpdf(s/max_shed|alpha1, beta1);
      inc_density = gamma_lpdf(x - s|alpha2, beta2);
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    //print("si = ", x);
    //print("out = ", out);
    return out;
}
}
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real <lower = 0> max_si;
  real <lower = -100> min_si;
  real <lower = 0> width;

}
parameters{
  // simplex[2] theta;
  real <lower = 0, upper = 1> pinvalid;
  real <lower = alpha_invalid, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter
}
model{
  real valid;
  real invalid;
  pinvalid ~ beta(1.5, 5);
  for (n in 1:N) {
    //print("alpha1 = ", alpha1);
    //print("beta1 = ", beta1);    
    //print("valid pdf = ", valid);
    invalid = invalid_lpdf(si[n] | max_si, min_si, alpha_invalid, beta_invalid);
    if (si[n] > 0) {
      valid = scenario1a_lpdf(si[n] | max_shed, alpha1, beta1, alpha2, beta2, width);
      target += log_mix(pinvalid, invalid, valid);    
    } else {
      target += invalid;
    }
  }
}
