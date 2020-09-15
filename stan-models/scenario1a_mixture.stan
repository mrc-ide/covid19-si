functions{
  real scenario1a_lpdf(real x,
                      real max_shed,
                      real alpha1,      
                      real beta1,   
                      real alpha2,    
                      real beta2) {   

    // out = 0;
    // for i in 0:x
    // out = out + exp(log(f(i)) + log(g(x - i)))
    // return log(out)
    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real width;
    // width should be smaller than the smallest SI otherwise
    // x - width will be -ve and gamma_lpdf will be 0, messing up
    // everything.
    width = 0.01;
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // s is the time of infection
    inf_density = beta_lpdf(width/max_shed|alpha1, beta1) + log(width) - log(max_shed);
    inc_density = gamma_lpdf(x - width|alpha2, beta2) + log(width);
    out =  exp(inf_density + inc_density);
    s = 2 * width;
    while(s < ulim) {
      inf_density = beta_lpdf(s/max_shed|alpha1, beta1) + log(width) - log(max_shed);
      inc_density = gamma_lpdf(x - s|alpha2, beta2)+ log(width);
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
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
  real <lower = 0> min_si;
}
parameters{
  //simplex[2] theta;
  real <lower = 0, upper = 1> pinvalid;
  real <lower = 1, upper = 100> alpha1; // infectious profile parameter
  real <lower = 1, upper = 100> beta1;  // infectious profile parameter
}
model{
  //vector[2] log_theta = log(theta);
  //real max_si = max(si) + 0.001; // so that the max si is not mapped to 1
  
  for (n in 1:N) {
    //vector[2] lps = log_theta;
    target += log(pinvalid) +
      beta_lpdf(si[n]/ max_si| alpha_invalid, beta_invalid);
    //target += log(pinvalid) + uniform_lpdf(si[n]|min_si, max_si);
    target += log(1 - pinvalid) +
      scenario1a_lpdf(si[n] | max_shed, alpha1, beta1, alpha2, beta2);
  }
}

