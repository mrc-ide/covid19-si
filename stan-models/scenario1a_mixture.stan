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
    width = 0.1;
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // s is the time of infection
    inf_density = beta_lpdf(width|alpha1, beta1) + log(width);
    inc_density = gamma_lpdf(x - width|alpha2, beta2) + log(width);
    out =  exp(inf_density + inc_density);
    s = 0.2;
    while(s < ulim) {
      inf_density = beta_lpdf(s/max_shed|alpha1, beta1) + log(width);
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
}
parameters{
  real <lower = 1> alpha1; // infectious profile parameter
  real <lower = 1> beta1;  // infectious profile parameter
  real <lower = 0, upper = 1> p; //mixture probability
}
model{
  real alpha_invalid = 0.75;
  real beta_invalid = 0.75;  
  alpha1 ~ uniform(1, 1000);
  beta1 ~ uniform(1, 1000);
  for (n in 1:N) {
    target += log_sum_exp(log(p) + beta_lpdf(si[n] | alpha_invalid, beta_invalid),
                          log(1 - p) + scenario1a_lpdf(si[n] | max_shed, alpha1, beta1, alpha2, beta2));
  }
}
