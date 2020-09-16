functions{
  real scenario2a_lpdf(real x,
                       real nu, 
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
    width = 0.05;
    // time of infection can be no larger than the
    // max shedding time
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // and no larger than isolation
    if(ulim > nu) ulim = nu;

    // s is the time of infection
    inf_density = beta_lpdf(width/max_shed|alpha1, beta1) + log(width);
    inc_density = gamma_lpdf(x - width|alpha2, beta2) + log(width);
    out =  exp(inf_density + inc_density);
    s = 2 * width;
    while(s < ulim) {
      inf_density = beta_lpdf(s/max_shed|alpha1, beta1) + log(width);
      inc_density = gamma_lpdf(x - s|alpha2, beta2) + log(width);
      out = out + exp(inf_density + inc_density);
      s = s + width;      
    }
    out = log(out);
    // Now do -log(F(nu))
    if (nu>max_shed) out = out - beta_lcdf(0.999999|alpha1, beta1);
    else out = out - beta_lcdf(nu/max_shed|alpha1, beta1);
    return out;
 }
}
data{
  int N; // number of data points
  real si[N]; // Serial Interval
  real nu[N]; // Time of isolation
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter  
}
parameters{
  real <lower = 0> alpha1; // infectious profile parameter
  real <lower = 0> beta1; // // infectious profile parameter
} model{
  alpha1 ~ uniform(1, 2000);
  beta1 ~ uniform(1, 2000);
  for (n in 1:N) {
    si[n] ~ scenario2a(nu[n], max_shed, alpha1, beta1, alpha2, beta2);
  }
}
