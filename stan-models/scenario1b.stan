functions{
  real scenario1b_lpdf(real x,
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
    return out;    
}
}
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
