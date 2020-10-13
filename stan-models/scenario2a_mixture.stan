functions{
   real scenario2a_lpdf(real x,
                        real nu, 
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
    // time of infection can be no larger than the
    // max shedding time
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // and no larger than isolation
    if(ulim > nu) ulim = nu;

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
    // Now do -log(F(nu))
    if(max_shed < nu) out = out - beta_lcdf(nu/max_shed|alpha1, beta1);

    return out;
 }
}
data{
  int N; // number of data points
  real si[N];
  real nu[N];
  real max_shed;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter
  real <lower = 0> max_si;  
  real <lower = 0> alpha_invalid;
  real <lower = 0> beta_invalid;
  real width;
}
parameters{
  simplex[2] theta;   
  real <lower = 1, upper = 50> alpha1; // infectious profile parameter
  real <lower = 1, upper = 50> beta1;  // infectious profile parameter
}
model{
  vector[2] log_theta = log(theta);  
  for (n in 1:N) {
    vector[2] lps = log_theta;
    lps[1] = lps[1] + beta_lpdf(si[n]/max_si | alpha_invalid, beta_invalid);
    lps[2] = lps[2] + scenario2a_lpdf(si[n] | nu[n], max_shed, alpha1, beta1, alpha2, beta2, width);
    target += log_sum_exp(lps);
  }
}
