functions{
  real scenario4a_lpdf(real x, 
                       real nu,
                       real max_shed,
                       real offset,
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
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // and no larger than isolation
    if(ulim > nu) ulim = nu;
    // s is the time of infection of the 2ndry case relative to onset in the infector
    // in order to account for pre-symptomatic infection we need s to be negative
    // our beta/gamma distributions don't support this so we need to introduce an offset
    // the range of s in this model is -offset to max_shed
    // initially rather than 0.1, s = -offset + 0.1
    // so 0.1 in the line below is a simplification of -offset + 0.1 - offset
    inf_density = beta_lpdf(0.1/(max_shed + offset)|alpha1, beta1) + log(0.1) - log(max_shed + offset);
    inc_density = gamma_lpdf((x - (-offset + 0.1))|alpha2, beta2) + log(0.1);
    //print("inc_density= ", inc_density, "x = ", x, "offset= ", offset);
    out =  exp(inf_density + inc_density);
    s = -offset + 0.2;
    while(s < ulim) {
      inf_density = beta_lpdf((s + offset)/(max_shed + offset)|alpha1, beta1) + log(0.1) - log(max_shed + offset);
      inc_density = gamma_lpdf(x - s|alpha2, beta2) + log(0.1); //the incubation period is still the SI (x) - s
      out = out + exp(inf_density + inc_density);
      s = s + 0.1;
    }
    out = log(out);
    
    // Now do -log(F(nu))
    if(nu < max_shed) out = out - beta_lcdf((nu + offset)/(max_shed + offset)|alpha1, beta1);
    // else isolation has happened after max_shed
    // and therefore all of an individual's infectivity has been
    // captured. 


    return out;
  }
}
data{
  int N; // number of data points
  real si[N]; // Serial Interval 
  real nu[N]; // Time of isolation
  real max_shed;
  real offset;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter  
}
parameters{
  real <lower = 0, upper = 30> alpha1; // infectious profile parameter
  real <lower = 0, upper = 30> beta1;  // infectious profile parameter
}
model{
  //alpha1 ~ uniform(1, 2000);
  //beta1 ~ uniform(1, 2000);
  for (n in 1:N) {
    si[n] ~ scenario4a(nu[n], max_shed, offset, alpha1, beta1, alpha2, beta2);
  }
}
