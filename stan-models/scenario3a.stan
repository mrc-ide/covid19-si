functions{
  real scenario3a_lpdf(real x,
                       real max_shed,
                       real offset1,
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
    // s is the time of infection of the 2ndry case relative to onset in the infector
    // in order to account for pre-symptomatic infection we need s to be negative
    // our beta/gamma distributions don't support this so we need to introduce an offset1
    // the range of s in this model is -offset1 to max_shed
    // initially rather than 0.1, s = -offset1 + 0.1
    // so 0.1 in the line below is a simplification of -offset1 + 0.1 - offset1
    inf_density = beta_lpdf(0.1/(max_shed + offset1)|alpha1, beta1);
    inc_density = gamma_lpdf((x - (-offset1 + 0.1))|alpha2, beta2);
    //print("inc_density= ", inc_density, "x = ", x, "offset1= ", offset1);
    out =  exp(inf_density + inc_density);
    s = -offset1 + 0.2;
    while(s < ulim) {
      inf_density = beta_lpdf((s + offset1)/(max_shed + offset1)|alpha1, beta1);
      inc_density = gamma_lpdf(x - s|alpha2, beta2); //the incubation period is still the SI (x) - s
      out = out + exp(inf_density + inc_density);
      s = s + 0.1;
    }
    out = log(out);

    return out;
  }
}
data{
  int N; // number of data points
  real si[N];  
  real max_shed;
  real offset1;
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter  
}
parameters{
  real <lower = 1, upper = 30> alpha1; // infectious profile parameter
  real <lower = 1, upper = 30> beta1;  // infectious profile parameter
}
model{
  //alpha1 ~ uniform(1, 10);
  //beta1 ~ uniform(1, 10);
  for (n in 1:N) {
    si[n] ~ scenario3a(max_shed, offset1, alpha1, beta1, alpha2, beta2);
  }
}
