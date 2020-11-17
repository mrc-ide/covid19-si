functions{
  real scenario3a_lpdf(real x,
                       real max_shed,
                       real offset1,
                       real alpha1,      
                       real beta1,   
                       real alpha2,    
                       real beta2,
                       real width) {   

    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // s is the time of infection
    out =  0;
    s = -offset1 + width;
    while(s < ulim) {
      inf_density = beta_lpdf((s + offset1)/(max_shed + offset1)|alpha1, beta1);
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
