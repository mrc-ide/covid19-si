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
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // s is the time of infection
    inf_density = beta_lpdf(0.1/max_shed|alpha1, beta1) + log(0.1) - log(max_shed);
    inc_density = gamma_lpdf(x - 0.1|alpha2, beta2) + log(0.1);
    out =  exp(inf_density + inc_density);
    s = 0.2;
    while(s < ulim) {
      inf_density = beta_lpdf(s/max_shed|alpha1, beta1) + log(0.1) - log(max_shed);
      inc_density = gamma_lpdf(x - s|alpha2, beta2)+ log(0.1);
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
  real <lower = 0> alpha2; // incubation period parameter
  real <lower = 0> beta2; // incubation period parameter  
}
parameters{

  
  real <lower = 0, upper = 100> alpha1; // infectious profile parameter
  real <lower = 0, upper = 100> beta1;  // infectious profile parameter

}
model{
  //alpha1 ~ uniform(1, 10);
  //beta1 ~ uniform(1, 10);
  for (n in 1:N) {
    si[n] ~ scenario1a(max_shed, alpha1, beta1, alpha2, beta2);
  }
}
