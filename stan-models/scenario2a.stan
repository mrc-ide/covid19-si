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
    // s is the time of infection
    inf_density = gamma_lpdf(0.5|alpha1, beta1) + log(0.5);
    inc_density = gamma_lpdf(x - 0.5|alpha2, beta2) + log(0.5);
    out =  exp(inf_density + inc_density);
    s = 0.5;
    while(s <= x) {
      s = s + 0.5;
      inf_density = gamma_lpdf(s|alpha1, beta1) + log(0.5);
      inc_density = gamma_lpdf(x - s|alpha2, beta2) + log(0.5);
      out = out + exp(inf_density + inc_density);
    }
    out = log(out);
    // Now do -log(F(nu))
    out = out - gamma_lcdf(nu|alpha1, beta1);
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
  real <lower = 0> mu; // infectious profile parameter
  real <lower = 0> cv; // // infectious profile parameter
} transformed parameters {
  real<lower=0> alpha1;
  real<lower=0> beta1;
  alpha1 = 1 / (cv^2);
  beta1 = 1 / (mu * cv^2);
}  
model{
  mu ~ uniform(1, 40);
  cv ~ uniform(0, 2);
  for (n in 1:N) {
    si[n] ~ scenario2a(nu[n], max_shed, alpha1, beta1, alpha2, beta2);
  }
}
