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
    if(nu >= max_shed) out = out;
    else out = out - beta_lcdf(nu/max_shed|alpha1, beta1);
    return out;
 }
}
