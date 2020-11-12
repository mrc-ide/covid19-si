functions{
 real scenario1a_lpdf(real x,
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
    // width should be smaller than the smallest SI otherwise
    // x - width will be -ve and gamma_lpdf will be 0, messing up
    // everything.
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
    if(nu < max_shed) out = out - beta_lcdf(nu/max_shed|alpha1, beta1);
    // else isolation has happened after max_shed
    // and therefore all of an individual's infectivity has been
    // captured. 
    return out;
 }
  // recall is the exponent on recall bias
  // alpha1 and beta 1 parameters of the infectious profile
  real scenario2a_with_recall_lpdf(real x,
                                   real nu, 
                                   real max_shed,
                                   real alpha1,      
                                   real beta1,
                                   real recall,
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
    if(nu < max_shed) out = out - beta_lcdf(nu/max_shed|alpha1, beta1);
    // else isolation has happened after max_shed
    // and therefore all of an individual's infectivity has been
    // captured.
    // Finally add contribution of recall bias
    out = out - recall * fabs(x - nu);
    return out;
 }

  
  real invalid_lpdf(real x, real max_si, real min_si, real alpha_invalid, real beta_invalid) {
    real out;
    real y;
    //y = (x + fabs(min_si))/ (max_si - min_si);
    // Mapping (min_si, max_si) to (0, 1)
    real a;
    real b;
    a = -min_si / (max_si - min_si);
    b = 1 / (max_si - min_si);
    y = a + b * x;
    out = beta_lpdf(y| alpha_invalid, beta_invalid);
    return(out);
  }
  
  
  real scenario3a_nowidth_lpdf(real x,
                       real max_shed,
                       real offset1,
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
    // width should be smaller than the smallest SI otherwise
    // x - width will be -ve and gamma_lpdf will be 0, messing up
    // everything.
    if(x > max_shed) ulim = max_shed;
    else ulim = x;
    // s is the time of infection
    inf_density = beta_lpdf((offset1 + width)/(max_shed + offset1)|alpha1, beta1);//+ log(width/(max_shed + offset1));
    inc_density = gamma_lpdf(x - (offset1 - width)|alpha2, beta2);// + log(width);
    out =  exp(inf_density + inc_density);
    s = -offset1 + (2 * width);
    while(s < ulim) {
      inf_density = beta_lpdf((s + offset1)/(max_shed + offset1)|alpha1, beta1);// + log(width/(max_shed + offset1));
      inc_density = gamma_lpdf(x - s|alpha2, beta2);// + log(width);
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    return out;
}

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
