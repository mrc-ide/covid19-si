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
  


    real scenario3a_lpdf(real x,
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

    return out;
  }
  
  real scenario3a_v2_lpdf(real x,
                       real max_shed,
                       real offset,
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
    inf_density = beta_lpdf(width/max_shed + offset|alpha1, beta1)+ log(width/(max_shed + offset));
    inc_density = gamma_lpdf(x - (offset - width)|alpha2, beta2) + log(width);
    out =  exp(inf_density + inc_density);
    s = -offset + (2 * width);
    while(s < ulim) {
      inf_density = beta_lpdf((s + offset)/(max_shed + offset)|alpha1, beta1) + log(width/(max_shed + offset));
      inc_density = gamma_lpdf(x - s|alpha2, beta2) + log(width);
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    return out;
}
}
