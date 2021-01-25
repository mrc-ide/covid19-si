functions{
  // maps a point x in the interval (x_1, y_1) into (x_2, y_2)
  // useful when we use beta density below
  // We map the interval (offset, max_shed) into (0.01, 0.99)
  // 0.01 and 0.99 because mapping offset to 0 gives Inf as log-likelihood
  // Example: map_into_interval(-1, -1, 21, 0, 1)
  // map_into_interval(0, 0, 1, -1, 21)
  real map_into_interval2(real x, real x_1, real y_1, real x_2,
                         real y_2) {
    real slope;
    real intercept;
    real y;
    slope = (x_2 - y_2) / (x_1 - y_1);
    intercept = x_2 - x_1 * slope;
    y = slope * x + intercept;
    return y;
  }
  
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
    out =  0;
    s = width;
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

  
  real invalid_lpdf(real x, real min_si, real max_si, 
                    real alpha_invalid, real beta_invalid) {
    real out;
    real y;
    y = map_into_interval2(x, min_si, max_si, 0.01, 0.99);
    out = beta_lpdf(y| alpha_invalid, beta_invalid);
    return(out);
  }
  
  // Make sure to pass offset1 as a *negative* number
  real scenario3a_lpdf(real x, real max_shed, real offset1, 
                       real alpha1, real beta1, real alpha2, real beta2,
                       real width) {

    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real max_shed_shifted;
    real s_remapped;
    
    if (x > max_shed) ulim = max_shed;
    else ulim = x;
    out = 0;
    // the whole infectious profile is shifted
    // right if offset < 0  and left if offset > 0
    max_shed_shifted = max_shed - offset1;
    //print("max shed shifted = ", max_shed_shifted);
    // Mapping s which varies from -offset to ulim to
    // interval (0, 1). We want to integrate from offset to ulim
    s = offset1;
    // Now map it into (0, 1)
    //s = (s - offset1) / max_shed_shifted;
    while (s <= ulim) {
      //print("s = ", s);
      s_remapped = map_into_interval2(s, offset1, max_shed, 0.01, 0.99);
      inf_density = beta_lpdf(s_remapped|alpha1, beta1);
      //print("With new mapping function ", inf_density);
      //inf_density = beta_lpdf((s - offset1)/max_shed_shifted |alpha1, beta1);
      //print("With old mapping function ", inf_density);      
      inc_density = gamma_lpdf(x - s|alpha2, beta2);
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    return out;
  }
  
  // Make sure to pass offset1 as a *negative* number
  real scenario4a_lpdf(real x, real nu, real max_shed, real offset1, 
                       real alpha1, real beta1, real alpha2, real beta2,
                       real min_si, real max_si, real width) {

    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real max_shed_shifted;
    real nu_shifted;
    real s_remapped;    

    if (x > max_shed) ulim = max_shed;
    else ulim = x;
    if(ulim > nu) ulim = nu;    
    
    nu_shifted = nu - offset1;
    max_shed_shifted = max_shed - offset1;
    s = offset1 + width;
    out = 0;
    while (s < ulim) {
      s_remapped = map_into_interval2(s, min_si, max_shed, 0.01, 0.99);
      inf_density = beta_lpdf(s_remapped|alpha1, beta1);
      inc_density = gamma_lpdf(x - s|alpha2, beta2);
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    if(nu < max_shed) {
      out = out - beta_lcdf(nu_shifted / max_shed_shifted|alpha1, beta1);
    }
    return out;
  }

  real s4_normalising_constant(real[] y_vec, real nu, real max_shed, 
                               real offset1, real alpha1, real beta1,
                               real alpha2, real beta2, real min_si,
                               real max_si, real width) {

    real denominator = 0;
    // Smallest SI allowed is should be at least offset1 + width
    // But in fact when SI is offset1 + width, pdf is -Inf
    // So make it a tiny bit bigger
    //real s = offset1 + width + 0.001;
    // Add in natural scale, and then take lof because we want
    // log of integral
    // This should be multipied with width but to speed things up I
    // have set the width to 1.
    for (y in y_vec) {
      denominator +=
        exp(scenario4a_lpdf(y| nu, max_shed, offset1, alpha1, beta1,
                            alpha2, beta2, min_si, max_si, width));
    }
    denominator = log(denominator);
    return denominator;
  }

  real approx_normalising_constant(real x, real nu, real max_si, 
                                   real min_si, real recall,
                                   real width) {
    real si_inner = min_si;
    real denominator = 0;
    while (si_inner <= max_si) {
      denominator = denominator + exp(-recall * fabs(si_inner - nu));
      si_inner = si_inner + width; 
    }
    denominator = log(denominator);
    return denominator;
  }  
  // x : SI
  // nu: Delay from symptom onset to isolation
  // max_shed: Maximum possible time of infection of secondary case
  // offset1: Minimum possible time of infection of secondary case
  // Should be negative if time is before symptom onset
  // recall: coefficient of recall bias component
  // alpha 1, beta1 : shape1 and shape2 parameters of infectious profile
  // alpha 2, beta2 : shape and rate parameters of incubation profile
  // width: width of rectangle for aproximation of integral. Should be
  // smaller than the samllest SI.
  // min_si cannot be achieved without an incubation period of 0
  // so when x = min_si, this function will return -Inf
  real full_model_lpdf(real x, real nu, real max_shed, real offset1,
                       real recall, real alpha1, real beta1,
                       real alpha2, real beta2, real width,
                       real max_si, real min_si) {
    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real max_shed_shifted;
    real nu_shifted;
    real denominator;
    
    if (x > max_shed) ulim = max_shed;
    else ulim = x;
    if(ulim > nu) ulim = nu;
    out = 0;
    // the whole infectious profile is shifted
    // right if offset < 0  and left if offset > 0
    max_shed_shifted = max_shed - offset1;
    nu_shifted = nu - offset1;
    //print("nu shifted = ", nu_shifted);
    //print("max shed shifted = ", max_shed_shifted);
    // Mapping s which varies from -offset to ulim to
    // interval (0, 1). We want to integrate from offset to ulim
    s = offset1 + width;
    // Now map it into (0, 1)
    //s = (s - offset1) / max_shed_shifted;
    
    while (s < ulim) {
      inf_density = beta_lpdf((s - offset1)/max_shed_shifted |alpha1, beta1);
      inc_density = gamma_lpdf(x - s|alpha2, beta2); 
      out = out + exp(inf_density + inc_density);
      s = s + width;
    }
    out = log(out);
    if(nu < max_shed) {
      out = out - beta_lcdf(nu_shifted / max_shed_shifted|alpha1, beta1);
    }
    out = out - recall * fabs(x - nu);
    //denominator = normalising_constant(nu, max_si, min_si, recall);
    //print("denominator = ", denominator);
    //denominator = approx_normalising_constant(x, nu, max_si, min_si,
    //                                          recall, width);
    //out = out - denominator;
    //print("out = ", out);
    return out;
  }

  real normalising_constant(real[] y_vec, real nu, real max_shed, 
                            real offset1, real recall, real alpha1, 
                            real beta1, real alpha2, real beta2, 
                            real width, real max_si, real min_si) {

    real denominator = 0;
    // Smallest SI allowed is should be at least offset1 + width
    // But in fact when SI is offset1 + width, pdf is -Inf
    // So make it a tiny bit bigger
    //real s = offset1 + width + 0.001;
    // Add in natural scale, and then take log because we want
    // log of integral
    // This should be multipied with width but to speed things up I
    // have set the width to 1.
    for (y in y_vec) {
      denominator +=
        exp(full_model_lpdf(y| nu, max_shed, offset1, recall, alpha1,
                            beta1, alpha2, beta2, width, max_si, min_si));
      //print("denominator is now ", denominator);
      //print("y is now ", y);
      //print("s is now ", s);
    }
    //print("total denominator = ", denominator);
    denominator = log(denominator);
    return denominator;
  }
  
}
