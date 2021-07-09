functions{
  real invalid_lpdf(real x, real min_si, real max_si) {
    real out;
    real y;
    //y = map_into_interval2(x, min_si, max_si, 0.01, 0.99);
    //out = beta_lpdf(y| alpha_invalid, beta_invalid);
    out = uniform_lpdf(x | min_si, max_si);
    return(out);
  }
  
  real nf_lpdf(real t, real a, real b, real c) {
    // a is what stanc calls xi
    // b is what stan calls omega
    // c is what stan calls alpha
    real out = skew_normal_lpdf(t|a, b, c);
    return out;
  }

  // Set recall to 0 to get model without recall bias
  // Set nu to be greater than max_shed to get model
  // without conditioning on nu.
  real validnf_lpdf(real x, real nu, real max_shed, real a,
                    real b, real c, real recall, 
                    real alpha2, real beta2, real width) {
    // The smalltest time point at which infection can occur.
    real s = -20;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    
    if (x > max_shed) ulim = max_shed;
    else ulim = x;
    if(ulim > nu) ulim = nu;
    out = 0;
    while (s < ulim) {
      inf_density =  nf_lpdf(s | a, b, c);
      inc_density = gamma_lpdf(x - s|alpha2, beta2);
      out = out + (exp(inf_density + inc_density) * width);
      s = s + width;
    }
    out = log(out) - recall * fabs(x - nu);
    return out;
  }

  // Assume that nu_vec is sorted so that nu_vec[i] <= nu_vec[i + 1]
  // for all i.
  // Similarly si_vec
  // nus are running across columns and SIs are running down rows
  matrix pdf_matrix(real[] nu_vec, real[] si_vec, real max_shed, 
                    real a, real b, real c, real recall, 
                    real alpha2, real beta2, real width,
                    int first_valid_nu) {

    int num_nu = size(nu_vec);
    int num_si = size(si_vec);
    matrix[num_si, num_nu] pdf_mat;
    matrix[num_nu, num_si] pdf_mat_t;
    // fill the fist valid column
    for (row in 1:num_si) {
      pdf_mat[row, first_valid_nu] = exp(validnf_lpdf(si_vec[row]|
                                                      nu_vec[first_valid_nu],
                                                      max_shed,
                                                      a, b, c, 
                                                      recall, alpha2,
                                                      beta2, width));

      
    }
    for (row in 1:num_si) {
      for (col in (first_valid_nu + 1):num_nu) {
        // First check of the nu here is greater than the SI
        if (nu_vec[col] > si_vec[row]) {
          // then copy the value from a previously calculated cell
          // in the same row but from an earlier column
          pdf_mat[row, col] = pdf_mat[row, col - 1];
        } else {
          // Fill in the basic pdf as each cell will have to conditioned
          // on nu separately.
          pdf_mat[row, col] = exp(validnf_lpdf(si_vec[row]| nu_vec[col],
                                               max_shed, a, b, c, 
                                               recall, alpha2, beta2,
                                               width));
        }
      }
    }
    return(pdf_mat);
  }
}
