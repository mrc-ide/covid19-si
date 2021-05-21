functions{
  real invalid_lpdf(real x, real min_si, real max_si) {
    real out;
    real y;
    //y = map_into_interval2(x, min_si, max_si, 0.01, 0.99);
    //out = beta_lpdf(y| alpha_invalid, beta_invalid);
    out = uniform_lpdf(x | min_si, max_si);
    return(out);
  }
  
  real nf_lpdf(real t, real a, real b, real c, real tmax) {
    real out;
    real growing;
    real falling;
    falling = b * exp(-a * (t - tmax));
    growing = a * exp(b * (t - tmax));
    out = c * (log(a + b) - log(falling + growing));
    return out;
  }

  // Set recall to 0 to get model without recall bias
  // Set nu to be greater than max_shed to get model
  // without conditioning on nu.
  real validnf_lpdf(real x, real nu, real max_shed, real a,
                    real b, real c, real tmax, real recall, 
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
      inf_density =  nf_lpdf(s | a, b, c, tmax);
      inc_density = gamma_lpdf(x - s|alpha2, beta2);
      out = out + (exp(inf_density + inc_density) * width);
      s = s + width;
    }
    out = log(out) - recall * fabs(x - nu);
    return out;
  }

  // Set recall to 0 to get model without recall bias
  // Set nu to be greater than max_shed to get model
  // without conditioning on nu.
  // Allow infection to happen after isolation
  // with some probability pleak
  real validnf_leaky_lpdf(real x, real nu, real max_shed, real a,
                    real b, real c, real tmax, real recall, 
                          real alpha2, real beta2, real pleak, real width) {
    // The smalltest time point at which infection can occur.
    real s = -20;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    
    if (x > max_shed) ulim = max_shed;
    else ulim = x;

    out = 0;
    while (s < ulim) {
      inf_density =  nf_lpdf(s | a, b, c, tmax);
      if (s > nu) {
        inf_density =  pleak * inf_density;
      }
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
  matrix leaky_pdf_matrix(real[] nu_vec, real[] si_vec, real max_shed, 
                    real a, real b, real c, real tmax, real recall, 
                          real alpha2, real beta2, real pleak, real width,
                    int first_valid_nu) {

    int num_nu = size(nu_vec);
    int num_si = size(si_vec);
    matrix[num_si, num_nu] pdf_mat;
    matrix[num_nu, num_si] pdf_mat_t;
    // fill the fist valid column
    for (row in 1:num_si) {
      pdf_mat[row, first_valid_nu] = exp(validnf_leaky_lpdf(si_vec[row]|
                                                            nu_vec[first_valid_nu],
                                                            max_shed,
                                                            a, b, c, tmax,
                                                            recall, alpha2,
                                                            beta2, pleak, width));

      
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
          pdf_mat[row, col] = exp(validnf_leaky_lpdf(si_vec[row]| nu_vec[col],
                                                     max_shed, a, b, c, tmax,
                                                     recall, alpha2, beta2, pleak, 
                                                     width));
        }
      }
    }
    return(pdf_mat);
  }
  
  
  // Assume that nu_vec is sorted so that nu_vec[i] <= nu_vec[i + 1]
  // for all i.
  // Similarly si_vec
  // nus are running across columns and SIs are running down rows
  matrix pdf_matrix(real[] nu_vec, real[] si_vec, real max_shed, 
                    real a, real b, real c, real tmax, real recall, 
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
                                                      a, b, c, tmax,
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
                                               max_shed, a, b, c, tmax,
                                               recall, alpha2, beta2,
                                               width));
        }
      }
    }
    return(pdf_mat);
  }

  // Accounting for possible left truncation of data
  // For a case imported at time min_shed, any infections prior to
  // min_shed will not be observed
  real validnf_with_left_bias_lpdf(real x, real nu, real max_shed, real min_shed,
                                   real a, real b, real c, real tmax, real recall, 
                                   real alpha2, real beta2, real width) {

    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real s = min_shed;
    
    if (x > max_shed) ulim = max_shed;
    else ulim = x;
    if(ulim > nu) ulim = nu;
    out = 0;
    while (s < ulim) {
      inf_density =  nf_lpdf(s | a, b, c, tmax);
      inc_density = gamma_lpdf(x - s|alpha2, beta2);
      out = out + (exp(inf_density + inc_density) * width);
      s = s + width;
    }
    out = log(out) - recall * fabs(x - nu);
    return(out);
  }

  // Returns an array where the first dimension indexes time to isolation,
  // second dimension indexes 
  // time to importation since symptom onset and third dimentions indexes SIs
  // Each cell contains the pdf on natural scale for given nu, rho, si.
  real[, ,] pdf_matrix_with_left_bias(real[] nu_vec, real[] si_vec, real[] min_shed_vec,
                                      real max_shed, real a, real b,
                                      real c, real tmax, real recall, 
                                      real alpha2, real beta2,
                                      real width) {

    int num_nu = size(nu_vec);
    int num_si = size(si_vec);
    int num_rho = size(min_shed_vec);
    real out[num_nu, num_rho, num_si];
    real nu;
    real rho;
    real si;
    
    for (nu_index in 1:num_nu) {
      nu = nu_vec[nu_index];
      for (rho_index in 1:num_rho) {
        rho = min_shed_vec[rho_index];
        for (si_index in 1:num_si) {
          si = si_vec[si_index];
          out[nu_index, rho_index, si_index] =
            exp(validnf_with_left_bias_lpdf(si| nu, max_shed, rho, a,
                                            b, c, tmax, recall, alpha2,
                                            beta2, width));
        }
      }
    }
    return(out);
  }
}
