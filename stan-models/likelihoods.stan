functions{
  real invalid_lpdf(real x, real min_si, real max_si) {
    real out;
    out = uniform_lpdf(x | min_si, max_si);
    return(out);
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
  real valid_beta_lpdf(real x, real nu, real max_shed, real offset1,
                      real recall, real alpha1, real beta1,
                      real alpha2, real beta2, real width) {
    real s;
    real out;
    real inf_density;
    real inc_density;
    real ulim;
    real max_shed_shifted;
    real nu_shifted;
    
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
    out = log(out) - recall * fabs(x - nu);
    return out;
  }  
  
  // Assume that nu_vec is sorted so that nu_vec[i] <= nu_vec[i + 1]
  // for all i.
  // Similarly si_vec
  // nus are running across columns and SIs are running down rows
  matrix pdf_matrix_recall(real[] nu_vec, real[] si_vec, real max_shed, 
                           real offset1, real recall,
                           real alpha1, real beta1, 
                           real alpha2, real beta2, int first_valid_nu, real width) {

    int num_nu = size(nu_vec);
    int num_si = size(si_vec);
    matrix[num_si, num_nu] pdf_mat;
    matrix[num_nu, num_si] pdf_mat_t;
    real max_shed_shifted = max_shed - offset1;
    real nu_shifted;

    // fill the fist valid column
    for (row in 1:num_si) {
      pdf_mat[row, first_valid_nu] = exp(valid_beta_lpdf(si_vec[row]| nu_vec[first_valid_nu], max_shed, 
                                                     offset1, recall, alpha1, beta1, alpha2,
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
          pdf_mat[row, col] = exp(valid_beta_lpdf(si_vec[row]| nu_vec[col], max_shed, 
                                              offset1, recall, alpha1, beta1, alpha2,
                                                  beta2, width));
      
        }
      }
    }
    return(pdf_mat);
  }
}
