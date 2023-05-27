functions {
  real aniso_dist ( real phi_3,              // anisotropic ratio minus one
                    real phi_4,              // anisotropic angle 
                    int ii, int jj,          // indices: ii^th and jj^th observations
                    vector coords_x,         // x coordinates of observation locations
                    vector coords_y ) {      // y coordinates of observation locations
    return sqrt( square( cos(phi_4) * (coords_x[ii] - coords_x[jj]) + 
                         sin(phi_4) * (coords_y[ii] - coords_y[jj]) ) * ( phi_3 + 1 ) + 
                 square( -sin(phi_4) * (coords_x[ii] - coords_x[jj]) + 
                         cos(phi_4) * (coords_y[ii] - coords_y[jj]) ) * inv( (phi_3 + 1) ) );
  }
     
    real aniso_matern_oneandhalfnu ( real phi_1,            // spatial standard deviation
                                     real phi_2,            // range 
                                     real dist ) {          // transformed distance
    return square( phi_1 ) * ( 1 + sqrt(3) * dist / phi_2 ) * exp( -sqrt(3) * dist / phi_2 ); 

    }
  
}


data {
  int<lower=1> N;                            // number of observations
  vector[N] y_obs;                           // response: mercury at each observation site
  matrix[N,2] coords;                        // Coordinates
  
  int<lower=1> K;                            // number of predictors (inc. intercept)
  matrix[N,K] X;                             // regression matrix 

  real sigma;                                // s.d. of normal prior on theta_3, theta_4
  real lambda_S;                             // parameter for prior on theta_1
  real lambda_Y;                             // parameter for prior on theta_2
}


parameters {
  // Source: Page 348 of Stan Reference Manual 2.17 (for Alpha)
  vector[N] Alpha;    
  vector[K] Beta;

  // (Re-parameterized) anisotropic Matern parameters
  real theta_1;                              // theta_1 = 2 * log(phi_1) + log(phi_2) (Note: labelled theta_3 in paper)                
  real theta_2;                              // theta_2 = 2 * log(phi_1) - log(phi_2) (Note: labelled theta_4 in paper)
  real theta_3;                              // Cartesian coordinate: ...cos(2*psi) (Note: labelled theta_2 in paper)
  real theta_4;                              // Cartesian coordinate: ...sin(2*psi) (Note: labelled theta_1 in paper)
  
  // Hyper-parameters
  real<lower=0> sqrt_inv_shape;               // Gamma distribution: Coefficient of variation
}


transformed parameters {
  matrix[N,N] Sigma;                                   // Matern covariance matrix 
  real<lower=0> shape = inv(square(sqrt_inv_shape));   // Gamma distribution: shape parameter 

  // Natural parameters ( phi_1 = sigma, phi_2 = gamma, phi_3 = rho - 1, phi_4 = 2 * psi )
  real<lower=0> phi_1 = exp( ( theta_1 + theta_2 ) / 4 );           // phi_1 ~ Exp(lambda_S)
  real<lower=0> phi_2 = exp( ( theta_1 - theta_2 ) / 2 );           // 1/phi_2 ~ Exp(lambda_Y)
  real<lower=0> phi_3 = square( theta_3 ) + square( theta_4 );      // phi_3 = rho - 1 ~ Exp(1/lambda_1^2)
  real<lower=-pi(),upper=pi()> phi_4 = atan2( theta_4, theta_3 );   // phi_4 = 2 * psi ~ Unif(-pi,pi]
    
  // Define Matern covariance matrix Sigma
  for (ii in 1:(N-1)) {
    for (jj in (ii+1):N) {
      Sigma[ii,jj] = aniso_matern_oneandhalfnu( phi_1, phi_2, 
                                                aniso_dist( phi_3, 0.5 * phi_4, ii, jj, 
                                                            coords[,1], coords[,2] ) );
      Sigma[jj,ii] = Sigma[ii,jj];       
    }
  }
  
  for (kk in 1:N) {                        
    Sigma[kk,kk] = square(phi_1); 
  }

} 


model {   
  // Priors
  sqrt_inv_shape ~ exponential(0.5);           // Gamma distribution: Coefficient of variation
  
  // Prior on theta_1 and theta_2
  target += log(lambda_Y) + log(lambda_S) - 
            2 * ( theta_1 - theta_2 ) / 2 -
            lambda_Y * inv( exp( ( theta_1 - theta_2 ) / 2 ) ) - 
            lambda_S * exp( ( theta_1 + theta_2 ) / 4 ); 
  // Jacobian adjustment on theta_1 and theta_2
  target += -log(4) +                    
             ( theta_1 + theta_2 ) / 4 + 
             ( theta_1 - theta_2 ) / 2;
  // Priors on theta_3 and theta_4           
  theta_3 ~ normal(0,sigma);                  // theta_3: sigma = sqrt(2) => rho-1 ~ Exp(1/4)
  theta_4 ~ normal(0,sigma);                  // theta_4:             and => phi_A ~ Unif[-pi/2,pi/2)
  
  // Regression parameters 
  Alpha ~ normal(0,1);                        // => U_s ~ MVN(0, Sigma)
  Beta ~ normal(0,2);                
  
  
  // Likelihood (Gamma distribution)
  y_obs ~ gamma( shape,                       // log-linked gamma response    
                 rep_vector(shape,N) ./ exp(X * Beta + cholesky_decompose(Sigma) * Alpha) );

}


generated quantities { 
  // Legacy variables (not necessary)
  real phi_X = phi_2 * ( phi_3 + 1 );
  real phi_Y = phi_2;  
}
