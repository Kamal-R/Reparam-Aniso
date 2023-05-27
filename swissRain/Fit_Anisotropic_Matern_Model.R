###
### Setup
###

### Load libraries (model)
library(geostatsp)
library(rstan)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(MASS)
library(hdrcde)
library(LogConcDEAD)

### Setup
source( "Helper_Functions.R" )
data( swissRain )
set.seed( 7586 ) 

### TRUE: Save the Stan models in the 'Results' folder
### FALSE: Do not save the Stan models
save_Stan_models <- TRUE


##########################################
###
### Isotropic Simulation
###
##########################################

phi = c(1, 60, 1, 0, 2) # (sigma, gamma, rho, 2*psi, nu)
somerc_coord <- swissRain@coords / 1000 # 100 x 2 matrix of coordinates
swissAltitude_vector <- raster::extract(swissAltitude, swissRain) / 1000 # altitude

Beta = c(2, 0.10) # regression coefficients
X <- cbind( rep(1,length(swissRain@data$ID)),
            swissAltitude_vector ) # design matrix
shape = 7 # Gamma distribution: shape parameter

### Construct anisotropic Matern covariance matrix
rstan::expose_stan_functions( "Anisotropic_Matern_Model_swissRain.stan" )
N = length(swissRain@data$ID)
Sigma = matrix(NA, N, N) 

## Off-diagonal entries
for ( ii in 1:(N-1) ) {
  for ( jj in (ii+1):N ) {

    Sigma[ii,jj] <- aniso_matern( phi[5], phi[1], phi[2],
                                  aniso_dist( phi[3]-1, 0.5 * phi[4], ii, jj,
                                              somerc_coord[,1], somerc_coord[,2] ) )
    Sigma[jj,ii] <- Sigma[ii,jj]

  }
}

### Diagonal entries 
for (kk in 1:N) {

  Sigma[kk,kk] <- phi[1]^2

}

# Simulate response variable
Us = MASS::mvrnorm(n = 1, mu = rep(0,N), Sigma = Sigma )
swissRain_y = rgamma( n = 100, shape = shape,
                      rate = shape / exp( X %*% Beta + Us ) )

### Prepare data for and run Stan model 
N = length( swissRain@data$rain )
data_iso <- list( N = N, K = 2, # Indices
                  y = swissRain_y, # simulated response variable
                  coords = cbind( somerc_coord[,1], somerc_coord[,2] ), # coordinates
                  swissAltitude = swissAltitude_vector, # spatial variate
                  sigma = sqrt(2.0), # std. dev. for prior on theta_1 and theta_2
                  nu = 2, # order of Matern (Bessel_K) function
                  lambda_S = 1/2, # rate parameter for exponential prior on sigma
                  lambda_Y = 50 ) # rate parameter for exponential prior on 1/gamma

swissRain.model <- stan_model( file = "Anisotropic_Matern_Model_swissRain.stan" ) 
swissRain_samples_iso <- sampling( object = swissRain.model, data = data_iso, refresh = 250, 
                                   chains = 4, cores = 4, iter = 3000, warmup = 2000,
                                   control = list( adapt_delta = 1-2.5e-8, max_treedepth = 12, 
                                                   metric = "diag_e" ) ) 

### Save Stan model 
if (save_Stan_model) {
  saveRDS( swissRain_samples_iso, file = "Results/Isotropic_Simulation_Fit.RDS" ) 
}




##########################################  
###
### Anisotropic Simulation
###
##########################################

### Simulate Matern random field
phi = c(1, 60, 4.5, 68 * pi/180, 2) # (sigma, gamma, rho, 2*psi, nu)
somerc_coord <- swissRain@coords / 1000 # 100 x 2 matrix of coordinates
swissAltitude_vector <- raster::extract(swissAltitude, swissRain) / 1000 # altitude

Beta = c(2, 0.10) # regression coefficients
X <- cbind( rep(1,length(swissRain@data$ID)),
            swissAltitude_vector ) # design matrix
shape = 7 # Gamma distribution: shape parameter


# Construct anisotropic Matern covariance matrix
rstan::expose_stan_functions( "Anisotropic_Matern_Model_swissRain.stan" )
N = length(swissRain@data$ID)
Sigma = matrix(NA, N, N) 

## Off-diagonal entries
for ( ii in 1:(N-1) ) {
  for ( jj in (ii+1):N ) {

    Sigma[ii,jj] <- aniso_matern( phi[5], phi[1], phi[2],
                                  aniso_dist( phi[3]-1, 0.5 * phi[4], ii, jj, 
                                              somerc_coord[,1], somerc_coord[,2] ) )
    Sigma[jj,ii] <- Sigma[ii,jj]

  }
}

### Diagonal entries 
for (kk in 1:N) {

  Sigma[kk,kk] <- phi[1]^2

}

# Simulate response variable
Us = MASS::mvrnorm(n = 1, mu = rep(0,N), Sigma = Sigma )
swissRain_y = rgamma( 100, shape = shape, rate = shape / exp( X %*% Beta + Us ) )

### Prepare data for and run Stan model
N = length( swissRain@data$rain )
data_aniso <- list( N = N, K = 2, # Indices
                    y = swissRain_y, # response variable
                    coords = somerc_coord, # matrix of spatial coordinates
                    swissAltitude = swissAltitude_vector, # spatial covariate
                    sigma = sqrt(2.0), # std. dev. of priors on theta_1 and theta_2
                    nu = 2, # order of the Matern (Bessel_K) function
                    lambda_S = 1/2, # rate parameter for exponential prior on sigma 
                    lambda_Y = 50 ) # rate parameter for exponential prior on 1/gamma 
                  
swissRain.model <- stan_model( file = "Anisotropic_Matern_Model_swissRain.stan") 
swissRain_samples_aniso <- sampling( object = swissRain.model, data = data_aniso, refresh = 250,
                                     chains = 4, cores = 4, iter = 3000, warmup = 2000,
                                     control = list( adapt_delta = 1-1e-2, max_treedepth = 12,
                                                     metric = "diag_e" ) )
  
### Save Stan model
if (save_Stan_model) {
  saveRDS( swissRain_samples_aniso, file = "Results/Anisotropic_Simulation_Fit.RDS" ) 
}




##########################################
###
### Switzerland Rainfall Data
###
##########################################

### Design matrix and response variable
somerc_coord <- swissRain@coords / 1000 # 100 x 2 matrix of coordinates
swissAltitude_vector <- raster::extract(swissAltitude, swissRain) / 1000 # altitude 
X <- cbind( rep(1,length(swissRain@data$ID)), swissAltitude_vector ) # design matrix
swissRain_y <- swissRain@data$rain # observed rainfall

### Prepare data for and run Stan model
N = length( swissRain@data$rain )
data_real <- list( N = N, K = 2, # Indices
              y = swissRain_y, # response variable
              coords = somerc_coord, # spatial coordinates
              swissAltitude = swissAltitude_vector,  # spatial covariate
              sigma = sqrt(2.0), # std. dev. of theta_1 and theta_2 priors
              nu = 2, # order of the Matern function
              lambda_S = 1/2, # rate parameter for exponential prior on sigma 
              lambda_Y = 50 ) # rate parameter for exponential prior on 1/gamma 

swissRain.model <- stan_model( file = "Anisotropic_Matern_Model_swissRain.stan")
swissRain_samples <- sampling( object = swissRain.model, data = data_real, thin = 1, refresh = 250, 
                               chains = 4, cores = 4, iter = 3000, warmup = 2000,
                               control = list( adapt_delta = 1-1e-4, max_treedepth = 12, 
                                               metric = "diag_e" ) )

### Save Stan model
if (save_Stan_model) {
  saveRDS( swissRain_samples, file = "Results/Switzerland_Rainfall_Fit.RDS" ) 
}

# format(Sys.time(), format="%H:%M:%S")