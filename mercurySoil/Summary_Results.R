###########################
###
### Setup
###
###########################
library(dplyr)
RandomFields::RFoptions(install="no")




###########################
###
### Load cmdstanr files 
### (from the 'Results" folder
### in the working directory)
###
###########################
### This assigns to output_files the .csv file names the 
### cmdstanr output has been saved to. 
### It requires that 'Fit_Anisotropic_Matern_Model.Rmd' has been
### run and mercurySoil_samples is in the global environment.
output_files <- mercurySoil_samples$output_files()

### Others, set output_files by entering the file names of the .csv's output by 
### cmdstanr (saved to the 'Results' folder by the 'Fit_Anisotropic_Matern_Model.Rmd'. 
### For example (update to the file names in your 'Results' folder):
# output_files <- c( "Results/Anisotropic_Matern_Model_MercurySoil_Chain_1.csv",
#                    "Results/Anisotropic_Matern_Model_MercurySoil_Chain_2.csv",
#                    "Results/Anisotropic_Matern_Model_MercurySoil_Chain_3.csv",
#                    "Results/Anisotropic_Matern_Model_MercurySoil_Chain_4.csv" )

cmdstanfiles <- list()
for (f in output_files) {
  cmdstanfiles[[f]] <- data.table::fread( cmd = paste0("grep -v '^#' ", f) )
}

cmdstanfiles_allchains <- cmdstanfiles %>% data.table::rbindlist()
### Remove cmdstanfiles list after constructing cmdstanfiles_allchains to free system memory
rm(cmdstanfiles); gc()



###########################
###
### Summary Table
###
###########################

### Select variables 
Beta <- cmdstanfiles_allchains %>% dplyr::select( starts_with("Beta") ) 
phi <- cmdstanfiles_allchains %>% dplyr::select( starts_with("phi") ) 
theta <- cmdstanfiles_allchains %>% dplyr::select( starts_with("theta") )
sqrt_inv_alpha <- cmdstanfiles_allchains %>% dplyr::select( starts_with("sqrt_inv_shape") ) 

### Summary: Regression
Beta_summary <- c( round( c( mean( Beta$Beta.1 ), 
                             quantile( Beta$Beta.1, probs = c(0.025,0.50,0.975) ) ), 2 ),
                   round( c( mean( Beta$Beta.2 ), 
                             quantile( Beta$Beta.2, probs = c(0.025,0.50,0.975) ) ), 2 ),
                   round( c( mean( Beta$Beta.3 ), 
                             quantile( Beta$Beta.3, probs = c(0.025,0.50,0.975) ) ), 2 ),
                   round( c( mean( Beta$Beta.4 ), 
                             quantile( Beta$Beta.4, probs = c(0.025,0.50,0.975) ) ), 2 ) )

### Summary: Matern parameters
Phi_summary <- c( round( c( mean( phi$phi_1 ), quantile( phi$phi_1, probs = c(0.025,0.50,0.975) ) ), 2 ),
                  round( c( mean( phi$phi_2 ), quantile( phi$phi_2, probs = c(0.025,0.50,0.975) ) ), 2 ),
                  round( c( mean( phi$phi_3+1 ), 
                            bayestestR::hdi( phi$phi_3+1, ci = 0.95 )$CI_low,
                            quantile( phi$phi_3+1, probs = 0.50 ),
                            bayestestR::hdi( phi$phi_3+1, ci = 0.95 )$CI_high ), 2 ),
                  round( -c( mean( 180/pi * phi$phi_4 ) / 2, 
                             quantile( 180/pi * phi$phi_4 / 2, 
                                       probs = c(0.025,0.50,0.975) ) ), 2 ) )

### Summary: Transformed Matern parameters 
Theta_summary <- c( round( c( mean( theta$theta_1 ), 
                              quantile( theta$theta_1, probs = c(0.025,0.50,0.975) ) ), 2 ),
                    round( c( mean( theta$theta_2 ), 
                              quantile( theta$theta_2, probs = c(0.025,0.50,0.975) ) ), 2 ),
                    round( c( mean( theta$theta_3 ), 
                              quantile( theta$theta_3, probs = c(0.025,0.50,0.975) ) ), 2 ),  
                    round( -c( mean( theta$theta_4 ), 
                               quantile( theta$theta_4, probs = c(0.025,0.50,0.975) ) ), 2 ) )

### Summary: Axis-Aligned range parameters
phi_Y <- phi$phi_2 * sqrt( (phi$phi_3+1) ) # Scaling by > 1
phi_X <- phi$phi_2 / sqrt( (phi$phi_3+1) ) # Scaling by < 1 

phi_Y_summary <- round( c( mean( phi_Y ), 
                           quantile( phi_Y, probs = c(0.025,0.50,0.975) ) ), 2 )
phi_X_summary <- round( c( mean( phi_X ), 
                           quantile( phi_X, probs = c(0.025,0.50,0.975) ) ), 2 )


### Summary: Coefficient of variation
alpha_summary <- round( c( mean( sqrt_inv_alpha$sqrt_inv_shape ), 
                           quantile( sqrt_inv_alpha$sqrt_inv_shape, probs = c(0.025,0.50,0.975) ) ), 2 )


### Summary table
knitr::kable( t( tibble(  "$\\rho$" = Phi_summary[9:12], 
                          "$\\psi$" = Phi_summary[13:16],
                          "$\\theta_{1}$" = Theta_summary[13:16], 
                          "$\\theta_{2}$" = Theta_summary[9:12],
                          "$\\phi_Y$" = phi_Y_summary, 
                          "$\\phi_X$" = phi_X_summary ) ) ) %>% 
  kableExtra::kable_styling()

### Remove cmdstanfiles_allchains after extracting relevant variables to free system memory
rm(cmdstanfiles_allchains); gc()




###########################
###
### Bivariate Log-Concave
### Density Estimation
###
###########################

include_legend <- TRUE

### Load libraries 
library(rstan)
library(ggplot2)
library(latex2exp)
library(MASS)
library(hdrcde)
library(LogConcDEAD)



#######################################
### Density estimation of (theta_1,theta_2)
#######################################

theta_XY <- cbind( -theta$theta_4,
                   theta$theta_3 )
theta_XY_rounded <- round( theta_XY, digits = 2 )
theta_XY_weights <- getweights( theta_XY_rounded )
theta_XY_lcd <- mlelcd( x = theta_XY_weights$x,
                        w = theta_XY_weights$w )

### Obtain grid points
N = 100
h_theta <- c( MASS::bandwidth.nrd( -theta$theta_4 ),
              MASS::bandwidth.nrd( theta$theta_3 ) )
real_grid_points_theta <- MASS::kde2d( -theta$theta_4,
                                       theta$theta_3,
                                       h = h_theta, n = N )

### Evaluate density at grid points 
log_concave_density_theta <- expand.grid("x" = real_grid_points_theta$x,
                                         "y" = real_grid_points_theta$y)
theta_XY_theta_eval <- dlcd( x = cbind( log_concave_density_theta$x,
                                        log_concave_density_theta$y ),
                             lcd = theta_XY_lcd )


### Calculate credible region
log_concave_density_theta$fhat <- as.vector(theta_XY_theta_eval)
log_concave_density_theta$fhat_discretized <- log_concave_density_theta$fhat /
  sum(log_concave_density_theta$fhat)
credible_region_levels_theta <- ggdensity:::find_cutoff(
  df = log_concave_density_theta,
  conf = rev( c( 0.25, 0.50, 0.80, 0.90, 0.95, 0.99 ) ) )


### Plot bivariate posterior density of (theta_1,theta_2)
png( file = "Figures/mercurySoil_Polar_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

par( mar = c(5.1, 5.1, 2.1, 2.1) )
hdrcde:::hdrcde.filled.contour(
  x = real_grid_points_theta$x,
  y = real_grid_points_theta$y,
  z = matrix( theta_XY_theta_eval, N, N ),
  levels = c( credible_region_levels_theta, max(theta_XY_theta_eval) ),
  color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
  xlim = c(-2,2), ylim = c(-2,2), 
  plot.axes = { title( main = "", 
                       xlab = expression( theta[1] ), 
                       ylab = expression( theta[2] ),
                       cex.lab = 1.80 )
    axis(1, seq(from = -2.0, to = 2.0, by = 1.0), cex.axis = 1.40) 
    axis(2, seq(from = -2.0, to = 2.0, by = 1.0), cex.axis = 1.40) } )
    if ( include_legend ) {
      legend( x = 0.75, y = 1.75, bg = "white",
              legend = rev( scales::percent(
                c(0.25, 0.50, 0.80, 0.90, 0.95, 0.99),
                accuracy = 1 ) ),
              fill = rev( hcl.colors(6, "YlOrRd") ),
              title = "Credible Region", box.col = "white",
              ncol = 2, cex = 1.10 )
    }
par( mar = c(5.1, 4.1, 4.1, 2.1) )

dev.off()




#######################################
### Density estimation of (rho,psi)
#######################################

### Obtain grid points 
N = 100
h_phi <- c( MASS::bandwidth.nrd( -180/pi * phi$phi_4 / 2 ),
            MASS::bandwidth.nrd( phi$phi_3 ) )
real_grid_points_phi <- MASS::kde2d( -180/pi * phi$phi_4 / 2,
                                     phi$phi_3,
                                     h = h_phi, n = N )


### Transform grid points from phi-space to theta-space
x2 = list()
x2$x = matrix( NA, ncol = N, nrow = N )
x2$y = matrix( NA, ncol = N, nrow = N )

for (ii in 1:N) {
  for (jj in 1:N) {
    
    # theta_1 
    x2$x[ii,jj] <- sqrt( real_grid_points_phi$y[ii] ) *
      sin( 2 * pi/180 * real_grid_points_phi$x[jj] ) 
    # theta_2
    x2$y[ii,jj] <- sqrt( real_grid_points_phi$y[ii] ) *
      cos( 2 * pi/180 * real_grid_points_phi$x[jj] ) 
    
  }
}


### Evaluate density at grid points in theta-space
theta_XY_eval_on_theta_grid <- dlcd( x = cbind( as.vector(x2$x),
                                                as.vector(x2$y) ),
                                     lcd = theta_XY_lcd )


### Transform f_theta to f_phi
x3 = list()
x3$x = matrix( NA, ncol = N, nrow = N )
x3$y = matrix( NA, ncol = N, nrow = N )
x3$z = matrix( NA, ncol = N, nrow = N )

theta_XY_eval_on_theta_grid_matrix <- matrix( theta_XY_eval_on_theta_grid,
                                              nrow = N, ncol = N )

for (ii in 1:N) {
  for (jj in 1:N) {
    
    ### Transform from theta-space to phi-space
    x3$x[ii,jj] <- x2$x[ii,jj]^2 + x2$y[ii,jj]^2
    x3$y[ii,jj] <- 0.5 * atan2( x2$x[ii,jj], x2$y[ii,jj] )
    x3$z[ii,jj] <- theta_XY_eval_on_theta_grid_matrix[ii,jj] * 1
    
  }
}


### Calculate credible regions
log_concave_density_phi <- expand.grid("x" = x3$x, "y" = 180/pi * x3$y)
log_concave_density_phi$fhat <- as.vector(x3$z)
log_concave_density_phi$fhat_discretized <- log_concave_density_phi$fhat /
  sum(log_concave_density_phi$fhat)
credible_region_levels_phi <- ggdensity:::find_cutoff( 
  df = log_concave_density_phi,
  conf = rev( c( 0.25, 0.50, 0.80, 0.90, 0.95, 0.99 ) ) )


### Plot bivariate posterior density of (rho,psi)
png( file = "Figures/mercurySoil_Angle_Ratio_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

par( mar = c(5.1, 5.1, 2.1, 2.1) )
hdrcde:::hdrcde.filled.contour(
  x = x3$x[,1]+1, y = 180/pi * x3$y[1,], z = x3$z,
  levels = c( credible_region_levels_phi, max(x3$z) ),
  color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
  xlim = c(1,2.5), ylim = c(-90,90), 
  plot.axes = { title( main = "", 
                       xlab = expression( rho ), 
                       ylab = expression( psi~~"(Degrees)" ),
                       cex.lab = 1.80 )
    axis(1, c(1.0,1.5,2.0,2.5), cex.axis = 1.40) 
    axis(2, seq(from = -90, to = 90, by = 30), cex.axis = 1.40) } )
    if ( include_legend ) { 
      legend( x = 2.0, y = 82.5, bg = "white",
              legend = rev( scales::percent(
                c(0.25, 0.50, 0.80, 0.90, 0.95, 0.99),
                accuracy = 1 ) ),
              fill = rev( hcl.colors(6, "YlOrRd") ),
              title = "Credible Region", box.col = "white",
              ncol = 2, cex = 1.10 )
    }
par( mar = c(5.1, 4.1, 4.1, 2.1) )

dev.off()



#######################################
### Density estimation of (theta_3,theta_4)
#######################################

theta_XY <- cbind( theta$theta_1,
                   theta$theta_2 )
theta_XY_rounded <- round( theta_XY, digits = 2 )
theta_XY_weights <- getweights( theta_XY_rounded )
theta_XY_lcd <- mlelcd(x = theta_XY_weights$x, w = theta_XY_weights$w)



### Obtain grid points
N = 100
h_theta <- c( MASS::bandwidth.nrd( theta$theta_1 ),
              MASS::bandwidth.nrd( theta$theta_2 ) )
real_grid_points_theta <- MASS::kde2d( theta$theta_1,
                                       theta$theta_2,
                                       h = h_theta, n = N )

### Evaluate density at grid points 
log_concave_density_theta <- expand.grid("x" = real_grid_points_theta$x, 
                                         "y" = real_grid_points_theta$y) 
theta_XY_theta_eval <- dlcd( x = cbind( log_concave_density_theta$x, 
                                        log_concave_density_theta$y ),
                             lcd = theta_XY_lcd )


### Calculate credible regions
log_concave_density_theta$fhat <- as.vector(theta_XY_theta_eval)
log_concave_density_theta$fhat_discretized <- log_concave_density_theta$fhat / 
  sum(log_concave_density_theta$fhat)
credible_region_levels_theta <- ggdensity:::find_cutoff( 
  df = log_concave_density_theta, 
  conf = rev( c( 0.25, 0.50, 0.80, 0.90, 0.95, 0.99 ) ) )

### Plot bivariate posterior density of (theta_3,theta_4)
png( file = "Figures/mercurySoil_TransIso_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

par( mar = c(5.1, 5.1, 2.1, 2.1) )
hdrcde:::hdrcde.filled.contour(
  x = real_grid_points_theta$x, 
  y = real_grid_points_theta$y, 
  z = matrix( theta_XY_theta_eval, N, N ),
  levels = c( credible_region_levels_theta, max(theta_XY_theta_eval) ),
  color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE), 
  xlim = c(3.0,6.0), ylim = c(-7.0,-5.0), 
  plot.axes = { title( main = "", 
                       xlab = expression( theta[3] ), 
                       ylab = expression( theta[4] ),
                       cex.lab = 1.80 )
    axis(1, seq(from = 3.0, to = 6.0, by = 1.0), cex.axis = 1.40)
    axis(2, seq(from = -7.0, to = -5.0, by = 1.0), cex.axis = 1.40) } )
    if ( include_legend ) {
      legend( x = 4.0, y = -4.5, bg = "white",
              legend = rev( scales::percent(
                c(0.25, 0.50, 0.80, 0.90, 0.95, 0.99),
                accuracy = 1 ) ),
              fill = rev( hcl.colors(6, "YlOrRd") ),
              title = "Credible Region", box.col = "white",
              ncol = 2, cex = 1.10 )
    }
par( mar = c(5.1, 4.1, 4.1, 2.1) )

dev.off()




#######################################
### Density estimation of (sigma,gamma)
#######################################
N = 100
h_phi <- c( MASS::bandwidth.nrd( phi$phi_1 ),
            MASS::bandwidth.nrd( phi$phi_2 ) )
real_grid_points_phi <- MASS::kde2d( phi$phi_1, # sigma 
                                     phi$phi_2, # gamma
                                     h = h_phi, n = N )

# Transform grid points from phi space to theta space
x2 = list()
x2$x = matrix( NA, ncol = N, nrow = N )
x2$y = matrix( NA, ncol = N, nrow = N )


for (ii in 1:N) {
  for (jj in 1:N) {
    
    # theta_3
    x2$x[ii,jj] <- 2 * log( real_grid_points_phi$x[ii] ) + 
      log( real_grid_points_phi$y[jj] ) 
    # theta_4
    x2$y[ii,jj] <- 2 * log( real_grid_points_phi$x[ii] ) - 
      log( real_grid_points_phi$y[jj] )
    
  }
}


### Evaluate density at the theta grid points
theta_XY_eval_on_theta_grid <- dlcd( x = cbind( as.vector(x2$x), 
                                                as.vector(x2$y) ),
                                     lcd = theta_XY_lcd )


### Transform f_theta to f_phi
x3 = list()
x3$x = matrix( NA, ncol = N, nrow = N )
x3$y = matrix( NA, ncol = N, nrow = N )
x3$z = matrix( NA, ncol = N, nrow = N )

theta_XY_eval_on_theta_grid_matrix <- matrix( theta_XY_eval_on_theta_grid, 
                                              nrow = N, ncol = N )

for (ii in 1:N) {
  for (jj in 1:N) {
    
    ### Transform from theta space to phi space 
    x3$x[ii,jj] <- exp( ( x2$x[ii,jj] + x2$y[ii,jj] ) / 4 ) # sigma 
    x3$y[ii,jj] <- exp( ( x2$x[ii,jj] - x2$y[ii,jj] ) / 2 ) # gamma
    x3$z[ii,jj] <- theta_XY_eval_on_theta_grid_matrix[ii,jj] * ( 4 / ( x3$x[ii,jj] * x3$x[ii,jj] ) )
    
  }
}


### Calculate credible regions
log_concave_density_phi <- expand.grid("x" = x3$x, "y" = x3$y) 
log_concave_density_phi$fhat <- as.vector(x3$z)
log_concave_density_phi$fhat_discretized <- log_concave_density_phi$fhat / 
  sum(log_concave_density_phi$fhat)

credible_region_levels_phi <- ggdensity:::find_cutoff( 
  df = log_concave_density_phi,
  conf = rev( c( 0.25, 0.50, 0.80, 0.90, 0.95, 0.99 ) ) )


png( file = "Figures/mercurySoil_Sigma_Gamma_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

par( mar = c(5.1, 5.1, 2.1, 2.1) )
hdrcde:::hdrcde.filled.contour(
  x = x3$x[,1], y = x3$y[1,], z = x3$z,
  levels = c( credible_region_levels_phi, max(x3$z) ),
  col = rev( hcl.colors(6, "YlOrRd") ),
  xlim = c(0.40,1), ylim = c(50,350), 
  plot.axes = { title( main = "", 
                       xlab = expression( sigma ), 
                       ylab = expression( gamma ),
                       cex.lab = 1.80 )
    axis(1, seq(from = 0, to = 1.10, by = 0.15), cex.axis = 1.40)
    axis(2, seq(from = 0, to = 300, by = 50), cex.axis = 1.40) } )
    if ( include_legend ) { 
      legend( x = 2.10, y = 180, bg = "white",
              legend = rev( scales::percent(
                c(0.25, 0.50, 0.80, 0.90, 0.95, 0.99),
                accuracy = 1 ) ),
              fill = rev( hcl.colors(6, "YlOrRd") ),
              title = "Credible Region", box.col = "white",
              ncol = 2, cex = 1.10 )
    }
par( mar = c(5.1, 4.1, 4.1, 2.1) )

dev.off()





####################################
###
### Process Sigma and Alpha
###
####################################

### Load posterior samples from cmdstanr
### Note: output_files is defined on lines 20-27
### (Reloading cmdstanfiles to extract Sigma and Alpha. These
### objects are not needed until this point in the script)
cmdstanfiles <- list()
for (f in output_files) {
  cmdstanfiles[[f]] <- data.table::fread( cmd = paste0("grep -v '^#' ", f) )
}


####################################
### Sigma
####################################

Sigma_C1 <- cmdstanfiles[[1]] %>% dplyr::select( starts_with("Sigma") )
Sigma_C2 <- cmdstanfiles[[2]] %>% dplyr::select( starts_with("Sigma") )
Sigma_C3 <- cmdstanfiles[[3]] %>% dplyr::select( starts_with("Sigma") )
Sigma_C4 <- cmdstanfiles[[4]] %>% dplyr::select( starts_with("Sigma") )

Sigma_C1_list <- vector( "list", length = dim(Sigma_C1)[1] )
Sigma_C2_list <- vector( "list", length = dim(Sigma_C2)[1] )
Sigma_C3_list <- vector( "list", length = dim(Sigma_C3)[1] )
Sigma_C4_list <- vector( "list", length = dim(Sigma_C4)[1] )

for ( ii in 1:dim(Sigma_C1)[1] ) {
  
  ### Note: Sigma is imported as N x N length vector -- convert to an N x N matrix.
  Sigma_C1_list[[ii]] <- matrix( unlist( unname( Sigma_C1[ii,] ) ),
                                 nrow = sqrt( dim(Sigma_C1)[2] ),
                                 ncol = sqrt( dim(Sigma_C1)[2] ) )
  
  
  Sigma_C2_list[[ii]] <- matrix( unlist( unname( Sigma_C2[ii,] ) ),
                                 nrow = sqrt( dim(Sigma_C2)[2] ),
                                 ncol = sqrt( dim(Sigma_C2)[2] ) )
  
  
  Sigma_C3_list[[ii]] <- matrix( unlist( unname( Sigma_C3[ii,] ) ),
                                 nrow = sqrt( dim(Sigma_C3)[2] ),
                                 ncol = sqrt( dim(Sigma_C3)[2] ) ) 
  
  
  Sigma_C4_list[[ii]] <- matrix( unlist( unname( Sigma_C4[ii,] ) ),
                                 nrow = sqrt( dim(Sigma_C4)[2] ),
                                 ncol = sqrt( dim(Sigma_C4)[2] ) ) 
  
  # if ( ii %% 25 == 0 ) { print( paste0("Iteration ", ii, " is complete.") ) }
  
}

Sigma <- c( Sigma_C1_list, Sigma_C2_list, Sigma_C3_list, Sigma_C4_list )



####################################
### Alpha
####################################

Alpha_C1 <- cmdstanfiles[[1]] %>% dplyr::select( starts_with("Alpha") )
Alpha_C2 <- cmdstanfiles[[2]] %>% dplyr::select( starts_with("Alpha") )
Alpha_C3 <- cmdstanfiles[[3]] %>% dplyr::select( starts_with("Alpha") )
Alpha_C4 <- cmdstanfiles[[4]] %>% dplyr::select( starts_with("Alpha") )

Alpha_C1_list <- vector( "list", length = dim(Alpha_C1)[1] )
Alpha_C2_list <- vector( "list", length = dim(Alpha_C2)[1] )
Alpha_C3_list <- vector( "list", length = dim(Alpha_C3)[1] )
Alpha_C4_list <- vector( "list", length = dim(Alpha_C4)[1] )

for ( ii in 1:dim(Alpha_C1)[1] ) {
  
  Alpha_C1_list[[ii]] <- Alpha_C1[ii,] 
  Alpha_C2_list[[ii]] <- Alpha_C2[ii,] 
  Alpha_C3_list[[ii]] <- Alpha_C3[ii,] 
  Alpha_C4_list[[ii]] <- Alpha_C4[ii,]
  
}

Alpha <- c( Alpha_C1_list, Alpha_C2_list, Alpha_C3_list, Alpha_C4_list )




####################################
###
### Process Sigma and Alpha
###
####################################
library('rgdal')
library('raster')
library('doParallel')
library('rnaturalearth')
library('sp')

sf::sf_use_s2(FALSE)
raster_cells <- 500
fixed_seed <- 2500



calculate_Us <- function(kk) { 
  
  # If Sigma is imported as an N x N vector (with N = 812)
  if ( length(Sigma[[kk]]) == 659344 ) { 
    
    # convert Sigma to an N x N matrix
    Sigma[[kk]] <- matrix( unlist( unname( Sigma[[kk]] ) ),
                           nrow = sqrt( length( Sigma[[kk]] ) ),
                           ncol = sqrt( length( Sigma[[kk]] ) ) ) 
  }
  
  return( chol( Sigma[[kk]] ) %*% unlist( Alpha[[kk]] ) )
  
}


####################################
### Load and process mercury soil data 
####################################
load(file = "Data/hgm.RData" ) 
rows_with_missing_values <- c( 2, 23, 24, 32, 35, 63, 401, 404, 457, 
                               566, 598, 607, 648, 672, 745, 771, 826 )

hgm@data <- hgm@data[-rows_with_missing_values,]
hgm@coords <- hgm@coords[-rows_with_missing_values,]
# dim(hgm@data); dim(hgm@coords)



####################################
### Setup
####################################
set.seed(fixed_seed) 
num_Us_samples <- 250
sample_idx <- sample(1:2000, 250) 

### Regopms to be plotted
theCountries = sp::spTransform( ne_countries( 
  country = c( "Norway", "France", "Sweden", 
               "Poland", "Austria", "Hungary", "Lithuania",
               "Latvia", "Estonia", "Germany",  "Greece", "Albania",
               "Croatia", "Switzerland", "Luxembourg", "Belgium", "Netherlands", 
               "Portugal", "Spain", "Ireland", "Italy", "Denmark", "United Kingdom",
               "Slovenia", "Finland", "Slovakia", "Czechia" ),
  type = 'countries', continent = "Europe"),
  projection(hgm), scale = "medium" )

### Regions to be plotted: approximate coordinates for labels
theCountries_sf = ne_countries( 
  country = c( "Norway", "France", "Sweden", 
               "Poland", "Austria", "Hungary", "Lithuania",
               "Latvia", "Estonia", "Germany",  "Greece", "Albania",
               "Croatia", "Switzerland", "Luxembourg", "Belgium", "Netherlands", 
               "Portugal", "Spain", "Ireland", "Italy", "Denmark", "United Kingdom",
               "Slovenia", "Finland", "Slovakia", "Czechia" ),
  type = 'countries', continent = "Europe", returnclass = "sf" )
temp <- sf::st_transform( sf::st_centroid(theCountries_sf)$geometry, crs = 3034 )
temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )



####################################
### Simulate random fields  
####################################
Us_conditional <- matrix( NA, ncol = num_Us_samples, nrow = 287500 ) 

for ( ii in 1:num_Us_samples ) {
  
  mercurySoil_stan = c( range = 1000 * phi$phi_2[sample_idx[ii]],
                        shape = rep(1.5,2*1000)[sample_idx[ii]],
                        variance = phi$phi_1[sample_idx[ii]]^2,
                        anisoAngleDegrees = -180/pi * phi$phi_4[sample_idx[ii]] / 2,
                        anisoRatio = phi$phi_3[sample_idx[ii]]+1 )

  hgm_Us <- hgm
  hgm_Us@bbox[1,1] <- hgm@bbox[1,1]-0.1*hgm@bbox[1,1]
  hgm_Us@data <- as.data.frame( calculate_Us(ii) )
  
  ### Simulate random fields using RFSimulate
  hgm_2 = RandomFields::RFsimulate(
    model = RandomFields::RMmatern( 
      var = mercurySoil_stan['variance'][[1]], 
      scale = mercurySoil_stan['range'][[1]], 
      nu = mercurySoil_stan['shape'][[1]], 
      Aniso = matrix( ncol = 2, c( sqrt( mercurySoil_stan['anisoRatio'][[1]] ), 0,
                                   0, 1 / sqrt( mercurySoil_stan['anisoRatio'][[1]] ) ) ) %*%
        matrix( ncol = 2, c( cos(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]]), sin(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]] ), 
                             -sin(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]]), cos(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]] ) ) ) ),
    x = raster::rasterToPoints( geostatsp::squareRaster(hgm_Us, raster_cells) )[,1],
    y = raster::rasterToPoints( geostatsp::squareRaster(hgm_Us, raster_cells) )[,2],
    data = hgm_Us, seed = fixed_seed )
  
  hgm_2 <- raster::raster( hgm_2, vals = hgm_2$V1,
                           nrows = geostatsp::squareRaster(hgm_Us, raster_cells)@nrows, 
                           ncol = geostatsp::squareRaster(hgm_Us, raster_cells)@ncols ) 
  
  
  ### Adjust location of country labels to fit space available on the screen (for the map)
  countries_to_label <- c(1,2,3,4,5,8,9,10,11,19,21,22,23,25)
  temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
  temp_matrix <- temp_matrix[c(1,2,3,4,5,8,9,10,11,19,21,22,23,25),]
  
  UK_idx <- which( theCountries$name_en[countries_to_label] == "United Kingdom")
  temp_matrix[UK_idx,] <- temp_matrix[UK_idx,] - 0.03*temp_matrix[UK_idx,]
  
  SWE_idx <- which( theCountries$name_en[countries_to_label] == "Sweden")
  temp_matrix[SWE_idx,] <- temp_matrix[SWE_idx,] - 0.022*temp_matrix[SWE_idx,]
  
  FIN_idx <- which( theCountries$name_en[countries_to_label] == "Finland")
  temp_matrix[FIN_idx,1] <- temp_matrix[FIN_idx,1] + 0.01*temp_matrix[FIN_idx,1]
  temp_matrix[FIN_idx,2] <- temp_matrix[FIN_idx,2] - 0.02*temp_matrix[FIN_idx,2]
  
  AUST_idx <- which( theCountries$name_en[countries_to_label] == "Austria")
  temp_matrix[AUST_idx,1] <- temp_matrix[AUST_idx,1] + 0.0075*temp_matrix[AUST_idx,1]
  
  GRE_idx <- which( theCountries$name_en[countries_to_label] == "Greece")
  temp_matrix[GRE_idx,1] <- temp_matrix[GRE_idx,1] - 0.016*temp_matrix[GRE_idx,1]
  
  FRA_idx <- which( theCountries$name_en[countries_to_label] == "France")
  temp_matrix[FRA_idx,] <- temp_matrix[FRA_idx,] + 0.15*temp_matrix[FRA_idx,]
  temp_matrix[FRA_idx,2] <- temp_matrix[FRA_idx,2] + 0.075*temp_matrix[FRA_idx,2]
  
  NOR_idx <- which( theCountries$name_en[countries_to_label] == "Norway")
  temp_matrix[NOR_idx,1] <- temp_matrix[NOR_idx,1] - 0.07*temp_matrix[NOR_idx,1]
  temp_matrix[NOR_idx,2] <- temp_matrix[NOR_idx,2] - 0.20*temp_matrix[NOR_idx,2]
  
  LAT_idx <- which( theCountries$name_en[countries_to_label] == "Latvia")
  temp_matrix[LAT_idx,1] <- temp_matrix[LAT_idx,1] + 0.015*temp_matrix[LAT_idx,1]
  
  
  ### Plot 3 posterior samples of the Matern random field, labelled (A), (B), and (C)
  if ( ii %in% c( which(sample_idx == 250),
                  which(sample_idx == 1132), 
                  which(sample_idx == 1559) ) ) { 
    
    png( filename = here::here( paste0("Figures/mercurySoil_Posterior_Sample_",sample_idx[ii],".png") ),
         res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
    
    mapmisc::map.new(hgm_2, asp = 0.7)
    hgm_Col = mapmisc::colourScale( x = hgm_2, dec = 2, style = "fixed", rev = TRUE,
                                    col = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")),
                                    breaks = c( -6.0, -3.0, -1.5, -0.5, 0, 0.5, 1.5, 3.0, 6.0 ) )
    plot( raster::mask(hgm_2, theCountries), breaks = hgm_Col$breaks,
          col = hgm_Col$col, add = TRUE, legend = FALSE)
    
    mapmisc::legendBreaks( pos = 'topright', hgm_Col, cex = 1.60,
                           bty = 'n', y.intersp = 1.15, inset = c(0.0,0.175) )
    plot( theCountries, col = '#FFFFFF00', add = TRUE )
    
    text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 1.90,
          labels = sapply( theCountries$name_en[countries_to_label], 
                           function(xx) { as.expression( bquote (bold(.(xx)) ) ) } ) )
    
    if ( ii == which( sample_idx == 250 ) ) {
      
      text( x = temp_matrix[1,1]-1200*1000, y = temp_matrix[1,2]+250*1000,
            labels = as.expression( bquote (bold("(A)") ) ), cex = 1.90 )
      
    }
    
    if ( ii == which( sample_idx == 1132 ) ) {
      
      text( x = temp_matrix[1,1]-1200*1000, y = temp_matrix[1,2]+250*1000,
            labels = as.expression( bquote (bold("(B)") ) ), cex = 1.90 )
      
    }
    
    if ( ii == which( sample_idx == 1559 ) ) {
      
      text( x = temp_matrix[1,1]-1200*1000, y = temp_matrix[1,2]+250*1000,
            labels = as.expression( bquote (bold("(C)") ) ), cex = 1.90 )
      
    }
    
    dev.off()
    
  }
  
  Us_conditional[,ii] <- hgm_2@data@values
  
} # Close for loop



##########################################
### Mean of Posterior Samples: Random Field
##########################################
ii = 500
mercurySoil_stan = c( range = 1000 * phi$phi_2[ii],
                      shape = rep(1.5,2*1000)[ii],
                      variance = phi$phi_1[ii]^2,
                      anisoAngleDegrees = -180/pi * phi$phi_4[ii] / 2,
                      anisoRatio = phi$phi_3[ii]+1 )

hgm_Us <- hgm
hgm_Us@bbox[1,1] <- hgm@bbox[1,1]-0.1*hgm@bbox[1,1]
hgm_Us@data <- as.data.frame( calculate_Us(ii) )

hgm_2 = RandomFields::RFsimulate(
  model = RandomFields::RMmatern( 
    var = mercurySoil_stan['variance'][[1]], 
    scale = mercurySoil_stan['range'][[1]], 
    nu = mercurySoil_stan['shape'][[1]], 
    Aniso = matrix( ncol = 2, c( sqrt( mercurySoil_stan['anisoRatio'][[1]] ), 0,
                                 0, 1 / sqrt( mercurySoil_stan['anisoRatio'][[1]] ) ) ) %*%
      matrix( ncol = 2, c( cos(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]]), sin(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]] ), 
                           -sin(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]]), cos(pi/180*mercurySoil_stan['anisoAngleDegrees'][[1]] ) ) ) ),
  x = raster::rasterToPoints( geostatsp::squareRaster(hgm_Us, raster_cells) )[,1],
  y = raster::rasterToPoints( geostatsp::squareRaster(hgm_Us, raster_cells) )[,2],
  data = hgm_Us, seed = fixed_seed )

hgm_2 <- raster::raster( hgm_2, vals = hgm_2$V1,
                         nrows = geostatsp::squareRaster(hgm_Us, raster_cells)@nrows, 
                         ncol = geostatsp::squareRaster(hgm_Us, raster_cells)@ncols ) 

### Adjust location of country labels to fit space available on the screen (for the map)
countries_to_label <- c(1,2,3,4,5,8,9,10,11,19,21,22,23,25)
temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
temp_matrix <- temp_matrix[c(1,2,3,4,5,8,9,10,11,19,21,22,23,25),]

UK_idx <- which( theCountries$name_en[countries_to_label] == "United Kingdom")
temp_matrix[UK_idx,] <- temp_matrix[UK_idx,] - 0.03*temp_matrix[UK_idx,]

SWE_idx <- which( theCountries$name_en[countries_to_label] == "Sweden")
temp_matrix[SWE_idx,] <- temp_matrix[SWE_idx,] - 0.022*temp_matrix[SWE_idx,]

FIN_idx <- which( theCountries$name_en[countries_to_label] == "Finland")
temp_matrix[FIN_idx,1] <- temp_matrix[FIN_idx,1] + 0.01*temp_matrix[FIN_idx,1]
temp_matrix[FIN_idx,2] <- temp_matrix[FIN_idx,2] - 0.02*temp_matrix[FIN_idx,2]

AUST_idx <- which( theCountries$name_en[countries_to_label] == "Austria")
temp_matrix[AUST_idx,1] <- temp_matrix[AUST_idx,1] + 0.0075*temp_matrix[AUST_idx,1]

GRE_idx <- which( theCountries$name_en[countries_to_label] == "Greece")
temp_matrix[GRE_idx,1] <- temp_matrix[GRE_idx,1] - 0.016*temp_matrix[GRE_idx,1]

FRA_idx <- which( theCountries$name_en[countries_to_label] == "France")
temp_matrix[FRA_idx,] <- temp_matrix[FRA_idx,] + 0.15*temp_matrix[FRA_idx,]
temp_matrix[FRA_idx,2] <- temp_matrix[FRA_idx,2] + 0.075*temp_matrix[FRA_idx,2]

NOR_idx <- which( theCountries$name_en[countries_to_label] == "Norway")
temp_matrix[NOR_idx,1] <- temp_matrix[NOR_idx,1] - 0.07*temp_matrix[NOR_idx,1]
temp_matrix[NOR_idx,2] <- temp_matrix[NOR_idx,2] - 0.20*temp_matrix[NOR_idx,2]

LAT_idx <- which( theCountries$name_en[countries_to_label] == "Latvia")
temp_matrix[LAT_idx,1] <- temp_matrix[LAT_idx,1] + 0.015*temp_matrix[LAT_idx,1]

### Calculate and use the mean of the posterior samples of the Matern random field
hgm_3 <- hgm_2 
hgm_3@data@values <- apply( Us_conditional, 1, mean )

### Plot the mean of the Matern random field 
png( filename = here::here( "Figures/mercurySoil_Posterior_Mean.png" ),
     res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )

  mapmisc::map.new(hgm_3, asp = 0.7)
  hgm_Col = mapmisc::colourScale( x = hgm_3, dec = 2, style = "fixed", rev = TRUE,
                                  col = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")),
                                  breaks = c( -6.0, -3.0, -1.5, -0.5, 0, 0.5, 1.5, 3.0, 6.0 ) )
  
  plot( raster::mask(hgm_3, theCountries), breaks = hgm_Col$breaks,
        col = hgm_Col$col, add = TRUE, legend = FALSE)
  
  mapmisc::legendBreaks( pos = 'topright', hgm_Col, cex = 1.60,
                         bty = 'n', y.intersp = 1.15, inset = c(0.0,0.175) )
  plot( theCountries, col = '#FFFFFF00', add = TRUE )
  
  text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 1.90,
        labels = sapply( theCountries$name_en[countries_to_label], 
                         function(xx) { as.expression( bquote (bold(.(xx)) ) ) } ) )
  text( x = temp_matrix[1,1]-1200*1000, y = temp_matrix[1,2]+250*1000,
        labels = as.expression( bquote (bold("(D)") ) ), cex = 1.90 )

dev.off()