###########################
###
### Setup
###
###########################

### Load libraries 
library(rstan)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(MASS)
library(hdrcde)
library(LogConcDEAD)

source( "Helper_Functions.R" )

swissRain_samples_aniso <- readRDS( "Results/Anisotropic_Simulation_Fit.RDS" ) 
swissRain_stanfit_aniso <- rstan::extract( swissRain_samples_aniso )
swissRain_Stan_data_aniso_plot <- swissRain_plotting_tibble( swissRain_stanfit_aniso )

# TRUE: Load saved density estimates
# FALSE: Perform density estimation
load_saved_LCD <- FALSE

### Rounding to perform on input to LCD estimator
LCD_digits <- 2

### TRUE: Add legends to plots of posterior densities
### FALSE: Do not include legends in plots of posterior densities
include_legend <- TRUE



#######################################
###
### Anisotropic Angle and Ratio
###
#######################################

#######################################
### Density estimation of (theta_1,theta_2)
#######################################

### Perform density estimation (load_saved_LCD <- FALSE)
if ( !load_saved_LCD ) {
  theta_XY <- cbind( -swissRain_Stan_data_aniso_plot$Theta_Y,
                     swissRain_Stan_data_aniso_plot$Theta_X )
  theta_XY_rounded <- round( theta_XY, digits = LCD_digits )
  theta_XY_weights <- getweights( theta_XY_rounded )
  theta_XY_lcd <- mlelcd( x = theta_XY_weights$x,
                          w = theta_XY_weights$w )
  saveRDS( theta_XY_lcd, file = "Results/Aniso_Sim_Polar_LCD.RDS" )
}

### Load saved density estimates (load_saved_LCD == TRUE)
if ( load_saved_LCD ) {
  theta_XY_lcd <- readRDS( "Results/Aniso_Sim_Polar_LCD.RDS" )
}
 

### Obtain grid points 
N = 100
h_theta <- c( MASS::bandwidth.nrd( -swissRain_Stan_data_aniso_plot$Theta_Y ),
              MASS::bandwidth.nrd( swissRain_Stan_data_aniso_plot$Theta_X ) )
aniso_grid_points_theta <- MASS::kde2d( -swissRain_Stan_data_aniso_plot$Theta_Y,
                                        swissRain_Stan_data_aniso_plot$Theta_X,
                                        h = h_theta, n = N )


### Evaluate density at grid points
log_concave_density_theta <- expand.grid("x" = aniso_grid_points_theta$x,
                                         "y" = aniso_grid_points_theta$y)
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


### Plot bivariate posterior density of (theta_1,theta_2)
png( file = "Figures/Aniso_Sim_Polar_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

  par( mar = c(5.1, 5.1, 2.1, 2.1) )
  hdrcde:::hdrcde.filled.contour(
    x = aniso_grid_points_theta$x,
    y = aniso_grid_points_theta$y,
    z = matrix( theta_XY_theta_eval, N, N ),
    levels = c( credible_region_levels_theta, max(theta_XY_theta_eval) ),
    color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
    xlim = c(-3.0,0), ylim = c(-1.0,2.0), 
    plot.axes = { title( main = "", 
                         xlab = expression( theta[1] ), 
                         ylab = expression( theta[2] ),
                         cex.lab = 1.80 )
      axis(1, seq(from = -3.0, to = 0, by = 0.5), cex.axis = 1.40) 
      axis(2, seq(from = -1.0, to = 2.0, by = 0.25), cex.axis = 1.40) } )
  if ( include_legend ) {
    legend( x = -0.875, y = -0.125, bg = "white",
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
h_phi <- c( MASS::bandwidth.nrd( swissRain_Stan_data_aniso_plot$Angle ),
            MASS::bandwidth.nrd( swissRain_Stan_data_aniso_plot$Radius ) )
aniso_grid_points_phi <- MASS::kde2d( swissRain_Stan_data_aniso_plot$Angle,
                                      swissRain_Stan_data_aniso_plot$Radius,
                                      h = h_phi, n = N )


### Transform grid points from phi-space to theta-space
x2 = list()
x2$x = matrix( NA, ncol = N, nrow = N )
x2$y = matrix( NA, ncol = N, nrow = N )

for (ii in 1:N) {
  for (jj in 1:N) {
    
    # theta_1
    x2$x[ii,jj] <- sqrt( aniso_grid_points_phi$y[ii] ) *
      sin( 2 * pi/180 * aniso_grid_points_phi$x[jj] ) 
    # theta_2
    x2$y[ii,jj] <- sqrt( aniso_grid_points_phi$y[ii] ) *
      cos( 2 * pi/180 * aniso_grid_points_phi$x[jj] ) 
    
  }
}

### Evaluate density on grid points in theta-space
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
png( file = "Figures/Aniso_Sim_Angle_Ratio_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

  par( mar = c(5.1, 5.1, 2.1, 2.1) )
  hdrcde:::hdrcde.filled.contour(
    x = x3$x[,1]+1, y = 180/pi * x3$y[1,], z = x3$z,
    levels = c( credible_region_levels_phi, max(x3$z) ),
    color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE),
    xlim = c(1,11), ylim = c(-55,-15), 
    plot.axes = { title( main = "", 
                         xlab = expression( rho ), 
                         ylab = expression( psi~~"(Degrees)" ),
                         cex.lab = 1.80 )
      axis(1, seq(from = 1, to = 11, by = 3), cex.axis = 1.40)
      axis(2, seq(from = -55, to = -15, by = 10), cex.axis = 1.40) } )
  if ( include_legend ) { 
    legend( x = 8, y = -17.5, bg = "white",
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
###
### Range and Std. Dev. of Random Field
###
#######################################

#######################################
### Density estimation of (theta_3,theta_4)
#######################################

### Perform density estimation (load_saved_LCD <- FALSE)
if ( !load_saved_LCD ) {
  theta_XY <- cbind( swissRain_stanfit_aniso$theta_1,
                     swissRain_stanfit_aniso$theta_2 )
  theta_XY_rounded <- round( theta_XY, digits = LCD_digits )
  theta_XY_weights <- getweights( theta_XY_rounded )
  theta_XY_lcd <- mlelcd(x = theta_XY_weights$x, w = theta_XY_weights$w)
  saveRDS( theta_XY_lcd, file = "Results/Aniso_Sim_TransIso_LCD.RDS" )
}

### Load saved density estimates (load_saved_LCD <- TRUE)
if ( load_saved_LCD ) {
  theta_XY_lcd <- readRDS( "Results/Aniso_Sim_TransIso_LCD.RDS" )  
}


### Obtain grid points 
N = 100
h_theta <- c( MASS::bandwidth.nrd( swissRain_stanfit_aniso$theta_1 ),
              MASS::bandwidth.nrd( swissRain_stanfit_aniso$theta_2 ) )
aniso_grid_points_theta <- MASS::kde2d( swissRain_stanfit_aniso$theta_1,
                                        swissRain_stanfit_aniso$theta_2,
                                        h = h_theta, n = N )

### Evaluate density at grid points 
log_concave_density_theta <- expand.grid("x" = aniso_grid_points_theta$x, 
                                         "y" = aniso_grid_points_theta$y) 
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
png( file = "Figures/Aniso_Sim_TransIso_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

  par( mar = c(5.1, 5.1, 2.1, 2.1) )
  hdrcde:::hdrcde.filled.contour(
    x = aniso_grid_points_theta$x, 
    y = aniso_grid_points_theta$y, 
    z = matrix( theta_XY_theta_eval, N, N ),
    levels = c( credible_region_levels_theta, max(theta_XY_theta_eval) ),
    color.palette = function(n) hcl.colors(n, "YlOrRd", rev = TRUE), 
    xlim = c(1.5,5.5), ylim = c(-7.5,-4.0), 
    plot.axes = { title( main = "", 
                         xlab = expression( theta[3] ), 
                         ylab = expression( theta[4] ),
                         cex.lab = 1.80 )
      axis(1, seq(from = 1.5, to = 5.5, by = 1.0), cex.axis = 1.40)
      axis(2, seq(from = -7.5, to = -4.0, by = 1.0), cex.axis = 1.40) } )
  if ( include_legend ) { 
    legend( x = 4.0, y = -6.5, bg = "white",
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

### Obtain grid points 
N = 100
h_phi <- c( MASS::bandwidth.nrd( swissRain_stanfit_aniso$phi_1 ),
            MASS::bandwidth.nrd( swissRain_stanfit_aniso$phi_2 ) )
aniso_grid_points_phi <- MASS::kde2d( swissRain_stanfit_aniso$phi_1, # sigma 
                                      swissRain_stanfit_aniso$phi_2, # gamma
                                      h = h_phi, n = N )


### Transform grid points from phi-space to theta-space
x2 = list()
x2$x = matrix( NA, ncol = N, nrow = N )
x2$y = matrix( NA, ncol = N, nrow = N )


for (ii in 1:N) {
  for (jj in 1:N) {
    
    # theta_3
    x2$x[ii,jj] <- log( aniso_grid_points_phi$x[ii] ) / sqrt(2) + 
                   log( aniso_grid_points_phi$y[jj] ) / sqrt(2)
    # theta_4
    x2$y[ii,jj] <- log( aniso_grid_points_phi$x[ii] ) * sqrt(2) - 
                   log( aniso_grid_points_phi$y[jj] ) * sqrt(2)
    
  }
}

### Evaluate density at grid points in theta space
theta_XY_eval_on_theta_grid <- dlcd( x = cbind( as.vector(x2$x), 
                                                as.vector(x2$y) ),
                                         lcd = theta_XY_lcd )


### Transform the sampled theta points to phi space
x3 = list()
x3$x = matrix( NA, ncol = N, nrow = N )
x3$y = matrix( NA, ncol = N, nrow = N )
x3$z = matrix( NA, ncol = N, nrow = N )

theta_XY_eval_on_theta_grid_matrix <- matrix( theta_XY_eval_on_theta_grid, 
                                              nrow = N, ncol = N )

for (ii in 1:N) {
  for (jj in 1:N) {
    
    ### Transform from theta space to phi space 
    x3$x[ii,jj] <- exp( ( x2$x[ii,jj] + x2$y[ii,jj]/2 ) / sqrt(2) ) # sigma 
    x3$y[ii,jj] <- exp( ( x2$x[ii,jj] - x2$y[ii,jj]/2 ) / sqrt(2) ) # gamma
    x3$z[ii,jj] <- theta_XY_eval_on_theta_grid_matrix[ii,jj] * ( 2 / ( x3$x[ii,jj] * x3$x[ii,jj] ) )
    
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


### Plot bivariate posterior density of (sigma,gamma)
png( file = "Figures/Aniso_Sim_Sigma_Gamma_LC2D.png",
     res = 300, width = 1.25 * 1928, height = 1.25 * 1514, units = "px" )

  par( mar = c(5.1, 5.1, 2.1, 2.1) )
  hdrcde:::hdrcde.filled.contour(
    x = x3$x[,1], y = x3$y[1,], z = x3$z,
    levels = c( credible_region_levels_phi, max(x3$z) ),
    col = rev( hcl.colors(6, "YlOrRd") ),
    xlim = c(0,3.5), ylim = c(0,200), 
    plot.axes = { title( main = "", 
                         xlab = expression( sigma ), 
                         ylab = expression( gamma ),
                         cex.lab = 1.80 )
      axis(1, seq(from = 0, to = 3.0, by = 0.5), cex.axis = 1.40)
      axis(2, seq(from = 0, to = 200, by = 50), cex.axis = 1.40) } )
  if ( include_legend ) {
    legend( x = 0.40, y = 180, bg = "white",
            legend = rev( scales::percent(
              c(0.25, 0.50, 0.80, 0.90, 0.95, 0.99),
              accuracy = 1 ) ),
            fill = rev( hcl.colors(6, "YlOrRd") ),
            title = "Credible Region", box.col = "white",
            ncol = 2, cex = 1.10 )
  }
  par( mar = c(5.1, 4.1, 4.1, 2.1) )

dev.off()
