###
### Helper-Functions
###

#####################
### Tibble Functions
#####################
swissRain_plotting_tibble <- function(stanfit_object) { 
  
  return( tibble( "Radius" = stanfit_object$phi_3, # was "Ratio"
                  "Angle" = -180/pi * stanfit_object$phi_4 / 2,
                  "Theta_X" = stanfit_object$theta_3,
                  "Theta_Y" = stanfit_object$theta_4,
                  "Sigma" = stanfit_object$phi_1, 
                  "Gamma" = stanfit_object$phi_2, 
                  "Theta_1" = stanfit_object$theta_1, 
                  "Theta_2" = stanfit_object$theta_2 ) ) 
}


#######################
### Traceplot Functions
#######################
rename_traceplots <- function(a_list) {
  
  colnames(a_list) <- c( 
    paste( expression(beta[0]) ), 
    paste( expression(beta[1]) ), 
    paste( expression(1/sqrt(alpha)) ), 
    paste( expression(sigma) ),
    paste( expression(gamma) ), 
    # paste( expression(gamma %.% sqrt(rho)) ), 
    # paste( expression(gamma / sqrt(rho)) ), 
    paste( expression(rho) ), 
    paste( expression(psi) ),
    paste( expression(theta[1]) ),
    paste( expression(theta[2]) ),
    paste( expression(theta[3]) ),
    paste( expression(theta[4]) ) )
  
  return(a_list)
}

named_mcmc_list <- function(stan_samples) {
  samples_list <- rstan::As.mcmc.list( stan_samples, 
                                       pars = c('Beta[1]', 'Beta[2]', 'sqrt_inv_shape', 
                                                'phi_1', 'phi_2', 'phi_3', 'phi_4', # 'phi_X', 'phi_Y', 
                                                'theta_4', 'theta_3', 'theta_1', 'theta_2') )
  
  for (ii in 1:4) { 
    
    ### Report anisotropic angle, not polar angle 
    samples_list[[ii]][,7] <- -180/ pi * c( samples_list[[ii]][,7] ) / 2
    ### Report theta_1 in the correct direction of rotation
    samples_list[[ii]][,8] <- -c( samples_list[[ii]][,8] ) 
    
    ### Report rho (the anisotropic ratio), not rho-1
    samples_list[[ii]][,6] <- c( samples_list[[ii]][,6]+1 ) 
    
    ### Report correct scaling for range in (rotated) X-axis and Y-axis directions
    # samples_list[[ii]][,5] <- samples_list[[ii]][,5] / sqrt( samples_list[[ii]][,7] + 1 )  
    # samples_list[[ii]][,6] <- samples_list[[ii]][,6] / sqrt( samples_list[[ii]][,7] + 1 ) 
  }
  
  samples_list <- lapply( samples_list, rename_traceplots )
  return( samples_list )
  
}

save_traceplots <- function( named_list, file_name ) {
  
  png( file = file_name, res = 300, units = "px",
       width = 1.5 * 1928, height = 1.5 * 1514 )
  traplot( named_list, greek = TRUE, style = "plain",
           col = c( "#E66101", "#998EC3", "#542788", "#F1A340" ),
           cex.axis = 1.10, cex.main = 1.80 )
  dev.off()
  
}



#######################
### Summary Table Functions
#######################
summary_table_natural_params <- function( stan_samples, extracted_samples ) {
  summary_stats_natural <- summary( stan_samples, 
                                    pars = c('Beta[1]', 'Beta[2]', 'phi_1', 'phi_2', 
                                             'phi_3', 'phi_4', 'sqrt_inv_shape') )$summary[,c(1,4,6,8)]
  
  colnames( summary_stats_natural ) <- c( "Mean", "2.5\\\\% Q", "Median", "97.5\\\\% Q" ) 
  rownames( summary_stats_natural ) <- c( "Intercept: $\\beta_{0}$",
                                          "Altitude: $\\beta_{1}$",
                                          "Spatial Std. Dev.: $\\phi_{1}$", 
                                          "Range: $\\phi_{2}$", 
                                          "Ratio: $\\phi_{3}$+1", 
                                          "Angle: $\\phi_{4}$",
                                          "Coefficient of Variation: $1/\\sqrt{\\alpha}$") 
  
  # Report anisotropic angle
  summary_stats_natural[6,1:4] <- -180/pi * summary_stats_natural[6,1:4] / 2 
  summary_stats_natural[6,c(2,4)] <- summary_stats_natural[6,c(4,2)]
  
  # Report HDI for shifted anisotropic ratio
  summary_stats_natural[5,c(1,3)] <- summary_stats_natural[5,c(1,3)] + 1
  summary_stats_natural[5,c(2,4)] <- c( bayestestR::hdi( extracted_samples$phi_3 )$CI_low+1,
                                        bayestestR::hdi( extracted_samples$phi_3 )$CI_high+1 )
  
  return( knitr::kable( summary_stats_natural,
                        digits = 2, format = "html" ) )
}


summary_table_polar_params <- function( stan_samples ) {
  summary_stats_polar <- summary( stan_samples, 
                                  pars =  c('theta_4', 'theta_3', 'theta_1', 'theta_2') )$summary[,c(1,4,6,8)]
  
  colnames( summary_stats_polar ) <- c( "Mean", "2.5\\\\% Q", "Median", "97.5\\\\% Q" ) 
  rownames( summary_stats_polar ) <- c( "Cartesian Coordinate: $\\theta_{1}$", 
                                        "Cartesian Coordinate: $\\theta_{2}$",
                                        "Sum: $\\theta_{3}$",
                                        "Difference: $\\theta_{4}$" )
  
  # Report theta_1 in the correct direction of rotation 
  summary_stats_polar[1,1:4] <- -summary_stats_polar[1,1:4] 
  summary_stats_polar[1,c(2,4)] <- summary_stats_polar[1,c(4,2)] 
  
  return( knitr::kable( summary_stats_polar,
                        digits = 2, format = "html" ) )
}
