
##########################
###
### Setup
###
##########################

library(geostatsp)
library(mcmcplots)
library(rstan)
library(dplyr)
library(actuar)
library(ggplot2)
library(sciscales)
library(ggnewscale)

source( "Helper_Functions.R" )
options( scipen = 999 )
data( swissRain ) 



##########################
##
## Load Stan models 
## Prepare data for plotting
##
##########################

#######################################
### Simulated Isotropic Scenario
#######################################
swissRain_stanfit_iso <- rstan::extract( readRDS(
  "Results/Isotropic_Simulation_Fit.RDS" ) )
swissRain_Stan_data_iso_plot <- swissRain_plotting_tibble( swissRain_stanfit_iso )


#######################################
### Simulated Anisotropic Scenario
#######################################
swissRain_stanfit_aniso <- rstan::extract( readRDS( 
  "Results/Anisotropic_Simulation_Fit.RDS" ) )
swissRain_Stan_data_aniso_plot <- swissRain_plotting_tibble( swissRain_stanfit_aniso )


#######################################
### Switzerland Rainfall Scenario 
#######################################
swissRain_stanfit_real <- rstan::extract( readRDS( 
  "Results/Switzerland_Rainfall_Fit.RDS" ) )
swissRain_Stan_data_plot <- swissRain_plotting_tibble( swissRain_stanfit_real )





##########################
##
## Traceplots
##
##########################


##########################################
### Isotropic Simulation  
##########################################
swissRain_samples_iso_list <- named_mcmc_list( readRDS(
  "Results/Isotropic_Simulation_Fit.RDS" ) )
save_traceplots( named_list = swissRain_samples_iso_list, 
                 file_name = "Figures/swissRain_Iso_Traceplot.png" ) 

##########################################
### Anisotropic Simulation 
##########################################
swissRain_samples_aniso_list <- named_mcmc_list( readRDS( 
  "Results/Anisotropic_Simulation_Fit.RDS" ) )
save_traceplots( named_list = swissRain_samples_aniso_list, 
                 file_name = "Figures/swissRain_Aniso_Traceplot.png" )

##########################################
### Switzerland Rainfall Data
##########################################
swissRain_samples_real_list <- named_mcmc_list( readRDS(
  "Results/Switzerland_Rainfall_Fit.RDS" ) )
save_traceplots( named_list = swissRain_samples_real_list, 
                 file_name = "Figures/swissRain_Traceplot.png" )



##########################
##
## Simulation Study
## Priors and Posteriors
## For Matern Parameters (Phi)
##
##########################


#####################
### Rho
#####################

prior_x <- seq(1, 10, length = 10 * 1000)
prior_y <- dexp(prior_x-1, rate = 1/4)
post_iso <- ks::kde.boundary(
  as.matrix(swissRain_Stan_data_iso_plot$Radius+1), 
  xmin = 1, xmax = 10, boundary.kernel = "beta",
  binned = FALSE, gridsize = 10 * 1000,
  h = 0.33 )
post_aniso <- density(swissRain_Stan_data_aniso_plot$Radius+1, 
                      bw = 0.25, n = 10 * 1000, 
                      from = 1, to = 10)
iso_true_phiR <- rep(1, 10 * 1000)
aniso_true_phiR <- rep(4.5, 10 * 1000)

Iso_Aniso_phiR <- tibble( "prior_x" = prior_x,
                          "prior_y" = prior_y, 
                          "post_iso_x" = post_iso$eval.points, 
                          "post_iso_y" = post_iso$estimate, 
                          "post_aniso_x" = post_aniso$x,
                          "post_aniso_y" = post_aniso$y, 
                          "true_iso_phiR" = iso_true_phiR, 
                          "true_aniso_phiR" = aniso_true_phiR )


png( filename = "Figures/swissRain_Iso_Aniso_PhiR.png", 
     res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )

  ggplot() + 
    geom_line( data = Iso_Aniso_phiR, linetype = "dotted", size = 1,
               aes( x = prior_x, y = prior_y, color = "black" ) ) +
    geom_line( data = Iso_Aniso_phiR, linetype = "solid", size = 1, alpha = 1,
               aes( x = post_iso_x, y = post_iso_y, color = "blue" ) ) +
    geom_line( data = Iso_Aniso_phiR, linetype = "solid", size = 1, alpha = 1,
               aes( x = true_iso_phiR, color = "red", 
                    y = seq(from = 0, to = 2.25, length = 10 * 1000) ) ) + 
    geom_line( data = Iso_Aniso_phiR, linetype = "solid", size = 1,
               aes( x = post_aniso_x, y = post_aniso_y, color = "green3" ) ) +
    geom_line( data = Iso_Aniso_phiR, linetype = "solid", size = 1, alpha = 0.75,
               aes( x = true_aniso_phiR, color = "purple", 
                    y = seq(from = 0, to = 2.25, length = 10 * 1000) ) ) + 
    scale_color_manual( values = c("black", "blue", "purple", "green3", "red" ) ) +
    xlab(expression(rho)) + ylab("Density") + xlim( c(1,10) ) +
    scale_x_continuous( breaks = c(1,4,7,10) ) +
    theme_classic() + 
    theme( panel.grid = element_blank(),
           axis.text = element_text(size = 20, color = "black"),
           axis.title = element_text(size = 22),
           axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
           axis.title.y = element_text(margin = margin(0, 12.5, 0, 0)),
           legend.title.align = 0.5,
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           legend.position = "none",
           legend.spacing.y = unit(0.50, 'cm') ) 

dev.off()


#####################
### Psi 
#####################
prior_x <- seq(-pi/2 * (180/pi), pi/2 * (180/pi), length = 512)
prior_y <- dunif(prior_x, min = -pi/2 * (180/pi), max = pi/2 * (180/pi))
post_iso <- density( swissRain_Stan_data_iso_plot$Angle, from = -90, to = 90 )
post_aniso <- density( swissRain_Stan_data_aniso_plot$Angle, bw = 1.00 )
iso_true_phiA <- rep(0,512)
aniso_true_phiA <- rep(-34,512)

Iso_Aniso_phiA <- tibble( "prior_x" = prior_x,
                          "prior_y" = prior_y, 
                          "post_iso_x" = post_iso$x, 
                          "post_iso_y" = post_iso$y, 
                          "post_aniso_x" = post_aniso$x,
                          "post_aniso_y" = post_aniso$y, 
                          "true_iso_phiA" = iso_true_phiA, 
                          "true_aniso_phiA" = aniso_true_phiA )



png( filename = "Figures/swissRain_Iso_Aniso_PhiA.png", 
     res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )

  ggplot() + 
    geom_line( data = Iso_Aniso_phiA, linetype = "dotted", size = 1,
               aes( x = prior_x, y = prior_y, color = "black" ) ) +
    geom_line( data = Iso_Aniso_phiA, linetype = "solid", size = 1,
               aes( x = post_iso_x, y = post_iso_y, color = "blue" ) ) +
    geom_line( data = Iso_Aniso_phiA, linetype = "solid", size = 1,
               aes( x = post_aniso_x, y = post_aniso_y, color = "green3" ) ) +
    geom_line( data = Iso_Aniso_phiA, linetype = "solid", size = 1, alpha = 0.75,
               aes( x = aniso_true_phiA, color = "purple", 
                    y = seq(from = 0,to = 0.15,length = 512) ) ) + 
    scale_color_manual( values = c("black", "blue", "purple", "green3", "red" ) ) +
    xlab(expression(psi~~"(Degrees)")) + ylab("Density") +
    theme_classic() + 
    theme( panel.grid = element_blank(),
           axis.text = element_text(size = 20, color = "black"),
           axis.title = element_text(size = 22),
           axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
           axis.title.y = element_text(margin = margin(0, 12.5, 0, 0)),
           legend.title.align = 0.5,
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           legend.position = "none",
           legend.spacing.y = unit(0.50, 'cm') ) 

dev.off()



#####################
### Sigma
#####################
prior_x <- seq(0, 4.5, length = 10 * 1000)
prior_y <- dexp(prior_x, rate = 1/2) 
post_iso <- density( swissRain_stanfit_iso$phi_1, 
                     bw = 0.15, n = 10 * 1000)
post_aniso <- density( swissRain_stanfit_aniso$phi_1, 
                       bw = 0.15, n = 10 * 1000)
iso_true_sigma <- rep(1, 10 * 1000)
aniso_true_sigma <- rep(1, 10 * 1000)

Iso_Aniso_sigma <- tibble( "prior_x" = prior_x,
                           "prior_y" = prior_y, 
                           "post_iso_x" = post_iso$x, 
                           "post_iso_y" = post_iso$y, 
                           "post_aniso_x" = post_aniso$x,
                           "post_aniso_y" = post_aniso$y, 
                           "true_iso_sigma" = iso_true_sigma, 
                           "true_aniso_sigma" = aniso_true_sigma )


png( filename = "Figures/swissRain_Iso_Aniso_sigma.png", 
     res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )

  ggplot() + 
    geom_line( data = Iso_Aniso_sigma, linetype = "dotted", size = 1,
               aes( x = prior_x, y = prior_y, color = "black" ) ) +
    scale_color_manual( values = "black", label = "Prior Density",
                        guide = guide_legend(order = 1, title = "",
                                             override.aes = list(size = 1.2) ) ) +
    new_scale_color() +
    geom_line( data = Iso_Aniso_sigma, linetype = "solid", alpha = 1,
               aes( x = post_iso_x, y = post_iso_y, color = "blue",  size = 1 ) ) +
    geom_line( data = Iso_Aniso_sigma, linetype = "solid", alpha = 1, 
               aes( x = true_iso_sigma, color = "red", size = 4,
                    y = seq(from = 0, to = 2.25, length = 10 * 1000) ) ) + 
    scale_color_manual( values = c("blue","red"), 
                        label = c("Posterior Density", "True Value"), 
                        guide = guide_legend(order = 2, keyheight = 1.75,
                                             title = "Isotropic Simulation",
                                             override.aes = list(size = 1.2) ) ) +
    scale_size( range = c(1.2,5), guide = "none" ) +
    new_scale_color() +
    geom_line( data = Iso_Aniso_sigma, linetype = "solid", size = 1,
               aes( x = post_aniso_x, y = post_aniso_y, color = "green3" ) ) +
    geom_line( data = Iso_Aniso_sigma, linetype = "solid", size = 2.5, alpha = 1,
               aes( x = true_aniso_sigma, color = "purple", 
                    y = seq(from = 0, to = 2.25, length = 10 * 1000) ) ) + 
    scale_color_manual( values = c("purple", "green3"), 
                        label = c("Posterior Density", "True Value"),
                        guide = guide_legend(order = 3, keyheight = 1.75,
                                             title = "Anisotropic Simulation",
                                             override.aes = list(size = 1.2) ) ) + 
    xlab( expression(sigma) ) + ylab("Density") + 
    coord_cartesian( xlim = c(0,4), ylim = c(0,2.25) ) +
    theme_classic() + 
    theme( panel.grid = element_blank(),
           axis.text = element_text(size = 20, color = "black"),
           axis.title = element_text(size = 22),
           axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
           axis.title.y = element_text(margin = margin(0, 12.5, 0, 0)),
           legend.title.align = 0.5,
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           legend.position = c(0.70, 0.75), # "none",
           legend.spacing.y = unit(0.50, 'cm') ) 

dev.off()



#####################
### Gamma
#####################
prior_x <- seq(from = 0, to = 225, length = 512)
prior_y <- dinvexp(prior_x, scale = 50) 
weighted.mean(prior_x, prior_y)
post_iso <- density( swissRain_Stan_data_iso_plot$Gamma, bw = 5 )
post_aniso <- density( swissRain_Stan_data_aniso_plot$Gamma, bw = 5 )
iso_true_gamma <- rep(60,512)
aniso_true_gamma <- rep(60,512)

Iso_Aniso_gamma <- tibble( "prior_x" = prior_x,
                           "prior_y" = prior_y, 
                           "post_iso_x" = post_iso$x, 
                           "post_iso_y" = post_iso$y, 
                           "post_aniso_x" = post_aniso$x,
                           "post_aniso_y" = post_aniso$y, 
                           "true_iso_gamma" = iso_true_gamma, 
                           "true_aniso_gamma" = aniso_true_gamma )

png( filename = "Figures/swissRain_Iso_Aniso_gamma.png", 
     res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )

  ggplot() + 
    geom_line( data = Iso_Aniso_gamma, linetype = "dotted", 
               aes( x = prior_x, y = prior_y, color = "black", size = 1 ) ) +
    geom_line( data = Iso_Aniso_gamma, linetype = "solid", 
               aes( x = post_iso_x, y = post_iso_y, color = "blue", size = 1 ) ) +
    geom_line( data = Iso_Aniso_gamma, linetype = "solid", alpha = 1,
               aes( x = iso_true_gamma, color = "red", size = 4,
                    y = seq(from = 0, to = 0.04, length = 512) ) ) + 
    geom_line( data = Iso_Aniso_gamma, linetype = "solid", 
               aes( x = post_aniso_x, y = post_aniso_y, size = 1, color = "green3" ) ) +
    geom_line( data = Iso_Aniso_gamma, linetype = "solid", alpha = 1,
               aes( x = aniso_true_gamma, color = "purple", size = 1.33,  
                    y = seq(from = 0,to = 0.04, length = 512) ) ) + 
    scale_color_manual( values = c("black", "blue", "purple", "green3", "red" ),
                        guide = guide_legend( override.aes = list(size = 1.2) ) ) +
    scale_size( range = c(1.2,5), guide = "none" ) +
    xlab(expression(gamma)) + ylab("Density") + 
    coord_cartesian(xlim = c(0,200), ylim = c(0,0.04) ) +
    theme_classic() + 
    theme( panel.grid = element_blank(),
           axis.text = element_text(size = 20, color = "black"),
           axis.title = element_text(size = 22),
           axis.title.x = element_text(margin = margin(15, 0, 0, 0)),
           axis.title.y = element_text(margin = margin(0, 12.5, 0, 0)),
           legend.title.align = 0.5,
           legend.title = element_text(size = 20),
           legend.text = element_text(size = 18),
           legend.position = "none", # c(0.60,0.80), 
           legend.spacing.y = unit(0.50, 'cm') ) 

dev.off()




##########################
##
## Summary Tables
## Model Parameters
##
##########################


#######################################
### Isotropic Simulation
#######################################
summary_table_natural_params( stan_samples = readRDS( 
  "Results/Isotropic_Simulation_Fit.RDS" ), 
  extracted_samples = swissRain_stanfit_iso ) %>% 
  kableExtra::kable_styling()
summary_table_polar_params( stan_samples = readRDS( 
  "Results/Isotropic_Simulation_Fit.RDS" ) )  %>% 
  kableExtra::kable_styling()


#######################################
### Anisotropic Simulation
#######################################
summary_table_natural_params( stan_samples = readRDS( 
  "Results/Anisotropic_Simulation_Fit.RDS" ), 
  extracted_samples = swissRain_stanfit_aniso ) %>% 
  kableExtra::kable_styling()
summary_table_polar_params( stan_samples = readRDS( 
  "Results/Anisotropic_Simulation_Fit.RDS" ) ) %>% 
  kableExtra::kable_styling()


#######################################
### Switzerland Rainfall Data
#######################################
summary_table_natural_params( stan_samples = readRDS(
  "Results/Switzerland_Rainfall_Fit.RDS" ), 
  extracted_samples = swissRain_stanfit_real ) %>% 
  kableExtra::kable_styling()
summary_table_polar_params( stan_samples = readRDS(
  "Results/Switzerland_Rainfall_Fit.RDS" ) ) %>% 
  kableExtra::kable_styling()





##########################
##
## Range Parameters
## (Including Axis-Aligned)
##
##########################

###########################
### Range
###########################

### Isotropic Simulation
range_iso <- c( round( mean( swissRain_stanfit_iso$phi_2 ), 2 ), 
                round( quantile( swissRain_stanfit_iso$phi_2, probs = 0.025 ), 2 ), 
                round( quantile( swissRain_stanfit_iso$phi_2, probs = 0.50 ), 2 ), 
                round( quantile( swissRain_stanfit_iso$phi_2, probs = 0.975 ), 2 ) )

### Anisotropic Simulation 
range_aniso <- c( round( mean( swissRain_stanfit_aniso$phi_2 ), 2 ),
                  round( quantile( swissRain_stanfit_aniso$phi_2, probs = 0.025 ), 2 ), 
                  round( quantile( swissRain_stanfit_aniso$phi_2, probs = 0.50 ), 2 ), 
                  round( quantile( swissRain_stanfit_aniso$phi_2, probs = 0.975 ), 2 ) )

### Switzerland Rainfall data 
range_real <- c( round( mean( swissRain_stanfit_real$phi_2 ), 2 ), 
                 round( quantile( swissRain_stanfit_real$phi_2, probs = 0.025 ), 2 ),
                 round( quantile( swissRain_stanfit_real$phi_2, probs = 0.50 ), 2 ), 
                 round( quantile( swissRain_stanfit_real$phi_2, probs = 0.975 ), 2 ) ) 


###########################
### Direction: psi 
### Scaled by >= 1
###########################

### Isotropic Simulation
range_y_iso <- c( round( mean( swissRain_stanfit_iso$phi_2 * sqrt(swissRain_stanfit_iso$phi_3+1) ), 2 ),
                  round( quantile( swissRain_stanfit_iso$phi_2 * sqrt(swissRain_stanfit_iso$phi_3+1), probs = 0.025 ), 2 ), 
                  round( quantile( swissRain_stanfit_iso$phi_2 * sqrt(swissRain_stanfit_iso$phi_3+1), probs = 0.50 ), 2 ), 
                  round( quantile( swissRain_stanfit_iso$phi_2 * sqrt(swissRain_stanfit_iso$phi_3+1), probs = 0.975 ), 2 ) )

### Anisotropic Simulation
range_y_aniso <- c( round( mean( swissRain_stanfit_aniso$phi_2 * sqrt(swissRain_stanfit_aniso$phi_3+1) ), 2 ), 
                    round( quantile( swissRain_stanfit_aniso$phi_2 * sqrt(swissRain_stanfit_aniso$phi_3+1), probs = 0.025 ), 2 ),
                    round( quantile( swissRain_stanfit_aniso$phi_2 * sqrt(swissRain_stanfit_aniso$phi_3+1), probs = 0.50 ), 2 ), 
                    round( quantile( swissRain_stanfit_aniso$phi_2 * sqrt(swissRain_stanfit_aniso$phi_3+1), probs = 0.975 ), 2 ) )

### Switzerland Rainfall data 
range_y_real <- c( round( mean( swissRain_stanfit_real$phi_2 * sqrt(swissRain_stanfit_real$phi_3+1) ), 2 ),
                   round( quantile( swissRain_stanfit_real$phi_2 * sqrt(swissRain_stanfit_real$phi_3+1), probs = 0.025 ), 2 ),
                   round( quantile( swissRain_stanfit_real$phi_2 * sqrt(swissRain_stanfit_real$phi_3+1), probs = 0.50 ), 2 ), 
                   round( quantile( swissRain_stanfit_real$phi_2 * sqrt(swissRain_stanfit_real$phi_3+1), probs = 0.975 ), 2 ) )


###########################
### Direction: psi + 90 deg
### Scaled by <= 1
###########################

### Isotropic Simulation 
range_x_iso <- c( round( mean( swissRain_stanfit_iso$phi_2 / sqrt(swissRain_stanfit_iso$phi_3+1) ), 2 ), 
                  round( quantile( swissRain_stanfit_iso$phi_2 / sqrt(swissRain_stanfit_iso$phi_3+1), probs = 0.025 ), 2 ),
                  round( quantile( swissRain_stanfit_iso$phi_2 / sqrt(swissRain_stanfit_iso$phi_3+1), probs = 0.50 ), 2 ), 
                  round( quantile( swissRain_stanfit_iso$phi_2 / sqrt(swissRain_stanfit_iso$phi_3+1), probs = 0.975 ), 2 ) ) 

### Anisotropic Simulation 
range_x_aniso <- c( round( mean( swissRain_stanfit_aniso$phi_2 / sqrt(swissRain_stanfit_aniso$phi_3+1) ), 2 ),
                    round( quantile( swissRain_stanfit_aniso$phi_2 / sqrt(swissRain_stanfit_aniso$phi_3+1), probs = 0.025 ), 2 ), 
                    round( quantile( swissRain_stanfit_aniso$phi_2 / sqrt(swissRain_stanfit_aniso$phi_3+1), probs = 0.50 ), 2 ), 
                    round( quantile( swissRain_stanfit_aniso$phi_2 / sqrt(swissRain_stanfit_aniso$phi_3+1), probs = 0.975 ), 2 ) ) 

### Switzerland Rainfall  
range_x_real <- c( round( mean( swissRain_stanfit_real$phi_2 / sqrt(swissRain_stanfit_real$phi_3+1) ), 2 ), 
                   round( quantile( swissRain_stanfit_real$phi_2 / sqrt(swissRain_stanfit_real$phi_3+1), probs = 0.025 ), 2 ),
                   round( quantile( swissRain_stanfit_real$phi_2 / sqrt(swissRain_stanfit_real$phi_3+1), probs = 0.50 ), 2 ), 
                   round( quantile( swissRain_stanfit_real$phi_2 / sqrt(swissRain_stanfit_real$phi_3+1), probs = 0.975 ), 2 ) ) 


###########################
### Summary Tables
###########################

### Isotropic Simulation
tibble( "Range" = range_iso,
        "Range (Y)" = range_y_iso,
        "Range (X)" = range_x_iso ) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling()

### Anisotropic Simulation
tibble( "Range" = range_aniso,
        "Range (Y)" = range_y_aniso,
        "Range (X)" = range_x_aniso ) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling()

### Switzerland Rainfall 
tibble( "Range" = range_real,
        "Range (Y)" = range_y_real,
        "Range (X)" = range_x_real ) %>% 
  knitr::kable() %>% 
  kableExtra::kable_styling()



##########################
##
## Rhat Values 
##
##########################

swissRain_stanfit_iso_2 <- readRDS( "Results/Isotropic_Simulation_Fit.RDS" )
round( summary( swissRain_stanfit_iso_2,
                pars = c( "sqrt_inv_shape", "Beta[1]", "Beta[2]",
                          "phi_1", "phi_2", "phi_3", "phi_4", 
                          "theta_1", "theta_2", "theta_3", "theta_4") )$summary[,10], 2 )
rm(swissRain_stanfit_iso_2); gc()


swissRain_stanfit_aniso_2 <- readRDS( "Results/Anisotropic_Simulation_Fit.RDS" )
round( summary( swissRain_stanfit_aniso_2,
                pars = c( "sqrt_inv_shape", "Beta[1]", "Beta[2]",
                          "phi_1", "phi_2", "phi_3", "phi_4", 
                          "theta_1", "theta_2", "theta_3", "theta_4") )$summary[,10], 2 )
rm(swissRain_stanfit_aniso_2); gc()


swissRain_stanfit_2 <- readRDS( "Results/Switzerland_Rainfall_Fit.RDS" )
round( summary( swissRain_stanfit_2,
                pars = c( "sqrt_inv_shape", "Beta[1]", "Beta[2]",
                          "phi_1", "phi_2", "phi_3", "phi_4", 
                          "theta_1", "theta_2", "theta_3", "theta_4") )$summary[,10], 2 )
rm(swissRain_stanfit_2); gc()