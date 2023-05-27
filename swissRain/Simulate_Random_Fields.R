  #######################################
  ###
  ### Setup
  ### 
  #######################################
  
  library(foreach)
  library(doParallel)
  library(dplyr)
  library(geostatsp)
  library(rnaturalearth)
  
  library(rgdal)
  library(raster)
  library(sp)
  library(sf)
  
  data("swissRain")
  source( "Helper_Functions.R" )
  
  RandomFields::RFoptions(install="no")
  sf::sf_use_s2(FALSE)
  
  fixed_seed <- 2500
  set.seed(fixed_seed)
  
  raster_cells <- 7.5 * 100
  num_samples <- 2.5 * 100
  
  
  #######################################
  ###
  ### Simulate Random Field
  ### Isotropic Model
  ### 
  #######################################
  
  swissRain_samples_iso <- readRDS( "Results/Isotropic_Simulation_Fit.RDS" )
  swissRain_stanfit_iso <- rstan::extract( swissRain_samples_iso )
  swissRain_Stan_data_iso_plot <- swissRain_plotting_tibble( swissRain_stanfit_iso )
  
  calculate_Us_iso <- function(kk) { 
    return( chol( swissRain_stanfit_iso$Sigma[kk,,] ) %*% 
              swissRain_stanfit_iso$Alpha[kk,] ) 
  }
  
  
  swissAltitude_vector <- raster::extract(swissAltitude, swissRain) / 1000 
  
  cl <- makeCluster( parallel::detectCores() / 2 )
  registerDoParallel(cl)
  
  swissMap = sp::spTransform( rnaturalearth::ne_countries( country = 'Switzerland', 
                                                           type = "countries", scale = "medium" ),
                              projection(swissRain) )
  
  num_Us_samples <- num_samples
  sample_idx <- sample(1:4000, num_samples) 
  
  Us_iso_conditional <- foreach ( ii = sample_idx, .combine = cbind ) %dopar% {
    
    iso_stan_2 = c( range = 1000 * swissRain_stanfit_iso$phi_2[ii],
                    shape = rep(2,4*1000)[ii],
                    variance = swissRain_stanfit_iso$phi_1[ii]^2,
                    anisoAngleDegrees = swissRain_Stan_data_iso_plot$Angle[ii],
                    anisoRatio = swissRain_Stan_data_iso_plot$Radius[ii]+1 )
    
    swissRain_Us <- swissRain 
    swissRain_Us@data <- as.data.frame( calculate_Us_iso(ii) )
    
    swissSim_2 = RandomFields::RFsimulate(
      model = RandomFields::RMmatern( 
        var = iso_stan_2['variance'][[1]], 
        scale = iso_stan_2['range'][[1]], 
        nu = iso_stan_2['shape'][[1]], 
        Aniso = matrix( ncol = 2, c( sqrt( iso_stan_2['anisoRatio'][[1]] ), 0,
                                     0, 1 / sqrt( iso_stan_2['anisoRatio'][[1]] ) ) ) %*%
          matrix( ncol = 2, c( cos(pi/180*iso_stan_2['anisoAngleDegrees'][[1]]), sin(pi/180*iso_stan_2['anisoAngleDegrees'][[1]] ), 
                               -sin(pi/180*iso_stan_2['anisoAngleDegrees'][[1]]), cos(pi/180*iso_stan_2['anisoAngleDegrees'][[1]] ) ) ) ),
      x = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,1],
      y = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,2],
      data = swissRain_Us, seed = fixed_seed )
    
    swissSim_2 <- raster::raster( swissSim_2, vals = swissSim_2$V1,
                                  nrows = geostatsp::squareRaster(swissMap, raster_cells)@nrows, 
                                  ncol = geostatsp::squareRaster(swissMap, raster_cells)@ncols ) 
    
    swissMap_sf = rnaturalearth::ne_countries( country = 'Switzerland', type = "countries", 
                                               scale = "medium", returnclass = 'sf' )
    
    temp <- sf::st_transform( sf::st_centroid(swissMap_sf)$geometry, crs = swissRain@proj4string )
    temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
    
    
    if (ii %in% c(754,2468,2915) ) {
    
      png( filename = paste0("Figures/swissRain_Iso_Posterior_Sample_",ii,".png"),
           res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
    
        mapmisc::map.new(swissBorder)
        swissCol = mapmisc::colourScale( x = swissSim_2, dec = 1, style = "fixed", rev = TRUE,
                                         col = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")),
                                         breaks = c( -8, -6, -4, -2, -1, 
                                                      0, 1, 2, 4, 6 ) )
        
        plot(raster::mask(swissSim_2, swissMap), breaks = swissCol$breaks,
             col = swissCol$col, add = TRUE, legend = FALSE)
        plot(swissMap, add = TRUE)
        mapmisc::legendBreaks( 'topleft', swissCol, cex = 1.60, inset = c(0.02,0.00),
                               bty = 'n', y.intersp = 1.15 )
        
        text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 2.50,
              labels = as.expression(bquote(bold(.(swissMap$sovereignt)))) ) 
        
        if (ii == 754) {
          
          text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
                labels = as.expression(bquote(bold("(A)"))) )
          
        }
        
        if (ii == 2468) {
          
          text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
                labels = as.expression(bquote(bold("(B)"))) )
          
        }
        
        if (ii == 2915) {
          
          text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
                labels = as.expression(bquote(bold("(C)"))) )
          
        }
  
    dev.off()
    
    }
  
    return( swissSim_2@data@values )
    
  } # Close (parallel) for loop
  
  stopCluster(cl)
  
  
  ##########################################
  ### Mean of Posterior Samples: Random Field
  ##########################################
  ii = 100
  iso_stan_2 = c( range = 1000 * swissRain_stanfit_iso$phi_2[ii],
                  shape = rep(2,4*1000)[ii],
                  variance = swissRain_stanfit_iso$phi_1[ii]^2,
                  anisoAngleDegrees = swissRain_Stan_data_iso_plot$Angle[ii],
                  anisoRatio = swissRain_Stan_data_iso_plot$Radius[ii]+1 )
  
  swissRain_Us <- swissRain 
  swissRain_Us@data <- as.data.frame( calculate_Us_iso(ii) )
  
  swissSim_2 = RandomFields::RFsimulate(
    model = RandomFields::RMmatern( 
      var = iso_stan_2['variance'][[1]], 
      scale = iso_stan_2['range'][[1]], 
      nu = iso_stan_2['shape'][[1]], 
      Aniso = matrix( ncol = 2, c( sqrt( iso_stan_2['anisoRatio'][[1]] ), 0,
                                   0, 1 / sqrt( iso_stan_2['anisoRatio'][[1]] ) ) ) %*%
        matrix( ncol = 2, c( cos(pi/180*iso_stan_2['anisoAngleDegrees'][[1]]), sin(pi/180*iso_stan_2['anisoAngleDegrees'][[1]] ), 
                             -sin(pi/180*iso_stan_2['anisoAngleDegrees'][[1]]), cos(pi/180*iso_stan_2['anisoAngleDegrees'][[1]] ) ) ) ),
    x = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,1],
    y = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,2],
    data = swissRain_Us, seed = fixed_seed )
  
  swissSim_2 <- raster::raster( swissSim_2, vals = swissSim_2$V1,
                                nrows = geostatsp::squareRaster(swissMap, raster_cells)@nrows, 
                                ncol = geostatsp::squareRaster(swissMap, raster_cells)@ncols ) 
  
  swissSim_3 <- swissSim_2 
  swissSim_3@data@values <- apply( Us_iso_conditional, 1, mean )
  
  swissMap_sf = rnaturalearth::ne_countries( country = 'Switzerland', type = "countries", 
                                             scale = "medium", returnclass = 'sf' )
  
  temp <- sf::st_transform( sf::st_centroid(swissMap_sf)$geometry, crs = swissRain@proj4string )
  temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
  
  png( filename = paste0("Figures/swissRain_Iso_Posterior_Mean.png"),
       res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
  
    mapmisc::map.new(swissBorder)
    swissCol = mapmisc::colourScale( x = swissSim_3, dec = 1, style = "fixed", rev = TRUE,
                                     col = colorRampPalette( RColorBrewer::brewer.pal(11, "Spectral") ),
                                     breaks = c( -8, -6, -4, -2, -1,
                                                  0, 1, 2, 4, 6) )
    
    plot(raster::mask(swissSim_3, swissMap), breaks = swissCol$breaks,
         col = swissCol$col, add = TRUE, legend = FALSE)
    plot(swissMap, add = TRUE)
    mapmisc::legendBreaks( 'topleft', swissCol, cex = 1.60, inset = c(0.02,0.00),
                           bty = 'n', y.intersp = 1.15 )
    
    text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 2.50,
          labels = as.expression(bquote(bold(.(swissMap$sovereignt)))) ) 
    text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
          labels = as.expression(bquote(bold("(D)"))) )
    
  dev.off()
  
  
  rm(swissRain_samples_iso)
  rm(swissRain_stanfit_iso)
  rm(swissRain_Stan_data_iso_plot)
  gc()
  
  
  
  
  ##########################################
  ###
  ### Anisotropic Model 
  ### Posterior Samples: Matern Random Field
  ###
  ##########################################
  
  swissRain_samples_aniso <- readRDS( "Results/Anisotropic_Simulation_Fit.RDS" ) 
  swissRain_stanfit_aniso <- rstan::extract( swissRain_samples_aniso )
  swissRain_Stan_data_aniso_plot <- swissRain_plotting_tibble( swissRain_stanfit_aniso )
  
  calculate_Us_aniso <- function(kk) { 
    return( chol( swissRain_stanfit_aniso$Sigma[kk,,] ) %*% 
              swissRain_stanfit_aniso$Alpha[kk,] ) 
  }
  
  
  swissAltitude_vector <- raster::extract(swissAltitude, swissRain) / 1000 
  
  cl <- makeCluster( parallel::detectCores() / 2 )
  registerDoParallel(cl)
  
  swissMap = sp::spTransform( rnaturalearth::ne_countries( country = 'Switzerland', 
                                                           type = "countries", scale = "medium" ),
                              projection(swissRain) )
  
  
  Us_aniso_conditional <- foreach ( ii = sample_idx, .combine = cbind ) %dopar% {
    
    aniso_stan_2 = c( range = 1000 * swissRain_stanfit_aniso$phi_2[ii],
                      shape = rep(2,4*1000)[ii],
                      variance = swissRain_stanfit_aniso$phi_1[ii]^2,
                      anisoAngleDegrees = swissRain_Stan_data_aniso_plot$Angle[ii],
                      anisoRatio = swissRain_Stan_data_aniso_plot$Radius[ii]+1 )
    
    swissRain_Us <- swissRain 
    swissRain_Us@data <- as.data.frame( calculate_Us_aniso(ii) )
    
    swissSim_2 = RandomFields::RFsimulate(
      model = RandomFields::RMmatern( 
        var = aniso_stan_2['variance'][[1]], 
        scale = aniso_stan_2['range'][[1]], 
        nu = aniso_stan_2['shape'][[1]], 
        Aniso = matrix( ncol = 2, c( sqrt( aniso_stan_2['anisoRatio'][[1]] ), 0,
                                     0, 1 / sqrt( aniso_stan_2['anisoRatio'][[1]] ) ) ) %*%
          matrix( ncol = 2, c( cos(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]]), sin(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]] ), 
                               -sin(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]]), cos(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]] ) ) ) ),
      x = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,1],
      y = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,2],
      data = swissRain_Us, seed = fixed_seed )
    
    swissSim_2 <- raster::raster( swissSim_2, vals = swissSim_2$V1,
                                  nrows = geostatsp::squareRaster(swissMap, raster_cells)@nrows, 
                                  ncol = geostatsp::squareRaster(swissMap, raster_cells)@ncols ) 
    
    swissMap_sf = rnaturalearth::ne_countries( country = 'Switzerland', type = "countries",
                                               scale = "medium", returnclass = 'sf' )
  
    temp <- sf::st_transform( sf::st_centroid(swissMap_sf)$geometry, crs = swissRain@proj4string )
    temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
    
    if (ii %in% c(1598,1784,1968) ) {
      
      png( filename = paste0("Figures/swissRain_Aniso_Posterior_Sample_",ii,".png"),
           res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
    
      mapmisc::map.new(swissBorder)
      swissCol = mapmisc::colourScale( x = swissSim_2, dec = 1, style = "fixed", rev = TRUE,
                                       col = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")),
                                       breaks = c( -8, -6, -4, -2, -1,
                                                    0, 1, 2, 4, 6) )
      plot(raster::mask(swissSim_2, swissMap), breaks = swissCol$breaks,
           col = swissCol$col, add = TRUE, legend = FALSE)
      plot(swissMap, add = TRUE)
      mapmisc::legendBreaks( 'topleft', swissCol, cex = 1.60, inset = c(0.02,0.00),
                             bty = 'n', y.intersp = 1.15 )
      text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 2.50,
            labels = as.expression(bquote(bold(.(swissMap$sovereignt)))) ) 
    
      if (ii == 1598) {
        
        text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
              labels = as.expression(bquote(bold("(A)"))) )
        
      }
      
      if (ii == 1784) {
        
        text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
              labels = as.expression(bquote(bold("(B)"))) )
        
      }
      
      if (ii == 1968) {
        
        text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
              labels = as.expression(bquote(bold("(C)"))) )
        
      }
      
  
    dev.off()
    
    }
    
    return( swissSim_2@data@values )
    
  } # Close (parallel) for loop
  
  stopCluster(cl)
    
    
  ##########################################
  ### Mean of Posterior Samples: Random Field
  ##########################################
  ii = 100
  aniso_stan_2 = c( range = 1000 * swissRain_stanfit_aniso$phi_2[ii],
                    shape = rep(2,4*1000)[ii],
                    variance = swissRain_stanfit_aniso$phi_1[ii]^2,
                    anisoAngleDegrees = swissRain_Stan_data_aniso_plot$Angle[ii],
                    anisoRatio = swissRain_Stan_data_aniso_plot$Radius[ii]+1 )
  
  swissRain_Us <- swissRain 
  swissRain_Us@data <- as.data.frame( calculate_Us_aniso(ii) )
  
  swissSim_2 = RandomFields::RFsimulate(
    model = RandomFields::RMmatern( 
      var = aniso_stan_2['variance'][[1]], 
      scale = aniso_stan_2['range'][[1]], 
      nu = aniso_stan_2['shape'][[1]], 
      Aniso = matrix( ncol = 2, c( sqrt( aniso_stan_2['anisoRatio'][[1]] ), 0,
                                   0, 1 / sqrt( aniso_stan_2['anisoRatio'][[1]] ) ) ) %*%
        matrix( ncol = 2, c( cos(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]]), sin(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]] ), 
                             -sin(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]]), cos(pi/180*aniso_stan_2['anisoAngleDegrees'][[1]] ) ) ) ),
    x = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,1],
    y = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,2],
    data = swissRain_Us, seed = fixed_seed )
  
  swissSim_2 <- raster::raster( swissSim_2, vals = swissSim_2$V1,
                                nrows = geostatsp::squareRaster(swissMap, raster_cells)@nrows, 
                                ncol = geostatsp::squareRaster(swissMap, raster_cells)@ncols ) 
  
  swissMap_sf = rnaturalearth::ne_countries( country = 'Switzerland', type = "countries", 
                                             scale = "medium", returnclass = 'sf' )
  
  temp <- sf::st_transform( sf::st_centroid(swissMap_sf)$geometry, crs = swissRain@proj4string )
  temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
  
  swissSim_3 <- swissSim_2 
  swissSim_3@data@values <- apply( Us_aniso_conditional, 1, mean )
  
  
  png( filename = paste0("Figures/swissRain_Aniso_Posterior_Mean.png"),
       res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
  
    mapmisc::map.new(swissBorder)
    swissCol = mapmisc::colourScale( x = swissSim_3, dec = 1, style = "fixed", rev = TRUE,
                                     col = colorRampPalette( RColorBrewer::brewer.pal(11, "Spectral") ),
                                     breaks = c( -8, -6, -4, -2, -1,
                                                 0, 1, 2, 4, 6 ) )
    
    plot(raster::mask(swissSim_3, swissMap), breaks = swissCol$breaks,
         col = swissCol$col, add = TRUE, legend = FALSE)
    plot(swissMap, add = TRUE)
    mapmisc::legendBreaks( 'topleft', swissCol, cex = 1.60, inset = c(0.02,0.00),
                           bty = 'n', y.intersp = 1.15 )
    
    text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 2.50,
          labels = as.expression(bquote(bold(.(swissMap$sovereignt)))) ) 
    text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
          labels = as.expression(bquote(bold("(D)"))) )
  
  dev.off()
  
  
  rm(swissRain_samples_aniso)
  rm(swissRain_stanfit_aniso)
  rm(swissRain_Stan_data_aniso_plot)
  gc()
  
  
  
  
  ##########################################
  ###
  ### Posterior Samples: Matern Random Field
  ### Switzerland Rainfall 
  ### 
  ##########################################
  
  swissRain_samples <- readRDS( "Results/Switzerland_Rainfall_Fit.RDS" )
  swissRain_stanfit <- rstan::extract( swissRain_samples )
  swissRain_Stan_data_real_plot <- swissRain_plotting_tibble( swissRain_stanfit )
  
  calculate_Us_real <- function(kk) { 
    return( chol( swissRain_stanfit$Sigma[kk,,] ) %*% 
              swissRain_stanfit$Alpha[kk,] ) 
  }
  
    
  swissAltitude_vector <- raster::extract(swissAltitude, swissRain) / 1000 
  
  cl <- makeCluster( parallel::detectCores() / 2 )
  registerDoParallel(cl)
  
  swissMap = sp::spTransform( rnaturalearth::ne_countries(  country = 'Switzerland', 
                                                            type = "countries", scale = "medium" ),
                              projection(swissRain) )
  
  num_Us_samples <- num_samples
  sample_idx <- sample(1:4000, num_samples)
  
  
  Us_real_conditional <- foreach ( ii = sample_idx, .combine = cbind ) %dopar% {
    
    real_stan_2 = c( range = 1000 * swissRain_stanfit$phi_2[ii],
                     shape = rep(2,4*1000)[ii],
                     variance = swissRain_stanfit$phi_1[ii]^2,
                     anisoAngleDegrees = swissRain_Stan_data_real_plot$Angle[ii],
                     anisoRatio = swissRain_Stan_data_real_plot$Radius[ii]+1 )
    
    swissRain_Us <- swissRain 
    swissRain_Us@data <- as.data.frame( calculate_Us_real(ii) )
    
  
    swissSim_2 = RandomFields::RFsimulate(
      model = RandomFields::RMmatern( 
        var = real_stan_2['variance'][[1]], 
        scale = real_stan_2['range'][[1]], 
        nu = real_stan_2['shape'][[1]], 
        Aniso = matrix( ncol = 2, c( sqrt( real_stan_2['anisoRatio'][[1]] ), 0,
                                     0, 1 / sqrt( real_stan_2['anisoRatio'][[1]] ) ) ) %*%
          matrix( ncol = 2, c( cos(pi/180*real_stan_2['anisoAngleDegrees'][[1]]), sin(pi/180*real_stan_2['anisoAngleDegrees'][[1]] ), 
                               -sin(pi/180*real_stan_2['anisoAngleDegrees'][[1]]), cos(pi/180*real_stan_2['anisoAngleDegrees'][[1]] ) ) ) ),
      x = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,1],
      y = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,2],
      data = swissRain_Us, seed = fixed_seed )
    
    swissSim_2 <- raster::raster( swissSim_2, vals = swissSim_2$V1,
                                  nrows = geostatsp::squareRaster(swissMap, raster_cells)@nrows, 
                                  ncol = geostatsp::squareRaster(swissMap, raster_cells)@ncols ) 
    
    swissMap_sf = rnaturalearth::ne_countries( country = 'Switzerland', type = "countries",
                                               scale = "medium", returnclass = 'sf' )
  
    temp <- sf::st_transform( sf::st_centroid(swissMap_sf)$geometry, crs = swissRain@proj4string )
    temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
    
    if ( ii %in% c(3774,3826,3903) ) {
      
      png( filename = paste0("Figures/swissRain_Posterior_Sample_",ii,".png"),
           res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
    
        mapmisc::map.new(swissBorder)
        swissCol = mapmisc::colourScale( x = swissSim_2, dec = 1, style = "fixed", rev = TRUE,
                                         col = colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")),
                                         breaks = c( -8, -6, -4, -2, -1, 
                                                     0, 1, 2, 4, 6 ) )
  
        plot(raster::mask(swissSim_2, swissMap), breaks = swissCol$breaks,
             col = swissCol$col, add = TRUE, legend = FALSE)
        plot(swissMap, add = TRUE)
        mapmisc::legendBreaks( 'topleft', swissCol, cex = 1.60, inset = c(0.02,0.00),
                               bty = 'n', y.intersp = 1.15 )
        
        text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 2.50,
              labels = as.expression(bquote(bold(.(swissMap$sovereignt)))) ) 
        
        if (ii == 3774) {
          
          text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
                labels = as.expression(bquote(bold("(A)"))) )
          
        }
        
        if (ii == 3826) {
          
          text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
                labels = as.expression(bquote(bold("(B)"))) )
          
        }
        
        if (ii == 3903) {
          
          text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
                labels = as.expression(bquote(bold("(C)"))) )
          
        }
        
      
      dev.off()
    
    }
  
    return( swissSim_2@data@values )
    
  } # Close (parallel) for loop
  
  stopCluster(cl)
  
    
  ##########################################
  ### Mean of Posterior Samples: Random Field
  ##########################################
  ii = 500
  real_stan_2 = c( range = 1000 * swissRain_stanfit$phi_2[ii],
                   shape = rep(2,4*1000)[ii],
                   variance = swissRain_stanfit$phi_1[ii]^2,
                   anisoAngleDegrees = swissRain_Stan_data_real_plot$Angle[ii],
                   anisoRatio = swissRain_Stan_data_real_plot$Radius[ii]+1 )
  
  swissRain_Us <- swissRain
  swissRain_Us@data <- as.data.frame( calculate_Us_real(ii) )
  
  swissSim_2 = RandomFields::RFsimulate(
    model = RandomFields::RMmatern( 
      var = real_stan_2['variance'][[1]], 
      scale = real_stan_2['range'][[1]], 
      nu = real_stan_2['shape'][[1]], 
      Aniso = matrix( ncol = 2, c( sqrt( real_stan_2['anisoRatio'][[1]] ), 0,
                                   0, 1 / sqrt( real_stan_2['anisoRatio'][[1]] ) ) ) %*%
        matrix( ncol = 2, c( cos(pi/180*real_stan_2['anisoAngleDegrees'][[1]]), sin(pi/180*real_stan_2['anisoAngleDegrees'][[1]] ), 
                             -sin(pi/180*real_stan_2['anisoAngleDegrees'][[1]]), cos(pi/180*real_stan_2['anisoAngleDegrees'][[1]] ) ) ) ),
    x = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,1],
    y = raster::rasterToPoints( geostatsp::squareRaster(swissMap, raster_cells) )[,2],
    data = swissRain_Us, seed = fixed_seed )
  
  swissSim_2 <- raster::raster( swissSim_2, vals = swissSim_2$V1,
                                nrows = geostatsp::squareRaster(swissMap, raster_cells)@nrows, 
                                ncol = geostatsp::squareRaster(swissMap, raster_cells)@ncols ) 
  
  swissMap_sf = rnaturalearth::ne_countries( country = 'Switzerland', type = "countries", 
                                             scale = "medium", returnclass = 'sf' )
  
  temp <- sf::st_transform( sf::st_centroid(swissMap_sf)$geometry, crs = swissRain@proj4string )
  temp_matrix <- matrix( unlist(temp), ncol = 2, byrow = TRUE )
  
  swissSim_3 <- swissSim_2 
  swissSim_3@data@values <- apply( Us_real_conditional, 1, mean )
  
  
  png( filename = paste0("Figures/swissRain_Posterior_Mean.png"),
       res = 300, width = 1.5 * 1928, height = 1.5 * 1514, unit = "px" )
  
    mapmisc::map.new(swissBorder)
    swissCol = mapmisc::colourScale( x = swissSim_3, dec = 1, style = "fixed", rev = TRUE,
                                     col = colorRampPalette( RColorBrewer::brewer.pal(11, "Spectral") ),
                                     breaks = c( -8, -6, -4, -2, -1,
                                                 0, 1, 2, 4, 6 ) )
    
    plot(raster::mask(swissSim_3, swissMap), breaks = swissCol$breaks,
         col = swissCol$col, add = TRUE, legend = FALSE)
    plot(swissMap, add = TRUE)
    mapmisc::legendBreaks( 'topleft', swissCol, cex = 1.60, inset = c(0.02,0.00),
                           bty = 'n', y.intersp = 1.15 )
    
    text( x = temp_matrix[,1], y = temp_matrix[,2], cex = 2.50,
          labels = as.expression(bquote(bold(.(swissMap$sovereignt)))) ) 
    text( x = temp_matrix[1,1]-110*1000 , y = temp_matrix[1,2]+110*1000, cex = 2.50,
          labels = as.expression(bquote(bold("(D)"))) )
    
  dev.off()
  
  
  
  rm(swissRain_samples)
  rm(swissRain_stanfit)
  rm(swissRain_Stan_data_real_plot)
  gc()