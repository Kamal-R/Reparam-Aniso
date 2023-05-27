This repository contains the files used the paper "A Parameter Transformation of the Anisotropic Matern Covariance Function."" All R scripts assume they are in R's working directory. In RStudio, this can be set using 'Session -> Set Working Directory -> To Source File Location' when the R script being run is open and selected in Rstudio. 

This repository contains two folders -- 'swissRain' and 'mercurySoil.' The 'swissRain' folder contains the files used to produce the results in section 4.2 (and the Appendix). The 'mercurySoil' folder contains the files used to produce the results in Section 4.3.

## The SwissRain Folder

The swissRain folder contains eight files and two sub-folders. They are described below. The first script to run in is 'Fit_Anisotropic_Matern_Model.R' as it fits (and optionally saves) Stan models to the Results folder.

### Files

-- Anisotropic_Matern_Model_swissRain.stan: The Stan file used to fit the spatial model with the re-parameterized Matern covariance function to the Switzerland rainfall data.

-- Fit_Anisotropic_Matern_Model.R: Run this R script to fit the re-parameterized Matern model to simulated isotropic and simulated anisotropic data, as well as the swissRain data. 

-- Univariate_Densities.R: Produces the summary tables, traceplots, and density plots seen in Section 4 and the Appendix. 

-- Isotropic_Simulation_LCD_Plots.R, Anisotropic_Simulation_LCD_Plots.R, Switzerland_Rainfall_LCD_Plots.R: Produces the bivariate density plots for the standard and transformed Matern parameters, for the isotropic, anisotropic, and Switzerland rainfall data (respectively) in Section 4. 

-- Simulate_Random_Fields.R: Produces plots of the posterior Matern random field for the Switzerland rainfall data (in Section 4) and for the simulated isotropic and anisotropic data (in the Appendix). 

-- Helper_Functions.R: Functions that pre-process the posterior samples from Stan into a format that is more suitable for plotting. 

### Sub-Folders

Results: The summary and plotting scripts assume the Stan models produced by Anisotropic_Matern_Model_swissRain.stan have been saved to this sub-folder. Optionally, the log-concave density estimates can also be saved to this sub-folder.

Figures: The figures produced by the various plotting scripts are saved to this sub-folder.


## The mercurySoil Folder

The mercurySoil folder contains four files and three folders. They are described below. The first script to run is 'Fit_Anisotropic_Matern_Model.R' as it fits and saves the Stan model to the Results folder. 

### Files

-- Anisotropic_Matern_Model_MercurySoil.stan: The Stan file used to fit the spatial model with the re-parameterized Matern covariance function to the mercury soil data. 

-- Fit_Anisotropic_Matern_Model.R: Run this R script to the fit the re-parameterized Matern model to the mercury soil data. 

-- Summary_Results.R: Produces the tables and figures (including the bivariate density plots) for the mercury soil data seen in Section 4.3. 

-- Simulate_Random_Fields.R: Produces the plots of the posterior Matern random field for the mercury soil data (seen in Section 4).


### Sub-Folders

-- Results: The 'Summary_Results.R' and 'Simulate_Random_Fields.R' scripts assume the results from Stan are saved to this sub-folder. 

-- Figures: The plots produces by both scripts are saved to this sub-folder. 

-- Data: The mercury soil data and covariate data are available in the 'Data' sub-folder.