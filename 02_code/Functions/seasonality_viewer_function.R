#--------------------------------------------------------------------------------------------------#
#---- seasonality_viewer() function -----#
#----------------------------------------#

#--- README ---#
#--------------#
##' The function below allows users to input seasonality parameters and generate a quick plot. The
##' function works by loading the rainfall() function from the mosquito_biology.cpp file and using
##' it to calculate the daily rainfall given a set of seasonal fourier parameters (g0, g, h) and a
##' number of years to produce daily rainfall values for. 
##' 
##' The function also uses the peak_season_offset() function to calculate and overlay the day on
##' which peak rainfall occurs each year. Once the rainfall and rainfall peak have been calculated,
##' the function plots rainfall over time to let you look at the meteorological monstrosity you've
##' managed to cook up for malaria this time!

# Load in the malariasimulation and Rcpp packages:
library(malariasimulation); library(Rcpp)

# Source in the rainfall() function (C++) - edit to wherever the file is saved)
sourceCpp("C:/Users/kem22/mosquito_biology.cpp")

# seasonality_viewer() function:
# Create the seasonality_viewer) function:
seasonality_viewer <- function(g0_in, g_in, h_in, years) {
  
  # Generate the simulation parameter list with the specified seasonality parameters
  simparams <- get_parameters(
    list(
      human_population = 1000,
      model_seasonality = TRUE,
      g0 = g0_in,
      g = g_in,
      h = h_in
    )
  )
  
  # Use the peak_season_offset() to calculate the yearly offset for peak rainfall season
  peak <- peak_season_offset(parameters = simparams)
  
  # Open a vector in which to store the rainfall
  daily_rain <- vector()
  
  # The rainfall() function:
  for(i in 1:(years*365)) {
    
    # Calculate the daily rainfall for each time step:
    daily_rain[i] <- rainfall(
      t = i,
      g0 = simparams$g0,
      g =  simparams$g,
      h = simparams$h,
      floor = simparams$rainfall_floor)
  }
  
  # Create a dataframe for plotting:
  rainfall.df <- data.frame(timestep = 1:(365 * years),
                            daily_rain = daily_rain)
  
  # Define a colour palette for plotting:
  plot_cols <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7","#F0E442", "#0072B2", "#D55E00")
  
  # Open a new plotting window,set the margins and plot the mosquito population
  plot.new(); par(new = TRUE, mar = c(4, 4, 1, 1))
  plot(x = rainfall.df$timestep, y = rainfall.df$daily_rain,
       xlab = "Time (days)", ylab = "Daily Rainfall",
       ylim = c(0, max(rainfall.df$daily_rain)*1.1),
       type = "l", lwd = 2, cex = 0.8,
       xaxs = "i", yaxs = "i",
       col = plot_cols[3])
  
  # Add grinlines
  #grid(lty = 1, col = "darkgrey", nx = 11, ny = NULL, lwd = 0.5); box()
  
  # Overlay the calculated seasonal peak
  abline(v =  seq(from = 0, to = (years-1), by = 1) * 365 + peak,
         lty = 2, lwd = 2.5, col = plot_cols[4])
  
}

# Example:
seasonality_viewer(g0_in = 0.284596,
                   g_in = c(-0.317878,-0.0017527,0.116455),
                   h_in = c(-0.331361,0.293128,-0.0617547), 
                   years = 3)

#--------------------------------------------------------------------------------------------------#












