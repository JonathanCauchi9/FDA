# Load necessary libraries
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("fda")
install.packages("imputeTS")

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(fda)
library(imputeTS)

# Load and preprocess the data
Radiation_Data <- read_excel("C:/Users/jonat/OneDrive/Desktop/Thesis/Data Files/Global Radiation.xlsx")
View(Radiation_Data)

# 2) Reshape from 'long' to 'wide'
Radiation_Data <- Radiation_Data %>%
  pivot_wider(
    # 'Station' is the column whose *unique values* become new column headers
    names_from = Station,
    # 'Daily mean temperature...' is the measurement that fills the cells
    values_from = `Global radiation (in J/cm2)`,
    # Keep 'Date' as our main row-identifier
    id_cols = Date
  )

# 3) Inspect the result
View(Radiation_Data)
Radiation_Data$Date <- as.Date(
  as.character(Radiation_Data$Date), 
  format = "%Y%m%d"
)
time_points_radiation <- as.numeric(
  Radiation_Data$Date - min(Radiation_Data$Date)
)

time_points_radiation <- as.numeric(Radiation_Data$Date - min(Radiation_Data$Date))
time_points_radiation
min_time <- min(time_points_radiation)
max_time <- max(time_points_radiation)
min_time
max_time

Radiation_Data <- Radiation_Data[,-c(1,5,26)]
View(Radiation_Data)
Radiation_numeric<-as.matrix(Radiation_Data)
View(Radiation_numeric)

na_by_col <- colSums(is.na(Radiation_Data))
print("Number of NA values per column:")
print(na_by_col)

#Imputation through spline smoothing
Radiation_Data_imp <- Radiation_Data %>%
  mutate(
    across(
      where(is.numeric),              # only numeric columns
      ~ na_interpolation(.x, option = "spline")  
      # or option = "linear" if you prefer a linear fill
    )
  )
View(Radiation_Data_imp)
View(Radiation_Data)

#Creating the Fourier Basis:
n_basis_fourier <- 7
basis_fourier <- create.fourier.basis(
  rangeval = c(min_time, max_time),
  nbasis   = n_basis_fourier
)

# Evaluate Fourier basis functions at the time points
PHI_fourier <- eval.basis(time_points_radiation, basis_fourier)
dim(PHI_fourier)

Radiation_Data_imp_numeric<- as.matrix(Radiation_Data_imp)
n_stations<-ncol(Radiation_Data_imp_numeric)
n_stations
n_time <- nrow(Radiation_Data_imp_numeric)            
n_time

# Checking the residuals to asses whether an AR approach is viable
n_order <- 4                             
Sec_derivative <- int2Lfd(max(0, n_order - 2))
# Initialize matrices to store residuals and fitted values for each station
residuals_matrix_radiation <- matrix(NA, nrow = n_time, ncol = n_stations)
fitted_values_matrix_radiation <- matrix(NA, nrow = n_time, ncol = n_stations)

# Loop over each station (column)
for (station in 1:n_stations) {
  cat(paste("Performing OLS for station", colnames(Radiation_Data_imp_numeric)[station], "\n"))
  
  # Extract data for the current station
  y <- Radiation_Data_imp_numeric[, station]
  
  # Construct the basis function matrix without regularization
  Phi <- eval.basis(time_points_radiation, basis_fourier)  # Evaluate basis at time points
  
  # Perform OLS estimation using the normal equation: c_hat = (Phi'Phi)^(-1) Phi'Y
  c_hat <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% y
  
  # Compute the fitted values: Y_hat = Phi * c_hat
  fitted_values <- Phi %*% c_hat
  
  # Calculate residuals: Observed - Fitted values
  residuals <- y - fitted_values
  
  # Store results in matrices
  residuals_matrix_radiation[, station] <- residuals  #Store residuals
  fitted_values_matrix_radiation[, station] <- fitted_values  #Store fitted values
}

# View residuals/fitted values for inspection
View(residuals_matrix_radiation)  # Residuals matrix: rows = time points, columns = stations
View(fitted_values_matrix_radiation)
View(Radiation_Data_imp_numeric)


#Checking whether there actually is auto-correlation between the residuals through numerous checks

#PACF Checks
# Load necessary library
library(forecast)

# Ensure residuals_matrix and Radiation_Data_imp_numeric are properly defined
n_stations <- ncol(Radiation_Data_imp_numeric)
n_stations
# Set the seed for reproducibility
set.seed(22)  # For consistent random sampling

# Randomly sample 6 stations
sampled_stations <- sample(1:n_stations, 6)

# Set up plotting window: one PACF plot per station
par(mfrow = c(3, 2))  # 3 rows, 2 columns for 6 plots

# Loop through each sampled station and plot the PACF
for (station in sampled_stations) {
  # Calculate the PACF for the current station's residuals
  pacf_result <- pacf(residuals_matrix_radiation[, station], plot = FALSE, lag.max = 50)
  
  # Extract the lag with the highest PACF value (excluding lag 0)
  lag_values <- pacf_result$lag
  pacf_values <- pacf_result$acf
  max_pacf_index <- which.max(abs(pacf_values))  # Find the index of the highest PACF value
  max_lag <- lag_values[max_pacf_index]
  max_pacf <- pacf_values[max_pacf_index]
  
  # Display the lag and PACF value
  cat(paste("station", station, "Highest PACF Value at Lag", max_lag, ":", round(max_pacf, 3), "\n"))
  
  # Plot the PACF for the current station
  main_title <- paste("PACF of Residuals for station", station, 
                      "\nMax PACF at Lag", max_lag, "=", round(max_pacf, 3))
  pacf(residuals_matrix_radiation[, station], main = main_title, lag.max = 50)
}

#Ljung Box Test


# Number of lags to check for autocorrelation
lag_k1 <- 1

# Loop through each station
for (station in 1:n_stations) {
  # Perform the Ljung-Box test on the residuals of the current station
  lb_test <- Box.test(residuals_matrix_radiation[, station], lag = lag_k1, type = "Ljung-Box")
  
  # Print the result for the current station
  cat(paste("Ljung-Box Test for station", station, "\n"))
  print(lb_test)
}
#P values are smaller than 0.05 for all stations so we accept H1: The auto-correlations up to lag 1 are not zero


#N.B.
#The PACF cuts off after a small number of lags (e.g., 1 or 2) because once you account for the correlation with the immediate past value, there's little to no direct correlation at further lags.
#However, the Ljung-Box test still detects significant autocorrelation across several lags (even if individual lags are small) because it looks at the cumulative autocorrelation over a range of lags.

# Reset plotting layout
par(mfrow = c(1, 1))

#Plotting a scatter plot of residuals at t+1 against residual at t to check for any correlation 
n_stations <- ncol(residuals_matrix_radiation)
n_time <- nrow(residuals_matrix_radiation)
n_stations
# Check autocorrelation in residuals for each station
for (station in 1:n_stations) {
  cat(paste("Checking autocorrelation for station", station, "\n"))
  
  # Scatter plot of residuals vs. lagged residuals (to visually inspect autocorrelation)
  plot(residuals_matrix_radiation[1:(n_time - 1), station], residuals_matrix_radiation[2:n_time, station], 
       xlab = "Residual(t)", ylab = "Residual(t+1)", 
       main = paste("Scatter plot of residuals for station", station))
  
}

#So fit an AR(1) Model to the residuals 

#Fitting the AR(1) model and estimating the variance covariance matrix of the residuals

# Load the necessary library for AR(1) model fitting
library(forecast)

# Assuming residuals_matrix is already prepared with residuals for each station (each column)
n_stations <- ncol(residuals_matrix_radiation)
n_time <- nrow(residuals_matrix_radiation)
n_stations
# Initialize storage for AR(1) parameters and the covariance matrix
rho_values <- numeric(n_stations)    # Autocorrelation parameter for each station
sigma2_values <- numeric(n_stations) # Variance of noise for each station

# Initialize an empty list to store the covariance matrices for each station
cov_matrix <- vector("list", n_stations)  # List to store the covariance matrix for each station

# Fit AR(1) model to each station's residuals and estimate parameters
wtvec_list <- list()  

for (station in 1:n_stations) {
  tryCatch({
    # Fit AR(1) model to the residuals of the current station
    ar1_model <- arima(residuals_matrix_radiation[, station], order = c(1, 0, 0))  # AR(1) model
    
    # Extract parameters from the AR(1) model
    rho_values[station] <- ar1_model$coef[1]  # Autocorrelation parameter (rho)
    sigma2_values[station] <- ar1_model$sigma2  # Variance of noise (sigma^2)
    
    # Print the estimated parameters for each station
    cat(paste("station", station, "AR(1) Model Parameters:\n"))
    cat("Autocorrelation parameter (rho):", rho_values[station], "\n")
    cat("Variance of noise (sigma^2):", sigma2_values[station], "\n")
    
    # Initialize the covariance matrix for this station
    station_cov_matrix <- matrix(NA, nrow = n_time, ncol = n_time)  # Covariance matrix for residuals
    
    # Fill in the covariance matrix using the AR(1) covariance structure
    for (t1 in 1:n_time) {
      for (t2 in 1:n_time) {
        # Compute the lag (h) between the two time points
        h <- abs(t1 - t2)
        
        # Calculate the covariance between the two time points
        station_cov_matrix[t1, t2] <- sigma2_values[station] * (rho_values[station]^h)
      }
    }
    
    # Store the covariance matrix for the current station
    cov_matrix[[station]] <- station_cov_matrix
    
    # Store the inverse of the covariance matrix as the weight vector
    wtvec_list[[station]] <- solve(station_cov_matrix)
    
    # Optionally, you can print the covariance matrix for inspection
    cat(paste("Covariance Matrix for station", station, ":\n"))
    print(station_cov_matrix)
  }, error = function(e) {
    # If there's an error, print a message and move on to the next station
    cat(paste("Error fitting AR(1) model for station", station, ": ", e$message, "\n"))
  })
}

parameter_estimates <- data.frame(
  station = 1:n_stations,
  Rho = rho_values,
  Sigma2 = sigma2_values
)
print(parameter_estimates)

lambda_range_weighted <- seq(1e-8, 11000, length.out = 50)  # Define the range of lambda values

# Initialize a data frame to store optimal lambda results for each station
optimal_lambdas_radiation_weighted <- data.frame(
  Station = colnames(Radiation_Data_imp_numeric),
  Optimal_lambda_weighted = NA,
  Min_GCV_weighted = NA
)

# Loop through each station
for (station in colnames(Radiation_Data_imp_numeric)) {
  cat(paste("Calculating metrics for station (weighted):", station, "\n"))
  
  # Extract the data and weight vector for the current station
  y <- Radiation_Data_imp_numeric[, station]
  wtvec <- wtvec_list[[station]]
  
  # Initialize storage for GCV values
  gcvsave_weighted <- numeric(length(lambda_range_weighted))
  
  # Loop through lambda values
  for (i in seq_along(lambda_range_weighted)) {
    lambda <- lambda_range_weighted[i]  
    fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = lambda)  
    
    # Perform smoothing
    smooth_result <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj, wtvec = wtvec)
    
    # Extract GCV values
    gcvsave_weighted[i] <- smooth_result$gcv  
  }
  
  # Identify the optimal lambda corresponding to the minimum GCV
  min_gcv_index <- which.min(gcvsave_weighted)
  optimal_lambda_weighted <- lambda_range_weighted[min_gcv_index]
  min_gcv_weighted <- gcvsave_weighted[min_gcv_index]
  
  # Store the optimal GCV results for this station
  optimal_lambdas_radiation_weighted[optimal_lambdas_radiation_weighted$Station == station, ] <- list(
    Station = station,
    Optimal_lambda_weighted = optimal_lambda_weighted,
    Min_GCV_weighted = min_gcv_weighted
  )
  
  # Plot GCV curve for the current station
  plot(lambda_range_weighted, gcvsave_weighted, type = "b", col = "blue", lwd = 2, pch = 19,
       xlab = "Lambda", ylab = "GCV", 
       main = paste("GCV Curve for station", station))
  
  # Highlight the optimal lambda value on the plot
  abline(v = optimal_lambda_weighted, col = "red", lty = 2, lwd = 1.5)
}

print(optimal_lambdas_radiation_weighted)

fitted_values_matrix_radiation_weight <- matrix(NA, nrow = n_time, ncol = n_stations)
# Access the first matrix in the list
# Access the first row of the first matrix
first_row <- cov_matrix[[1]][1,]
print(first_row)

# Smooth all stations with their respective optimal lambdas
smoothed_stations_weighted <- list()

for (station in 1:n_stations) {
  station_name <- colnames(Radiation_Data_imp_numeric)[station]
  cat(paste("Smoothing station:", station_name, "\n"))
  
  # Retrieve optimal lambda for the current station
  optimal_lambda <- optimal_lambdas_radiation_weighted$Optimal_lambda[station]
  
  # Extract radiation data for the current station
  y <- Radiation_Data_imp_numeric[, station]
  
  # Retrieve corresponding weight vector
  wtvec <- wtvec_list[[station_name]]
  
  # Create functional parameter object with optimal lambda and weights
  fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = optimal_lambda)
  
  # Perform smoothing with weights
  smooth_result_weight <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj, wtvec = wtvec)
  
  # Store the smoothed functional data object
  smoothed_stations_weighted[[station_name]] <- smooth_result_weight$fd
  
  # Extract fitted values for the station
  fitted_values_matrix_radiation_weight[, station] <- eval.fd(time_points_radiation, smooth_result_weight$fd)
}

# Convert the fitted values matrix to a dataframe for easier inspection
fitted_values_df <- as.data.frame(fitted_values_matrix_radiation_weight)
colnames(fitted_values_df) <- colnames(Radiation_Data_imp_numeric)

# View the fitted values
View(fitted_values_df)

# Reset plotting layout to default
par(mfrow = c(1, 1))


# Combine smoothed functions into a single functional data object for plotting
combined_fd_radiation <- fd(
  coef = do.call(cbind, lapply(smoothed_stations_weighted, function(fd) fd$coefs)),
  basisobj = basis_fourier
)

# Plot all smoothed functions
plot(
  combined_fd_radiation,
  xlab = "Time", ylab = "Smoothed Value",
  main = "Smoothed Functions for all stations"
)

#Performing Functional Principal Component Analysis (FPCA) 
weighted_fpca_results_radiation <- pca.fd(fdobj = combined_fd_radiation, nharm = 2)
weighted_fpca_results_radiation$harmonics #a functional data object for the harmonics or eigen-functions
weighted_fpca_results_radiation$scores #s matrix of scores on the principal components or harmonics
weighted_fpca_results_radiation$values #the complete set of eigenvalues
weighted_fpca_results_radiation$varprop #a vector giving the proportion of variance explained by each eigen-function
weighted_fpca_results_radiation$meanfd #a functional data object giving the mean function
plot.pca.fd(weighted_fpca_results_radiation)
weighted_fpca_results_radiation$varprop
initial_fpca_results_radiation$varprop

par(mfrow = c(1, 1))
#Plotting the first harmonic
plot(weighted_fpca_results_radiation$harmonics[1], main = "First Harmonic (PC1)", xlab = "Date", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

#Plotting the second harmonic
plot(weighted_fpca_results_radiation$harmonics[2], main = "Second Harmonic (PC2)", xlab = "Date", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

# Compute the first derivatives of the smoothed functions
combined_fd_deriv_rad <- deriv.fd(combined_fd_radiation, Lfdobj = 1)

# Plot the first derivatives for all assets in one plot
plot(combined_fd_deriv_rad, 
     xlab = "Time (Days)", 
     ylab = "First Derivative", 
     main = "First Derivative of Smoothed Asset Functions")

# Alternatively, to plot derivatives for a subset of assets (e.g., first 3 assets)
selected_indices <- 1:3  # Change this to select different assets

# Set up a multi-plot layout
par(mfrow = c(ceiling(length(selected_indices) / 2), 2))
for (i in selected_indices) {
  asset_name <- colnames(Prices_data_numeric)[i]
  plot(combined_fd_deriv_rad[i], 
       xlab = "Time (Days)", 
       ylab = "First Derivative", 
       main = paste("Derivative of", asset_name))
}
# Reset plotting layout
par(mfrow = c(1, 1))
