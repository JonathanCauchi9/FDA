install.packages("refund")
install.packages("tidyr")
install.packages("tseries")
install.packages("urca")
install.packages("FinTS")
install.packages("rugarch")
library(rugarch)
library(ggplot2)
library(FinTS)
library(dplyr)
library(fda)
library(refund)
library(tidyr)
library(ggplot2)
library(forecast)
View(pricesdata)
# Extract the SPX.L data along with the Date column
spx_data <- pricesdata %>% 
  dplyr::select(Date, SPX.L)

# View the extracted data (this opens the data viewer in RStudio)
View(spx_data)

Prices<-pricesdata
View(Prices)
dev.off()

#Data Preparation 
#Extracting the Date column and converting it to a numeric sequence for time indexing
Prices$Date <- as.Date(Prices$Date)
time_points <- as.numeric(Prices$Date - min(Prices$Date)) 
time_points
length(time_points)
max(time_points)

#Removing the Date column to isolate the Prices data
Prices_data <- Prices %>% dplyr::select(-Date)
View(Prices_data)



# Checking for Periodicity
set.seed(16)  
sample_size <- 6  
sampled_assets <- sample(colnames(Prices_data), sample_size)

par(mfrow = c(ceiling(sample_size / 2), 2))  

# Loop through the sampled assets to plot the ACF
for (asset in sampled_assets) {
  # Title for the plot
  main_title <- paste("ACF of", asset, "Prices")
  
  # Plot the ACF for the asset
  acf(Prices_data[[asset]], main = main_title, lag.max = 50)
}

# Reset plotting layout
par(mfrow = c(1, 1))



#Step 2: Creating the B-spline Basis
#Define the time range
time_points
min_time <- min(time_points)
max_time <- max(time_points)
min_time
max_time

#Define knots and basis parameters
knots <- c(seq(min_time,max_time,87)) #Set knots every 87 days
knots
n_knots <- length(knots)
n_knots
n_order <- 4                              #Order of the basis (cubic B-spline)
n_basis <- n_knots + n_order - 2          #Total number of basis functions
n_basis
basis <- create.bspline.basis(rangeval = c(min_time, max_time), nbasis = n_basis, norder = n_order, breaks = knots)
n_basis
is.fd(basis)

#Evaluating the basis functions at the times where our data curve was observed i.e on 1187 points 
PHI = eval.basis(time_points, basis) 
dim(PHI)
PHI #1187x24 Matrix 

par(mfrow = c(1, 1))
#Plotting the basis functions and placement of the knots 
matplot(time_points,PHI,type='l',lwd=1,lty=1, xlab='days',ylab='B-spline function values',cex.lab=1,cex.axis=1)
for (i in 1:n_knots)
{
  abline(v=knots[i], lty=2, lwd=1)
}

#Finding the optimal number of basis functions for the fourier basis
K_values <- seq(4, 100, by = 2)  

# We'll create an empty data frame to store the results
results_df <- data.frame(
  Asset = character(),
  K     = numeric(),
  s2    = numeric(),
  stringsAsFactors = FALSE
)

n <- length(time_points)
assets <- colnames(Prices_data)

for (asset_name in assets) {
  
  # 2a) Extract the asset's time series as a numeric vector Y
  #     Must match length = n
  Y <- as.numeric(Prices_data[[asset_name]])
  
  # 2b) Loop over each K in K_values
  for (k in K_values) {
    
    fourier_basis <- create.fourier.basis(
      rangeval = c(min(time_points), max(time_points)),
      nbasis   = k
    )
    
    Phi <- eval.basis(time_points, fourier_basis)
    c_hat <- solve(t(Phi) %*% Phi, t(Phi) %*% Y)
    
    Y_hat <- Phi %*% c_hat
    RSS <- sum((Y - Y_hat)^2)
    s2_val <- RSS / (n - k)
    
    results_df <- rbind(
      results_df,
      data.frame(
        Asset = asset_name,
        K     = k,
        s2    = s2_val,
        stringsAsFactors = FALSE
      )
    )
  }
}

print(results_df)
best_k_df <- results_df %>%
  group_by(Asset) %>%
  slice_min(s2)

print(best_k_df)
ggplot(
  data = subset(results_df, Asset == assets[1]),
  aes(x = K, y = s2)
) +
  geom_line() +
  geom_point() +
  labs(title = paste("Residual variance vs K for asset:", assets[14]))

#By visualizing the elbow point in the plots, we deduce that 22 basis functions is adequate 

#Creating the Fourier Basis:
n_basis_fourier <- 23
basis_fourier <- create.fourier.basis(
  rangeval = c(min_time, max_time),
  nbasis   = n_basis_fourier
)

# Evaluate Fourier basis functions at the time points
PHI_fourier <- eval.basis(time_points, basis_fourier)
dim(PHI_fourier)

# Prepare a data frame to store MSE for each asset
mse_df <- data.frame(
  Asset         = colnames(Prices_data),
  MSE_Bspline   = NA_real_,
  MSE_Fourier   = NA_real_,
  stringsAsFactors = FALSE
)

# Loop over columns (assets)
for (i in seq_along(colnames(Prices_data))) {
  
  # Extract the i-th asset's time series as a numeric vector
  Y <- as.numeric(Prices_data[[i]])
  
  # --- (A) Fit using B-spline basis ---
  # Solve for coefficients: b_bspline = (PHI'PHI)^(-1) PHI'Y
  coeff_bspline <- solve(t(PHI) %*% PHI, t(PHI) %*% Y)
  
  # Predicted values
  Yhat_bspline <- PHI %*% coeff_bspline
  
  # MSE (B-spline)
  mse_bspline <- mean((Y - Yhat_bspline)^2)
  
  # --- (B) Fit using Fourier basis ---
  coeff_fourier <- solve(t(PHI_fourier) %*% PHI_fourier, t(PHI_fourier) %*% Y)
  Yhat_fourier <- PHI_fourier %*% coeff_fourier
  
  mse_fourier <- mean((Y - Yhat_fourier)^2)
  
  # --- (C) Store MSE values in the data frame
  mse_df$MSE_Bspline[i] <- mse_bspline
  mse_df$MSE_Fourier[i] <- mse_fourier
}

print(mse_df)
# Suppose 'mse_df' has columns: Asset, MSE_Bspline, MSE_Fourier
mse_df$Winner <- ifelse(mse_df$MSE_Fourier < mse_df$MSE_Bspline, 
                        "Fourier", 
                        ifelse(mse_df$MSE_Fourier > mse_df$MSE_Bspline,
                               "B-spline","Tie"))
table(mse_df$Winner)

#Weight Vector due to correlated error terms - e_t=rho*e_t-1+epsilon
#Check whether the residuals are autocorrelated by Durbin Watson or other tests/plots
#Need to estimate rho and simgm^2 where rho is the autocorrelation paramter and sigma^2 is the variance of the noise terms epsilon assuming an AR(1) model for each asset
#Then we can find the estimated covariance matrix of the residuals for each asset 
#Then we can find its inverse and use it as the weight matrix wtmat
#A specified value of lambda is used initially to smooth the data and obtain the residuals to find the wtvec for each asset
#Then this new wtvec would be used as input argument to the smooth.basis function in order to find the optimal lambda value for each asset 

# Data preparation
Prices_data_numeric <- as.matrix(Prices_data) 
View(Prices_data_numeric)
n_assets <- ncol(Prices_data_numeric)          
n_assets
n_time <- nrow(Prices_data_numeric)            
n_time
#time_points <- seq(0, 1, length.out = n_time) 
#time_points
Sec_derivative <- int2Lfd(max(0, n_order - 2))

residuals_matrix <- matrix(NA, nrow = n_time, ncol = n_assets)  # Store residuals for each asset
fitted_values_matrix <- matrix(NA, nrow = n_time, ncol = n_assets)  # Store fitted values for each asset

# Loop over each asset (column)
for (asset in 1:n_assets) {
  cat(paste("Performing OLS for asset", asset, "\n"))
  
  # Extract data for the current asset
  y <- Prices_data_numeric[, asset]
  
  # Construct the basis function matrix without regularization
  Phi <- eval.basis(time_points, basis)  # Evaluate basis at time points
  
  # Perform OLS estimation using the normal equation: c_hat = (Phi'Phi)^(-1) Phi'Y
  c_hat <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% y
  
  # Compute the fitted values: Y_hat = Phi * c_hat
  fitted_values <- Phi %*% c_hat
  
  # Calculate residuals: Observed - Fitted
  residuals <- y - fitted_values
  
  # Store results in matrices
  residuals_matrix[, asset] <- residuals  # Store residuals
  fitted_values_matrix[, asset] <- fitted_values  # Store fitted values
}

# View residuals and fitted values for inspection
View(residuals_matrix)  # Residuals matrix: rows = time points, columns = assets
View(fitted_values_matrix)
View(Prices_data_numeric)


#Checking whether there actually is auto-correlation between the residuals through numerous checks

#PACF Checks
# Load necessary library
library(forecast)

# Ensure residuals_matrix and Prices_data_numeric are properly defined
n_assets <- ncol(Prices_data_numeric)

# Set the seed for reproducibility
set.seed(22)  # For consistent random sampling

# Randomly sample 6 assets
sampled_assets <- sample(1:n_assets, 6)

# Set up plotting window: one PACF plot per asset
par(mfrow = c(3, 2))  # 3 rows, 2 columns for 6 plots

# Loop through each sampled asset and plot the PACF
for (asset in sampled_assets) {
  asset_name <- colnames(Prices_data_numeric)[asset]  # Get the asset name
  
  # Calculate the PACF for the current asset's residuals
  pacf_result <- pacf(residuals_matrix[, asset], plot = FALSE, lag.max = 50)
  
  # Extract the lag with the highest PACF value (excluding lag 0)
  lag_values <- pacf_result$lag
  pacf_values <- pacf_result$acf
  max_pacf_index <- which.max(abs(pacf_values))  # Find the index of the highest PACF value
  max_lag <- lag_values[max_pacf_index]
  max_pacf <- pacf_values[max_pacf_index]
  
  # Display the lag and PACF value
  cat(paste("Asset", asset_name, "Highest PACF Value at Lag", max_lag, ":", round(max_pacf, 3), "\n"))
  
  # Plot the PACF for the current asset
  main_title <- paste("PACF of Residuals for", asset_name, 
                      "\nMax PACF at Lag", max_lag, "=", round(max_pacf, 3))
  pacf(residuals_matrix[, asset], main = main_title, lag.max = 50)
}



#The plots suggest fitting an AR(1) model as suspected

#Conducting the Ljung-Box Test to check the autocorrelation in the residuals of all assets
#Checking the null hypothesis of no autocorrelation up to a certain lag


# Number of lags to check for autocorrelation
lag_k <- 1

# Loop through each asset
# Loop through each asset
for (asset in 1:n_assets) {
  # Get the name of the current asset
  asset_name <- colnames(Prices_data_numeric)[asset]
  
  # Perform the Ljung-Box test on the residuals of the current asset
  lb_test <- Box.test(residuals_matrix[, asset], lag = lag_k, type = "Ljung-Box")
  
  # Print the result for the current asset with its name
  cat(paste("Ljung-Box Test for Asset:", asset_name, "\n"))
  print(lb_test)
}

#P values are smaller than 0.05 for all assets so we accept H1: The auto-correlations up to lag 1 are not zero


#N.B.
#The PACF cuts off after a small number of lags (e.g., 1 or 2) because once you account for the correlation with the immediate past value, there's little to no direct correlation at further lags.
#However, the Ljung-Box test still detects significant autocorrelation across several lags (even if individual lags are small) because it looks at the cumulative autocorrelation over a range of lags.

# Reset plotting layout
par(mfrow = c(1, 1))

#Plotting a scatter plot of residuals at t+1 against residual at t to check for any correlation 
n_assets <- ncol(residuals_matrix)
n_time <- nrow(residuals_matrix)

# Check autocorrelation in residuals for each asset
for (asset in 1:n_assets) {
  asset_name <- colnames(Prices_data_numeric)[asset]  # Get asset name
  
  cat(paste("Checking autocorrelation for Asset:", asset_name, "\n"))
  
  # Scatter plot of residuals vs. lagged residuals (to visually inspect autocorrelation)
  plot(residuals_matrix[1:(n_time - 1), asset], residuals_matrix[2:n_time, asset], 
       xlab = "Residual(t)", ylab = "Residual(t+1)", 
       main = paste("Scatter plot of residuals for", asset_name))
}


#So fit an AR(1) Model to the residuals 

#Fitting the AR(1) model and estimating the variance covariance matrix of the residuals

# Load the necessary library for AR(1) model fitting
library(forecast)

# Assuming residuals_matrix is already prepared with residuals for each asset (each column)
n_assets <- ncol(residuals_matrix)
n_time <- nrow(residuals_matrix)

# Initialize storage for AR(1) parameters and the covariance matrix
rho_values <- numeric(n_assets)    # Autocorrelation parameter for each asset
sigma2_values <- numeric(n_assets) # Variance of noise for each asset

# Initialize an empty list to store the covariance matrices for each asset
cov_matrix <- vector("list", n_assets)  # List to store the covariance matrix for each asset

# Fit AR(1) model to each asset's residuals and estimate parameters
wtvec_list <- list()  # List to store the inverse covariance matrices for each asset

for (asset in 1:n_assets) {
  tryCatch({
    # Fit AR(1) model to the residuals of the current asset
    ar1_model <- arima(residuals_matrix[, asset], order = c(1, 0, 0))  # AR(1) model
    
    # Extract parameters from the AR(1) model
    rho_values[asset] <- ar1_model$coef[1]  # Autocorrelation parameter (rho)
    sigma2_values[asset] <- ar1_model$sigma2  # Variance of noise (sigma^2)
    
    # Print the estimated parameters for each asset
    cat(paste("Asset", asset, "AR(1) Model Parameters:\n"))
    cat("Autocorrelation parameter (rho):", rho_values[asset], "\n")
    cat("Variance of noise (sigma^2):", sigma2_values[asset], "\n")
    
    # Initialize the covariance matrix for this asset
    asset_cov_matrix <- matrix(NA, nrow = n_time, ncol = n_time)  # Covariance matrix for residuals
    
    # Fill in the covariance matrix using the AR(1) covariance structure
    for (t1 in 1:n_time) {
      for (t2 in 1:n_time) {
        # Compute the lag (h) between the two time points
        h <- abs(t1 - t2)
        
        # Calculate the covariance between the two time points
        asset_cov_matrix[t1, t2] <- sigma2_values[asset] * (rho_values[asset]^h)
      }
    }
    
    # Store the covariance matrix for the current asset
    cov_matrix[[asset]] <- asset_cov_matrix
    
    # Store the inverse of the covariance matrix as the weight vector
    wtvec_list[[asset]] <- solve(asset_cov_matrix)
    
    # Optionally, you can print the covariance matrix for inspection
    cat(paste("Covariance Matrix for Asset", asset, ":\n"))
    print(asset_cov_matrix)
  }, error = function(e) {
    # If there's an error, print a message and move on to the next asset
    cat(paste("Error fitting AR(1) model for Asset", asset, ": ", e$message, "\n"))
  })
}

# Optionally, store the rho and sigma2 values for further analysis
parameter_estimates <- data.frame(
  Asset = 1:n_assets,
  Rho = rho_values,
  Sigma2 = sigma2_values
)

print(parameter_estimates)
#Checking
dim(residuals_matrix)  # Ensure this is (n_time, n_assets)
# Checking the dimensions of the covariance matrix
dim(cov_matrix[[106]])# Should return (n_time, n_time)
print(cov_matrix[[1]])
cov_matrix[[1]][2, 1]  # First element of the covariance matrix
wtvec_list[[6]][1,1]
#There is a problem with estimating rho and sigma^2 for asset 25 and asset 104 - Ask (Constant values of residuals?) 

# Define the range of lambda values
lambda_range <- seq(1e-8, 1e4, length.out = 50)

# Data frame to store optimal lambda values based on GCV
optimal_lambdas_gcv <- data.frame(
  Asset = colnames(Prices_data_numeric),
  Optimal_lambda_GCV = NA,
  Min_GCV = NA,
  stringsAsFactors = FALSE
)

# Loop through each asset
for (asset in colnames(Prices_data_numeric)) {
  cat(paste("Calculating GCV for asset:", asset, "\n"))
  
  # Extract the data for the current asset
  y <- Prices_data_numeric[, asset]
  wtvec <- wtvec_list[[asset]]
  
  # Initialize storage for GCV values
  gcvsave <- numeric(length(lambda_range))
  
  # Loop through lambda values
  for (i in seq_along(lambda_range)) {
    lambda <- lambda_range[i]  
    fdParobj <- fdPar(basis, Sec_derivative, lambda = lambda)  
    
    # Perform smoothing
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj,wtvec = wtvec)
    
    # Extract GCV value
    gcvsave[i] <- smooth_result$gcv  
  }
  
  # Identify the optimal lambda corresponding to the minimum GCV value
  min_gcv_index <- which.min(gcvsave)
  optimal_lambda_gcv <- lambda_range[min_gcv_index]
  min_gcv <- gcvsave[min_gcv_index]
  
  # Store the optimal results for this asset
  optimal_lambdas_gcv[optimal_lambdas_gcv$Asset == asset, ] <- list(
    Asset = asset,
    Optimal_lambda_GCV = optimal_lambda_gcv,
    Min_GCV = min_gcv
  )
  
  # Plot GCV curve for the current asset
  plot(lambda_range, gcvsave, type = "b", col = "blue", lwd = 2, pch = 19,
       xlab = "Lambda", ylab = "GCV", 
       main = paste("GCV Curve for Asset", asset))
  
  # Highlight the optimal lambda on the plot
  abline(v = optimal_lambda_gcv, col = "red", lty = 2, lwd = 1.5)
}

# Print optimal lambda values based on GCV
print(optimal_lambdas_gcv)


# Smoothing all assets with optimal lambdas
smoothed_assets <- list()

# Loop through each asset to apply smoothing
for (asset in optimal_lambdas_gcv$Asset) {
  cat(paste("Smoothing asset:", asset, "\n"))
  
  tryCatch({
    # Extract the optimal lambda for the current asset
    optimal_lambda <- optimal_lambdas_gcv[optimal_lambdas_gcv$Asset == asset, "Optimal_lambda_GCV"]
    
    # Extract the data and weight vector for the current asset
    y <- Prices_data_numeric[, asset]
    wtvec <- wtvec_list[[asset]]
    
    # Define the fdPar object using the optimal lambda
    fdParobj <- fdPar(basis, Sec_derivative, lambda = optimal_lambda)
    
    # Perform smoothing using the optimal lambda
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj, wtvec = wtvec)
    
    # Store the smoothed functional data object in the list
    smoothed_assets[[asset]] <- smooth_result$fd
  }, error = function(e) {
    cat(paste("Error smoothing asset:", asset, "- Skipping.\n"))
  })
}

# Combine all smoothed functional data into a single fd object
coef_matrix <- sapply(smoothed_assets, function(fd_obj) fd_obj$coefs)
combined_fd <- fd(coef = coef_matrix, basisobj = basis)

# Plot all smoothed functions using plot.fd
plot(
  combined_fd,
  xlab = "Date", ylab = "Smoothed Value",
  main = "Smoothed Functions for All Assets"
)
#Conducting FPCA
fpca_results_AR <- pca.fd(fdobj = combined_fd, nharm = 2)
fpca_results_AR$harmonics #a functional data object for the harmonics or eigen-functions
fpca_results_AR$scores #s matrix of scores on the principal components or harmonics
fpca_results_AR$values #the complete set of eigenvalues
fpca_results_AR$varprop #a vector giving the proportion of variance explained by each eigen-function
fpca_results_AR$meanfd #a functional data object giving the mean function
plot.pca.fd(fpca_results_AR)


#Plotting the harmonics individually (functional principal components)
#Plotting the first harmonic
plot(fpca_results_AR$harmonics[1], main = "First Harmonic (PC1)", xlab = "Date", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180), 
     labels = as.Date(seq(min(Prices$Date), max(Prices$Date), by = "180 days")), 
     las = 2)

#Plotting the second harmonic
plot(fpca_results_AR$harmonics[2], main = "Second Harmonic (PC2)", xlab = "Date", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180), 
     labels = as.Date(seq(min(Prices$Date), max(Prices$Date), by = "180 days")), 
     las = 2)

# Convert FPCA scores to a data frame
fpc_scores <- as.data.frame(fpca_results_AR$scores)

# Rename the columns (assuming you extracted 2 principal components/harmonics)
colnames(fpc_scores) <- c("PC1", "PC2")

# Add a column for the asset names. Here, we assume the order of scores 
# corresponds to the columns in `Prices_data_numeric`.
fpc_scores$Asset <- colnames(Prices_data_numeric)

# Reorder the columns to have the Asset name first
fpc_scores <- fpc_scores[, c("Asset", "PC1", "PC2")]

# Print the FPCA scores
print(fpc_scores)

# Optionally, view the scores in a neat table (if using RStudio)
View(fpc_scores)

# Plot the FPCA scores in a scatter plot
library(ggplot2)
ggplot(fpc_scores, aes(x = PC1, y = PC2, label = Asset)) +
  geom_point(color = "blue", size = 3) +             # Plot the points
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +    # Add labels near points
  labs(title = "FPCA Scores Scatter Plot",
       x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()

# Compute the first derivatives of the smoothed functions
combined_fd_deriv <- deriv.fd(combined_fd, Lfdobj = 1)

# Define the assets of interest
desired_assets <- c("ITM.L", "SPMR.MI")

# Identify the indices for these assets in the Prices_data_numeric matrix
selected_indices <- which(colnames(Prices_data_numeric) %in% desired_assets)

# Set up a multi-plot layout: two plots side by side
par(mfrow = c(1, length(selected_indices)))

# Loop through each selected asset and plot its derivative
for (i in selected_indices) {
  asset_name <- colnames(Prices_data_numeric)[i]
  
  # Plot the derivative for the current asset
  plot(combined_fd_deriv[i], 
       xlab = "Time (Days)", 
       ylab = "First Derivative", 
       main = paste("Derivative of", asset_name))
}

# Reset the plotting layout
par(mfrow = c(1, 1))
View(Prices_data_numeric)
View(Prices_data)
View(Prices)
