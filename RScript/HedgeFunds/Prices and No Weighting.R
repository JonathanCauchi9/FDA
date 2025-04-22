install.packages("refund")
install.packages("tidyr")
install.packages("tseries")
install.packages("urca")
install.packages("FinTS")
install.packages("rugarch")
install.packages("tseries")
install.packages("seastests")
library(rugarch)
library(ggplot2)
library(FinTS)
library(dplyr)
library(fda)
library(refund)
library(tidyr)
library(ggplot2)
library(forecast)
library(tseries)
library(seastests)
View(pricesdata)
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

stl_results <- list()  # Store the STL results for later inspection

par(mfrow = c(2, 2))  # Adjust rows/cols to your preference

for (asset in colnames(Prices_data)) {
  # Convert the data for 'asset' into a time series.
  # Here, we assume daily data with possible yearly seasonality => frequency=365
  # Adjust 'frequency' to reflect your data’s actual expected periodic cycle.
  ts_data <- ts(Prices_data[[asset]], frequency = 365)
  
  # Perform STL decomposition
  stl_decomp <- stl(ts_data, s.window = "periodic")
  
  # Store results in a list if you want to inspect them later
  stl_results[[asset]] <- stl_decomp
  
  # Plot the decomposition
  plot_title <- paste("STL Decomposition -", asset)
  plot(stl_decomp, main = plot_title)
}
par(mfrow = c(1, 1))

# Initialize a data frame to store the seasonality test results for asset prices
seasonality_results <- data.frame(
  Asset = character(),
  IsSeasonal = logical(),
  stringsAsFactors = FALSE
)

# Loop through each asset in Prices_data to test for seasonality
for (asset in colnames(Prices_data)) {
  # Convert the asset price series to a time series object with daily frequency
  # (adjust the 'frequency' argument if needed, e.g., 252 for trading days)
  ts_data <- ts(Prices_data[[asset]], frequency = 365)
  
  # Check for seasonality using isSeasonal
  seasonal_check <- isSeasonal(ts_data)
  
  # Provide immediate feedback in the console
  if (seasonal_check) {
    cat("Seasonality detected for asset:", asset, "\n")
  } else {
    cat("No significant seasonality detected for asset:", asset, "\n")
  }
  
  # Store the result in the seasonality_results data frame
  seasonality_results <- rbind(
    seasonality_results,
    data.frame(Asset = asset, IsSeasonal = seasonal_check, stringsAsFactors = FALSE)
  )
}



#Checking for stationarity
# Load necessary packages
library(tseries)

# Initialize vectors to store results
# Load necessary packages
library(tseries)
library(seastests)

# Initialize vectors to store initial stationarity results
assets <- colnames(Prices_data)
p_values <- numeric(length(assets))
is_stationary <- logical(length(assets))

# Create a new dataframe to store the adjusted (differenced) data
Prices_data_adjusted <- Prices_data

# ---- Step 1: Initial ADF Test and Differencing for Non-Stationary Assets ----
for (i in seq_along(assets)) {
  # Convert data to time series
  ts_data <- ts(Prices_data[[assets[i]]], frequency = 365)
  
  # Perform ADF test
  adf_test <- adf.test(ts_data, alternative = "stationary", k = 0)
  
  # Store p-value and stationarity status
  p_values[i] <- adf_test$p.value
  is_stationary[i] <- p_values[i] < 0.05
  
  # Apply first-order differencing if the asset is non-stationary
  if (!is_stationary[i]) {
    differenced_data <- diff(ts_data, differences = 1)
    
    # Pad the last row to maintain original length (copy last value)
    differenced_data <- c(differenced_data, tail(differenced_data, 1))
    
    # Replace the original column with the differenced data
    Prices_data_adjusted[[assets[i]]] <- differenced_data
  }
}

# Combine initial stationarity results
stationarity_results <- data.frame(Asset = assets, P_Value = p_values, IsStationary = is_stationary)

# Print initial stationarity results
print(stationarity_results)


# ---- Step 2: ADF Test on Adjusted Data ----
p_values_adjusted <- numeric(length(assets))
is_stationary_adjusted <- logical(length(assets))

for (i in seq_along(assets)) {
  ts_data_adjusted <- ts(Prices_data_adjusted[[assets[i]]], frequency = 365)
  adf_test_adjusted <- adf.test(ts_data_adjusted, alternative = "stationary", k = 0)
  p_values_adjusted[i] <- adf_test_adjusted$p.value
  is_stationary_adjusted[i] <- p_values_adjusted[i] < 0.05
}

# Combine results for the adjusted data
stationarity_results_adjusted <- data.frame(Asset = assets, 
                                            P_Value = p_values_adjusted, 
                                            IsStationary = is_stationary_adjusted)

# Print the stationarity results after differencing
print(stationarity_results_adjusted)


# ---- Step 3: Seasonality Test on Stationary Data ----
is_seasonal <- logical(length(assets))

for (i in seq_along(assets)) {
  ts_data_seasonal <- ts(Prices_data_adjusted[[assets[i]]], frequency = 365)
  is_seasonal[i] <- isSeasonal(ts_data_seasonal)
}

# Combine results into a data frame
seasonality_results <- data.frame(Asset = assets, IsSeasonal = is_seasonal)

# Print seasonality results
print(seasonality_results)

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
PHI #1187x25 Matrix 

par(mfrow = c(1, 1))
#Plotting the basis functions and placement of the knots 
matplot(time_points,PHI,type='l',lwd=1,lty=1, xlab='days',ylab='B-spline function values',cex.lab=1,cex.axis=1)
for (i in 1:n_knots)
{
    abline(v=knots[i], lty=2, lwd=1)}

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




Prices_data_numeric <- as.matrix(Prices_data) 

# Define lambda range
lambda_range <- seq(1e-8, 1e4, length.out = 50)  # Lambda values to test

# Initialize storage for results
optimal_lambdas <- data.frame(
  Asset = colnames(Prices_data_numeric),
  Optimal_lambda = NA,
  Min_GCV = NA
)
n_assets <- ncol(Prices_data_numeric)
# Loop through each asset
for (asset in 1:n_assets) {
  asset_name <- colnames(Prices_data_numeric)[asset]
  cat(paste("Processing Asset:", asset_name, "\n"))
  
  y <- Prices_data_numeric[, asset]
  
  # Initialize GCV storage
  gcvsave <- numeric(length(lambda_range))
  
  # Loop through lambda values
  for (i in seq_along(lambda_range)) {
    lambda <- lambda_range[i]
    fdParobj <- fdPar(basis, Sec_derivative, lambda)
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj)
    gcvsave[i] <- smooth_result$gcv  # GCV for this lambda
  }
  
  # Find the optimal lambda with the minimum GCV
  min_gcv_index <- which.min(gcvsave)
  optimal_lambda <- lambda_range[min_gcv_index]
  min_gcv <- gcvsave[min_gcv_index]
  
  # Store the results
  optimal_lambdas[asset, ] <- list(
    Asset = asset_name,
    Optimal_lambda = optimal_lambda,
    Min_GCV = min_gcv
  )
  
  # Plot the GCV curve for the current asset
  plot(
    lambda_range, gcvsave, type = "b", col = "blue", lwd = 2, pch = 19,
    xlab = "Lambda", ylab = "GCV",
    main = paste("GCV Curve for Asset", asset_name)
  )
  abline(v = optimal_lambda, col = "red", lty = 2, lwd = 1.5)  # Highlight optimal lambda
}

print(optimal_lambdas)
lowest_gcv_row <- optimal_lambdas[which.min(optimal_lambdas$Min_GCV), ]

# Extract individual elements if desired
lowest_gcv_asset  <- lowest_gcv_row$Asset
lowest_gcv_asset
highest_gcv_row <- optimal_lambdas[which.max(optimal_lambdas$Min_GCV), ]

# Extract individual elements if desired
highest_gcv_asset  <- highest_gcv_row$Asset
highest_gcv_asset
# Extract optimal lambda and minimum GCV values from the optimal_lambdas data frame
optimal_lambda_values <- as.numeric(optimal_lambdas$Optimal_lambda)
min_gcv_values <- as.numeric(optimal_lambdas$Min_GCV)

# Print descriptive statistics for optimal lambdas
cat("Descriptive statistics for Optimal Lambdas (Smoothing Analysis):\n")
print(summary(optimal_lambda_values))
cat("Standard Deviation:", sd(optimal_lambda_values), "\n\n")

# Print descriptive statistics for minimum GCV values
cat("Descriptive statistics for Minimum GCVs (Smoothing Analysis):\n")
print(summary(min_gcv_values))
cat("Standard Deviation:", sd(min_gcv_values), "\n")

# Compute and print the total number of time points
n_time <- length(time_points)  # or alternatively: nrow(Prices_data)
cat("Total number of time points:", n_time, "\n")


# Smooth all assets with their respective optimal lambdas
smoothed_assets <- list()
for (asset in 1:n_assets) {
  asset_name <- colnames(Prices_data_numeric)[asset]
  cat(paste("Smoothing Asset:", asset_name, "\n"))
  
  optimal_lambda <- optimal_lambdas$Optimal_lambda[asset]
  y <- Prices_data_numeric[, asset]
  fdParobj <- fdPar(basis, Sec_derivative, lambda = optimal_lambda)
  smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj)
  smoothed_assets[[asset_name]] <- smooth_result$fd
}

# Combine smoothed functions into a single functional data object for plotting
combined_fd <- fd(
  coef = do.call(cbind, lapply(smoothed_assets, function(fd) fd$coefs)),
  basisobj = basis
)

# Plot all smoothed functions
plot(
  combined_fd,
  xlab = "Time", ylab = "Smoothed Value",
  main = "Smoothed Functions for all Assets"
)

#Performing Functional Principal Component Analysis (FPCA) 
fpca_results <- pca.fd(fdobj = prices.fd, nharm = 2)
fpca_results$harmonics #a functional data object for the harmonics or eigen-functions
fpca_results$scores #s matrix of scores on the principal components or harmonics
fpca_results$values #the complete set of eigenvalues
fpca_results$varprop #a vector giving the proportion of variance explained by each eigen-function
fpca_results$meanfd #a functional data object giving the mean function
plot.pca.fd(fpca_results)


#Conducting VARIMAX Rotation
fpca_resultsrot <- varmx.pca.fd(fpca_results, nharm = 2)
fpca_resultsrot$harmonics #a functional data object for the harmonics or eigen-functions
fpca_resultsrot$scores #s matrix of scores on the principal components or harmonics
fpca_resultsrot$values #the complete set of eigenvalues
fpca_resultsrot$varprop #a vector giving the proportion of variance explained by each eigen-function
fpca_resultsrot$meanfd #a functional data object giving the mean function
plot.pca.fd(fpca_resultsrot)

#Plotting the harmonics individually (functional principal components)
#Plotting the first harmonic
plot(fpca_results$harmonics[1], main = "First Harmonic (PC1)", xlab = "Date", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180), 
     labels = as.Date(seq(min(Prices$Date), max(Prices$Date), by = "180 days")), 
     las = 2)

#Plotting the second harmonic
plot(fpca_results$harmonics[2], main = "Second Harmonic (PC2)", xlab = "Date", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180), 
     labels = as.Date(seq(min(Prices$Date), max(Prices$Date), by = "180 days")), 
     las = 2)




#Checking where there was the largest change in price
#Calculate the absolute difference in price for the JMAT.L Asset
abs(diff(Prices_data$JMAT.L))
Prices_data_change <- abs(diff(Prices_data$JMAT.L))
View(Prices_data_change)
#Find the index of the largest change
max_change_index <- which.max(Prices_data_change)
max_change_index
#Get the corresponding date and price change
date_of_max_change <- pricesdata$Date[max_change_index + 1] # +1 because `diff` reduces length by 1
max_change_value <- Prices_data_change[max_change_index]
max_change_value
#Print the results
cat("Date of largest change:", date_of_max_change, "\n")
cat("Largest price change:", max_change_value, "\n")

#Derivatives

# Select three assets (modify these names/indexes as needed)
subset_assets <- names(smoothed_assets)[1:3]
subset_fd <- smoothed_assets[subset_assets]

# Compute the first derivative for each selected functional data object
# Lfdobj = 1 specifies that we want the first derivative
derivative_fd_list <- lapply(subset_fd, function(fd_obj) {
  deriv.fd(fd_obj, Lfdobj = 1)
})

# Evaluate the derivative functions at the original time points
deriv_values_list <- lapply(derivative_fd_list, function(fd_obj) {
  eval.fd(time_points, fd_obj)
})

# Option 1: Plot each derivative in a separate panel
par(mfrow = c(1, length(subset_assets)))  # Arrange plots side-by-side
for (i in seq_along(deriv_values_list)) {
  plot(time_points, deriv_values_list[[i]],
       type = "l", col = "blue", lwd = 2,
       main = paste("1st Derivative of", subset_assets[i]),
       xlab = "Time", ylab = "1st Derivative")
}
par(mfrow = c(1, 1))  # Reset to default

# Option 2 (alternative): Overlay all derivatives on one plot
plot(time_points, deriv_values_list[[1]],
     type = "l", col = "blue", lwd = 2,
     main = "First Derivatives of Selected Assets",
     xlab = "Time", ylab = "1st Derivative")
cols <- c("red", "darkgreen")
for (i in 2:length(deriv_values_list)) {
  lines(time_points, deriv_values_list[[i]], col = cols[i-1], lwd = 2)
}
legend("topright", legend = subset_assets,
       col = c("blue", cols), lty = 1, bty = "n")

