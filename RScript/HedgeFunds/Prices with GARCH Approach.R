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
PHI #1187x25 Matrix 

par(mfrow = c(1, 1))
#Plotting the basis functions and placement of the knots 
matplot(time_points,PHI,type='l',lwd=1,lty=1, xlab='days',ylab='B-spline function values',cex.lab=1,cex.axis=1)
for (i in 1:n_knots)
{
  abline(v=knots[i], lty=2, lwd=1)
}

# GARCH Method
# Checking for ADF Test Results for All Assets
library(tseries)  # For ADF test

# Load necessary library
library(tseries)

# Initialize results dataframe for ADF test
adf_results <- data.frame(
  Asset = colnames(Prices_data),
  ADF_p_value = NA,
  ADF_Test_Statistic = NA
)

# Initialize a list to store differenced data
differenced_data <- list()

# Loop through each asset to perform ADF test and store differenced data
for (asset in colnames(Prices_data)) {
  cat(paste("Performing ADF test for asset:", asset, "\n"))
  
  # Extract the data for the current asset
  asset_data <- Prices_data[[asset]]
  
  # Perform Augmented Dickey-Fuller (ADF) Test
  adf_test <- adf.test(asset_data, alternative = "stationary")
  
  # Store the results
  adf_results[adf_results$Asset == asset, ] <- data.frame(
    Asset = asset,
    ADF_p_value = adf_test$p.value,
    ADF_Test_Statistic = adf_test$statistic
  )
  
  # Perform first-order differencing and store the result in the list
  differenced_data[[asset]] <- diff(asset_data)
}

# Convert the list of differenced data into a dataframe
differenced_data <- do.call(cbind, differenced_data)

# Convert to a dataframe with proper column names
differenced_data <- as.data.frame(differenced_data)
View(differenced_data)
# Display ADF results
print(adf_results)

# Display a preview of the differenced data
head(differenced_data)


#Ask about these assets which are stationary !
stationary_assets <- adf_results %>% filter(ADF_p_value < 0.05)
stationary_assets
adf_results
View(differenced_data)
# Perform ADF Test on Differenced Data

# Initialize results dataframe for differenced data
adf_diff_results <- data.frame(
  Asset = colnames(differenced_data), 
  ADF_p_value = NA,
  ADF_Test_Statistic = NA
)

# Loop through each asset in the differenced data
for (asset in colnames(differenced_data)) {  # Exclude 'Date' column
  cat(paste("Performing ADF test on differenced data for asset:", asset, "\n"))
  
  # Extract the differenced data for the current asset
  diff_data <- differenced_data[[asset]]
  
  # Perform Augmented Dickey-Fuller (ADF) Test
  adf_test <- adf.test(diff_data, alternative = "stationary")
  
  # Store the results
  adf_diff_results[adf_diff_results$Asset == asset, ] <- data.frame(
    Asset = asset,
    ADF_p_value = adf_test$p.value,
    ADF_Test_Statistic = adf_test$statistic
  )
}
adf_diff_results
#All p-values smaller than 0.05 so the data is no longer stationary so now we fit ARMA on the differenced data. 

#Fitting the ARMA model on each of the differenced assets
# Ask on auto.arima function 
library(forecast)  # For ARMA modeling

# Initialize lists to store results
arma_models <- list()
arma_residuals <- list()

# Fit ARMA model for each differenced column (asset)
# Ask to make sure that you need to fit the ARMA model on the differenced data 
for (asset in colnames(differenced_data)) {  
  cat(paste("Fitting ARMA model for asset:", asset, "\n"))
  
  # Extract the differenced data for the current asset
  diff_data <- differenced_data[[asset]]
  
  # Fit ARMA model using auto.arima() to select p and q (d = 0 because data is already differenced)
  arma_model <- auto.arima(diff_data, d = 0, seasonal = FALSE)
  
  # Store the model and residuals
  arma_models[[asset]] <- arma_model
  arma_residuals[[asset]] <- residuals(arma_model)
  
  # Print model summary
  print(summary(arma_model))
}

# Combine ARMA model summaries into a single table for review
arma_model_summaries <- data.frame(
  Asset = names(arma_models),
  AR_Order = sapply(arma_models, function(model) model$arma[1]),
  MA_Order = sapply(arma_models, function(model) model$arma[2]),
  AIC = sapply(arma_models, AIC),
  BIC = sapply(arma_models, BIC)
)
arma_model_summaries

#Obtaining the residuals in the form of a matrix
n_assets <- length(arma_models)# Number of assets (columns in differenced_data)
n_assets
n_time <- nrow(differenced_data) # Number of time points (rows in differenced_data)
n_time

# Create an empty residuals matrix
residuals_matrix <- matrix(NA, nrow = n_time, ncol = n_assets,
                           dimnames = list(NULL, names(arma_models)))

# Loop through each fitted ARMA model to extract residuals
for (asset in names(arma_models)) {
  cat(paste("Extracting residuals for asset:", asset, "\n"))
  
  # Get residuals from the fitted ARMA model
  residuals <- arma_residuals[[asset]]
  
  # Adjust the length to match the matrix dimensions (if needed)
  residuals_matrix[1:length(residuals), asset] <- residuals
}

# Convert the residuals matrix to a dataframe for easier handling
residuals_df <- as.data.frame(residuals_matrix)
View(residuals_df)

#Testing for ARCH effects
library(tseries)

# Loop through each asset's residuals
# Initialize results dataframe
garch_effects_results <- data.frame(
  Asset = names(arma_models),
  Ljung_Box_p_value = NA,
  ARCH_LM_p_value = NA
)

# Take lag=ln(1187)
# Loop through each asset's residuals
for (asset in names(arma_models)) {
  cat(paste("Testing for GARCH effects on asset:", asset, "\n"))
  
  # Extract residuals for the current asset
  residuals <- arma_residuals[[asset]]
  
  # Square the residuals to test for volatility clustering
  squared_residuals <- residuals^2
  
  # Perform Ljung-Box test on squared residuals
  ljung_box_test <- Box.test(squared_residuals, lag = 8, type = "Ljung-Box")
  
  # Perform ARCH LM test
  arch_lm_test <- ArchTest(residuals, lags = 8)
  
  # Store the results
  garch_effects_results[garch_effects_results$Asset == asset, ] <- list(
    Asset = asset,
    Ljung_Box_p_value = ljung_box_test$p.value,
    ARCH_LM_p_value = arch_lm_test$p.value
  )
}

# Identify assets with non-significant GARCH effects (p-value > 0.05)
non_significant_garch <- garch_effects_results %>% 
  filter(Ljung_Box_p_value > 0.05 | ARCH_LM_p_value > 0.05)

significant_garch <- garch_effects_results %>%
  filter(Ljung_Box_p_value <= 0.05 & ARCH_LM_p_value <= 0.05)

print(non_significant_garch)

# Initialize a list to store the fitted GARCH models
garch_models <- list()

# Loop through each asset and fit GARCH model for those with significant effects
for (asset in significant_garch$Asset) {
  cat(paste("Fitting GARCH(1,1) model for asset:", asset, "\n"))
  
  # Extract residuals for the current asset
  residuals <- arma_residuals[[asset]]
  
  # Define the GARCH(1,1) specification
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"  # Assuming normal distribution for residuals
  )
  
  # Fit the GARCH model
  garch_fit <- ugarchfit(spec = garch_spec, data = residuals)
  
  # Store the fitted model
  garch_models[[asset]] <- garch_fit
  
  # Print summary of the fitted model
  print(garch_fit)
}

# Extracting the conditional variances for each asset
garch_variances <- list()

for (asset in names(garch_models)) {
  cat(paste("Extracting conditional variances for asset:", asset, "\n"))
  
  garch_fit <- garch_models[[asset]]
  
  if (!is.null(garch_fit)) {
    tryCatch({
      conditional_variances <- sigma(garch_fit)^2
      garch_variances[[asset]] <- conditional_variances
    }, error = function(e) {
      cat(paste("Error extracting variances for asset:", asset, "-", e$message, "\n"))
    })
  } else {
    cat(paste("No valid GARCH model for asset:", asset, ". Skipping...\n"))
  }
}

# Assign equal variance (1) for assets with non-significant GARCH effects
for (asset in non_significant_garch$Asset) {
  cat(paste("Assigning equal weights for asset with no significant GARCH effects:", asset, "\n"))
  garch_variances[[asset]] <- rep(1, length(arma_residuals[[asset]]))
}

# Convert the list to a dataframe
garch_variances_df <- as.data.frame(garch_variances)
View(garch_variances_df)

# Obtain the reciprocal of the conditional variances for weighting
reciprocal_variances <- list()

for (asset in names(garch_variances)) {
  cat(paste("Calculating reciprocal of conditional variances for asset:", asset, "\n"))
  
  conditional_variances <- garch_variances[[asset]]
  
  reciprocal_variances[[asset]] <- ifelse(conditional_variances > 0, 
                                          1 / conditional_variances, 
                                          NA)
}

# Convert reciprocal variances to dataframe
reciprocal_variances_df <- as.data.frame(reciprocal_variances)
View(reciprocal_variances_df)

new_row <- reciprocal_variances_df[nrow(reciprocal_variances_df), ]  # Extract the last row
reciprocal_variances_df <- rbind(reciprocal_variances_df, new_row) # Append the new row
View(reciprocal_variances_df)

Prices_data_numeric <- as.matrix(Prices_data)
View(Prices_data_numeric)

# List of assets to update
assets_to_update <- c("NZYM-B.CO", "ALIV-SDB.ST", "ASSA-B.ST", "SKF-B.ST", "VOLV-B.ST")

# Replace only the first hyphen with a dot for the specific assets in Prices_data_numeric
colnames(Prices_data_numeric) <- sapply(colnames(Prices_data_numeric), function(col) {
  if (col %in% assets_to_update) {
    sub("-", ".", col)  # Replace only the first occurrence of hyphen with dot
  } else {
    col  # Keep other names unchanged
  }
})

# Obtaining the GCV across lambda values for the assets using the GARCH approach
Sec_derivative <- int2Lfd(max(0, n_order - 2))
n_assets <- ncol(reciprocal_variances_df)  # Total number of assets
n_assets
# Define a larger range of lambda values
lambda_range <- seq(1e-8, 1e4, length.out = 50)  # Adjust range and steps as needed

preliminary_results <- list()  # To store GCV results for all assets
optimal_lambdas <- data.frame(Asset = colnames(reciprocal_variances_df), 
                              Optimal_lambda_GCV = NA, 
                              Min_GCV = NA)

# Loop through each asset
for (asset in colnames(reciprocal_variances_df)) {
  cat(paste("Calculating metrics for asset:", asset, "\n"))
  
  # Extract the data and weight vector for the current asset
  y <- Prices_data_numeric[, asset]
  wtvec <- reciprocal_variances_df[[asset]]
  
  # Initialize storage for GCV values
  gcvsave <- numeric(length(lambda_range))
  
  # Loop through lambda values
  for (i in seq_along(lambda_range)) {
    lambda <- lambda_range[i]  
    fdParobj <- fdPar(basis, Sec_derivative, lambda = lambda)  
    
    # Perform smoothing
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj, wtvec = wtvec)
    
    # Extract GCV value
    gcvsave[i] <- smooth_result$gcv  
  }
  
  # Identify the optimal lambda and corresponding GCV value
  min_gcv_index <- which.min(gcvsave)
  optimal_lambda_gcv <- lambda_range[min_gcv_index]
  min_gcv <- gcvsave[min_gcv_index]
  
  # Store the optimal results for this asset
  optimal_lambdas[optimal_lambdas$Asset == asset, ] <- list(
    Asset = asset,
    Optimal_lambda_GCV = optimal_lambda_gcv,
    Min_GCV = min_gcv
  )
  
  # Store the full GCV curve for plotting or further analysis
  preliminary_results[[asset]] <- list(
    lambda_range = lambda_range,
    gcvsave = gcvsave
  )
  
  # Plot GCV curve for the current asset
  plot(lambda_range, gcvsave, type = "b", col = "blue", lwd = 2, pch = 19,
       xlab = "Lambda", ylab = "GCV", 
       main = paste("GCV Curve for Asset", asset))
  abline(v = optimal_lambda_gcv, col = "red", lty = 2, lwd = 1.5)  # Highlight optimal GCV lambda
  
}

print(optimal_lambdas)

# Smoothing all assets with optimal lambdas
smoothed_assets <- list()

# Loop through each asset to apply smoothing
for (asset in optimal_lambdas$Asset) {
  cat(paste("Smoothing asset:", asset, "\n"))
  
  # Extract the optimal lambda for the current asset
  optimal_lambda <- optimal_lambdas[optimal_lambdas$Asset == asset, "Optimal_lambda_GCV"]
  
  # Extract the data and weight vector for the current asset
  y <- Prices_data_numeric[, asset]
  wtvec <- reciprocal_variances_df[[asset]]
  
  # Define the fdPar object using the optimal lambda
  fdParobj <- fdPar(basis, Sec_derivative, lambda = optimal_lambda)
  
  # Perform smoothing using the optimal lambda
  smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj, wtvec = wtvec)
  
  # Store the smoothed functional data object in the list
  smoothed_assets[[asset]] <- smooth_result$fd
}
smoothed_assets
length(smoothed_assets)

# Combine smoothed functions into a single functional data object
combined_fd_garch <- fd(
  coef = do.call(cbind, lapply(smoothed_assets, function(fd) fd$coefs)), 
  basisobj = basis
)
combined_fd_garch$coefs[,1]
#Plotting the functions
plot(combined_fd_garch, xlab = "Time", ylab = "Value", main = "Smoothed Functions")
evaluated_values <- eval.fd(time_points, combined_fd_garch)
View(evaluated_values)
is.fd(combined_fd_garch)
#Conducting FPCA
fpca_results_garch <- pca.fd(fdobj = combined_fd_garch, nharm = 2)
fpca_results_garch$harmonics #a functional data object for the harmonics or eigen-functions
fpca_results_garch$scores #s matrix of scores on the principal components or harmonics
fpca_results_garch$values #the complete set of eigenvalues
fpca_results_garch$varprop #a vector giving the proportion of variance explained by each eigen-function
fpca_results_garch$meanfd #a functional data object giving the mean function
plot.pca.fd(fpca_results_garch)

#Plotting the harmonics individually (functional principal components)
#Plotting the first harmonic
plot(fpca_results_garch$harmonics[1], main = "First Harmonic (PC1)", xlab = "Date", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180), 
     labels = as.Date(seq(min(Prices$Date), max(Prices$Date), by = "180 days")), 
     las = 2)

#Plotting the second harmonic
plot(fpca_results_garch$harmonics[2], main = "Second Harmonic (PC2)", xlab = "Date", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180), 
     labels = as.Date(seq(min(Prices$Date), max(Prices$Date), by = "180 days")), 
     las = 2)
# Convert FPCA scores to a data frame
fpc_scores_garch <- as.data.frame(fpca_results_garch$scores)

# Rename the columns (assuming you extracted 2 principal components/harmonics)
colnames(fpc_scores_garch) <- c("PC1", "PC2")

# Add a column for the asset names. Here, we assume the order of scores 
# corresponds to the columns in `Prices_data_numeric`.
fpc_scores_garch$Asset <- colnames(Prices_data_numeric)

# Reorder the columns to have the Asset name first
fpc_scores_garch <- fpc_scores_garch[, c("Asset", "PC1", "PC2")]

# Print the FPCA scores
print(fpc_scores_garch)

# Optionally, view the scores in a neat table (if using RStudio)
View(fpc_scores_garch)# Plot the FPCA scores in a scatter plot
library(ggplot2)
ggplot(fpc_scores_garch, aes(x = PC1, y = PC2, label = Asset)) +
  geom_point(color = "blue", size = 3) +             
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +    
  labs(title = "FPCA Scores Scatter Plot",
       x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()























