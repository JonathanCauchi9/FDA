# Load necessary libraries
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("fda")
install.packages("imputeTS")
install.packages("ggplot2")
install.packages("forecast")
install.packages("seastests")
install.packages("tseries")


library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(fda)
library(imputeTS)
library(forecast)
library(seastests)
library(tseries)
rm(list = ls())

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
radts<-ts(Radiation_Data[,3])
plot.ts(radts)
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

# Checking for Periodicity 

set.seed(9)  
sample_size <- 6  
sampled_stations <- sample(colnames(Radiation_Data_imp), sample_size)

par(mfrow = c(ceiling(sample_size / 2), 2))  

# Loop through the sampled stations to plot the ACF
for (station in sampled_stations) {
  # Title for the plot
  main_title <- paste("ACF of station", station)
  
  # Plot the ACF for the station
  acf(Radiation_Data_imp[[station]], main = main_title, lag.max = 20)
}

stl_results <- list()  # Store the STL results for later inspection

par(mfrow = c(2, 2))  # Adjust rows/cols to your preference

for (station in colnames(Radiation_Data_imp)) {
  ts_data <- ts(Radiation_Data_imp[[station]], frequency = 365)
  
  # Perform STL decomposition
  stl_decomp <- stl(ts_data, s.window = "periodic")
  
  # Store results in a list if you want to inspect them later
  stl_results[[station]] <- stl_decomp
  
  # Plot the decomposition
  plot_title <- paste("STL Decomposition -", station)
  plot(stl_decomp, main = plot_title)
}
par(mfrow = c(1, 1))

#Checking stationarity of the data
stationarity_results <- data.frame(Station = character(), P_Value = numeric(), IsStationary = logical(), stringsAsFactors = FALSE)

# Loop through each station (column) to test for stationarity
for (station in colnames(Radiation_Data_imp)) {
  # Convert the data to a time series object
  ts_data <- ts(Radiation_Data_imp[[station]], frequency = 365)
  
  # Perform ADF test
  adf_test <- adf.test(ts_data, alternative = "stationary", k = 0)
  
  # Extract p-value
  p_value <- adf_test$p.value
  
  # Check if the data is stationary (p-value < 0.05 indicates stationarity)
  is_stationary <- p_value < 0.05
  
  # Append results to the data frame
  stationarity_results <- rbind(stationarity_results, 
                                data.frame(Station = station, 
                                           P_Value = p_value, 
                                           IsStationary = is_stationary))
}

# View the results
print(stationarity_results)

# Initialize a data frame to store the seasonality test results
seasonality_results <- data.frame(
  Station = character(),
  IsSeasonal = logical(),
  stringsAsFactors = FALSE
)

# Loop through each station to test for seasonality
for (station in colnames(Radiation_Data_imp)) {
  # Convert the station data to a time series object with daily frequency
  ts_data <- ts(Radiation_Data_imp[[station]], frequency = 365)
  
  # Check for seasonality using isSeasonal
  seasonal_check <- isSeasonal(ts_data)
  
  # Print a message for immediate feedback
  if (seasonal_check) {
    cat("Seasonality detected for station:", station, "\n")
  } else {
    cat("No significant seasonality detected for station:", station, "\n")
  }
  
  # Store the result in the data frame
  seasonality_results <- rbind(
    seasonality_results,
    data.frame(Station = station, IsSeasonal = seasonal_check, stringsAsFactors = FALSE)
  )
}

# View the results
print(seasonality_results)

#Step 2: Creating the B-spline Basis
#Define the time range
time_points_radiation
min_time <- min(time_points_radiation)
max_time <- max(time_points_radiation)
min_time
max_time

#Define knots and basis parameters
knots <- c(seq(min_time,max_time,137)) #Set knots every 137 days
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
PHI = eval.basis(time_points_radiation, basis) 
dim(PHI)
PHI #1097x11 Matrix 

par(mfrow = c(1, 1))
#Plotting the basis functions and placement of the knots 
matplot(time_points_radiation,PHI,type='l',lwd=1,lty=1, xlab='days',ylab='B-spline function values',cex.lab=1,cex.axis=1)
for (i in 1:n_knots)
{
  abline(v=knots[i], lty=2, lwd=1)}

#Finding the optimal number of basis functions for the fourier basis
K_values <- seq(4, 100, by = 2)  

# We'll create an empty data frame to store the results
results_df <- data.frame(
  station = character(),
  K     = numeric(),
  s2    = numeric(),
  stringsAsFactors = FALSE
)

n <- length(time_points_radiation)
stations <- colnames(Radiation_Data_imp)
View(Radiation_Data_imp)

for (station_name in stations) {
  Y <- as.numeric(Radiation_Data_imp[[station_name]])
  
  # 2b) Loop over each K in K_values
  for (k in K_values) {
    
    fourier_basis <- create.fourier.basis(
      rangeval = c(min(time_points_radiation), max(time_points_radiation)),
      nbasis   = k
    )
    
    Phi <- eval.basis(time_points_radiation, fourier_basis)
    c_hat <- solve(t(Phi) %*% Phi, t(Phi) %*% Y)
    
    Y_hat <- Phi %*% c_hat
    RSS <- sum((Y - Y_hat)^2)
    s2_val <- RSS / (n - k)
    
    results_df <- rbind(
      results_df,
      data.frame(
        Station = station_name,
        K     = k,
        s2    = s2_val,
        stringsAsFactors = FALSE
      )
    )
  }
}

print(results_df)
best_k_df <- results_df %>%
  group_by(Station) %>%
  slice_min(s2)

print(best_k_df)
ggplot(
  data = subset(results_df, Station == stations[1]),
  aes(x = K, y = s2)
) +
  geom_line() +
  geom_point() +
  labs(title = paste("Residual variance vs K for station:", stations[23]))


#By visualizing the elbow point in the plots, we deduce that 6 basis functions is adequate 

#Creating the Fourier Basis:
n_basis_fourier <- 7
basis_fourier <- create.fourier.basis(
  rangeval = c(min_time, max_time),
  nbasis   = n_basis_fourier
)

# Evaluate Fourier basis functions at the time points
PHI_fourier <- eval.basis(time_points_radiation, basis_fourier)
dim(PHI_fourier)

# Prepare a data frame to store MSE for each station
mse_df <- data.frame(
  Station         = colnames(Radiation_Data_imp),
  MSE_Bspline   = NA_real_,
  MSE_Fourier   = NA_real_,
  stringsAsFactors = FALSE
)

# Loop over columns (stations)
for (i in seq_along(colnames(Radiation_Data_imp))) {
  
  # Extract the i-th station's time series as a numeric vector
  Y <- as.numeric(Radiation_Data_imp[[i]])
  
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
# Suppose 'mse_df' has columns: station, MSE_Bspline, MSE_Fourier
mse_df$Winner <- ifelse(mse_df$MSE_Fourier < mse_df$MSE_Bspline, 
                        "Fourier", 
                        ifelse(mse_df$MSE_Fourier > mse_df$MSE_Bspline,
                               "B-spline","Tie"))
table(mse_df$Winner)

#So in this case, a fourier basis is more optimal - Hinting at the periodicity of the data. 

Radiation_Data_imp_numeric<- as.matrix(Radiation_Data_imp)
n_stations<-ncol(Radiation_Data_imp_numeric)
n_stations
# Define lambda range
# Unweighted Analysis
lambda_range_unweighted <- seq(1e-8, 11000, length.out = 50)  # Lambda values to test
Sec_derivative <- int2Lfd(max(0, n_order - 2))

# Initialize storage for results
optimal_lambdas_radiation_unweighted <- data.frame(
  Station = colnames(Radiation_Data_imp_numeric),
  Optimal_lambda_unweighted = NA,
  Min_GCV_unweighted = NA
)

par(mfrow = c(1, 1))

# Loop through each station for unweighted analysis
for (station in 1:n_stations) {
  station_name <- colnames(Radiation_Data_imp_numeric)[station]
  cat(paste("Processing station (unweighted):", station_name, "\n"))
  
  y <- Radiation_Data_imp_numeric[, station]
  
  # Initialize GCV storage
  gcvsave_unweighted <- numeric(length(lambda_range_unweighted))
  
  # Loop through lambda values
  for (i in seq_along(lambda_range_unweighted)) {
    lambda <- lambda_range_unweighted[i]
    fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda)
    smooth_result <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj)
    gcvsave_unweighted[i] <- smooth_result$gcv  # GCV for this lambda
  }
  
  # Find the optimal lambda with the minimum GCV
  min_gcv_index <- which.min(gcvsave_unweighted)
  optimal_lambda_unweighted <- lambda_range_unweighted[min_gcv_index]
  min_gcv_unweighted <- gcvsave_unweighted[min_gcv_index]
  
  # Store the results
  optimal_lambdas_radiation_unweighted[station, ] <- list(
    Station = station_name,
    Optimal_lambda_unweighted = optimal_lambda_unweighted,
    Min_GCV_unweighted = min_gcv_unweighted
  )
  
  # Plot the GCV curve for the current station
  plot(
    lambda_range_unweighted, gcvsave_unweighted, type = "b", col = "blue", lwd = 2, pch = 19,
    xlab = "Lambda", ylab = "GCV",
    main = paste("GCV Curve for\n", strwrap(station_name, width = 30), collapse = "\n")
  )
  abline(v = optimal_lambda_unweighted, col = "red", lty = 2, lwd = 1.5)  # Highlight optimal lambda
}

print(optimal_lambdas_radiation_unweighted)
optimal_lambda_values <- as.numeric(optimal_lambdas_radiation_unweighted$Optimal_lambda_unweighted)
min_gcv_values <- as.numeric(optimal_lambdas_radiation_unweighted$Min_GCV_unweighted)

# Print descriptive statistics for optimal lambdas
cat("Descriptive statistics for Optimal Lambdas (Unweighted Analysis):\n")
print(summary(optimal_lambda_values))
cat("Standard Deviation:", sd(optimal_lambda_values), "\n\n")

# Print descriptive statistics for minimum GCV values
cat("Descriptive statistics for Minimum GCVs (Unweighted Analysis):\n")
print(summary(min_gcv_values))
cat("Standard Deviation:", sd(min_gcv_values), "\n")
n_time <- nrow(Radiation_Data_imp_numeric)



# Store fitted values matrix for unweighted analysis
fitted_values_matrix_radiation_unweight <- matrix(NA, nrow = n_time, ncol = n_stations)

# Smooth all stations with their respective optimal lambdas
smoothed_stations_unweighted <- list()
for (station in 1:n_stations) {
  station_name <- colnames(Radiation_Data_imp_numeric)[station]
  cat(paste("Smoothing station:", station_name, "\n"))
  optimal_lambda <- optimal_lambdas_radiation_unweighted$Optimal_lambda[station]
  y <- Radiation_Data_imp_numeric[, station]
  fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = optimal_lambda)
  smooth_result_un <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj)
  smoothed_stations_unweighted[[station_name]] <- smooth_result_un$fd
  fitted_values_matrix_radiation_unweight[, station] <- eval.fd(time_points_radiation, smooth_result_un$fd)
}
par(mfrow = c(1, 1))

# Convert the fitted values matrix to a dataframe for easier inspection
fitted_values_df_un <- as.data.frame(fitted_values_matrix_radiation_unweight)
colnames(fitted_values_df_un) <- colnames(Radiation_Data_imp_numeric)

# View the fitted values
View(fitted_values_df_un)

# Combine smoothed functions into a single functional data object for plotting
combined_fd_radiation <- fd(
  coef = do.call(cbind, lapply(smoothed_stations_unweighted, function(fd) fd$coefs)),
  basisobj = basis_fourier
)

# Plot all smoothed functions
plot(
  combined_fd_radiation,
  xlab = "Time", ylab = "Smoothed Value",
  main = "Smoothed Functions for all stations"
)

#Performing Functional Principal Component Analysis (FPCA) 
unweighted_fpca_results_radiation <- pca.fd(fdobj = combined_fd_radiation, nharm = 3)
unweighted_fpca_results_radiation$harmonics #a functional data object for the harmonics or eigen-functions
unweighted_fpca_results_radiation$scores #s matrix of scores on the principal components or harmonics
unweighted_fpca_results_radiation$values #the complete set of eigenvalues
unweighted_fpca_results_radiation$varprop #a vector giving the proportion of variance explained by each eigen-function
unweighted_fpca_results_radiation$meanfd #a functional data object giving the mean function
plot.pca.fd(unweighted_fpca_results_radiation)
unweighted_fpca_results_radiation$varprop

par(mfrow = c(1, 1))
#Plotting the first harmonic
plot(unweighted_fpca_results_radiation$harmonics[1], main = "First Harmonic (PC1)", xlab = "Date", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

#Plotting the second harmonic
plot(unweighted_fpca_results_radiation$harmonics[2], main = "Second Harmonic (PC2)", xlab = "Date", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

fpc_scores <- as.data.frame(unweighted_fpca_results_radiation$scores)

# Rename the columns (assuming you have 2 principal components/harmonics)
colnames(fpc_scores) <- c("PC1", "PC2")

# Add a column for the station names. 
# Here, we assume that the order of the scores corresponds to the order of the stations 
# in your imputed data matrix.
fpc_scores$Station <- colnames(Radiation_Data_imp_numeric)

# Reorder the columns to have the Station name first
fpc_scores <- fpc_scores[, c("Station", "PC1", "PC2")]

# Print the FPCA scores with station numbers
print(fpc_scores)

# Optionally, view the scores in a neat table (if using RStudio)
View(fpc_scores)

ggplot(fpc_scores, aes(x = PC1, y = PC2, label = Station)) +
  geom_point(color = "blue", size = 3) +         # Plot the points
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) + # Add station labels near points
  labs(title = "FPCA Scores Scatter Plot",
       x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()


#Conducting VARIMAX Rotation
unweighted_fpca_resultsrot_radiation <- varmx.pca.fd(unweighted_fpca_results_radiation, nharm = 2)
unweighted_fpca_resultsrot_radiation$harmonics #a functional data object for the harmonics or eigen-functions
unweighted_fpca_resultsrot_radiation$scores #s matrix of scores on the principal components or harmonics
unweighted_fpca_resultsrot_radiation$values #the complete set of eigenvalues
unweighted_fpca_resultsrot_radiation$varprop #a vector giving the proportion of variance explained by each eigen-function
unweighted_fpca_resultsrot_radiation$meanfd #a functional data object giving the mean function
plot.pca.fd(unweighted_fpca_resultsrot_radiation)


#Plotting the first harmonic
plot(unweighted_fpca_resultsrot_radiation$harmonics[1], main = "First Harmonic (PC1) Rot", xlab = "Date", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

#Plotting the second harmonic
plot(unweighted_fpca_resultsrot_radiation$harmonics[2], main = "Second Harmonic (PC2) Rot", xlab = "Date", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

# Compute the first derivatives of the smoothed radiation functions
combined_fd_radiation_deriv <- deriv.fd(combined_fd_radiation, Lfdobj = 1)

# Plot the first derivatives for all stations in one plot
plot(combined_fd_radiation_deriv, 
     xlab = "Time (Days)", 
     ylab = "First Derivative", 
     main = "First Derivative of Smoothed Radiation Functions")

# Alternatively, plot derivatives for a subset of stations only
selected_indices <- c(1, 3, 5)  # Adjust these indices to select desired stations
par(mfrow = c(ceiling(length(selected_indices) / 2), 2))
for (i in selected_indices) {
  station_name <- colnames(Radiation_Data_imp_numeric)[i]
  plot(combined_fd_radiation_deriv[i], 
       xlab = "Time (Days)", 
       ylab = "First Derivative", 
       main = paste("Derivative of", station_name))
}
par(mfrow = c(1, 1))

