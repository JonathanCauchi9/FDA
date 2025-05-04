install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("fda")
install.packages("imputeTS")
install.packages("zoo")

rm(list = ls())

library(readxl)
library(dplyr)
library(tidyr)
library(fda)
library(imputeTS)
library(zoo)

Radiation_Data <- read_excel("C:/Users/jonat/OneDrive/Desktop/Thesis/Data Files/Global Radiation.xlsx")
View(Radiation_Data)

Radiation_Data <- Radiation_Data %>%
  pivot_wider(
    names_from = Station,
    values_from = `Global radiation (in J/cm2)`,
    id_cols = Date
  )

View(Radiation_Data)
Radiation_Data$Date <- as.Date(as.character(Radiation_Data$Date), format = "%Y%m%d")
time_points_radiation <- as.numeric(Radiation_Data$Date - min(Radiation_Data$Date))
min_time <- min(time_points_radiation)
max_time <- max(time_points_radiation)

Radiation_Data <- Radiation_Data[,-c(1,5,26)]
View(Radiation_Data)
Radiation_numeric <- as.matrix(Radiation_Data)
View(Radiation_numeric)

na_by_col <- colSums(is.na(Radiation_Data))
print("Number of NA values per column:")
print(na_by_col)

Radiation_Data_imp <- Radiation_Data %>%
  mutate(
    across(
      where(is.numeric),
      ~ na_interpolation(.x, option = "spline")
    )
  )
View(Radiation_Data_imp)

#Creating the Fourier Basis:
n_basis_fourier <- 7
basis_fourier <- create.fourier.basis(rangeval = c(min_time, max_time), nbasis = n_basis_fourier)

PHI_fourier <- eval.basis(time_points_radiation, basis_fourier)
dim(PHI_fourier)

Radiation_Data_imp_numeric <- as.matrix(Radiation_Data_imp)
n_stations <- ncol(Radiation_Data_imp_numeric)
n_time <- nrow(Radiation_Data_imp_numeric)
n_time

# We now define a weight vector based on a centered rolling variance.
# For each station, the weight vector will be computed from a 25-day window.
# (For a centered window of 25 days, we take 12 days before and 12 days after the current day.)

# Defining lambda range for the smoothing parameter search:
lambda_range_weighted <- seq(1e-8, 1000, length.out = 50)

optimal_lambdas_radiation_weighted <- data.frame(
  Station = colnames(Radiation_Data_imp_numeric),
  Optimal_lambda_weighted = NA,
  Min_GCV_weighted = NA
)
par(mfrow = c(3, 2))
# Prepare the smoothing parameters using the second derivative penalty:
n_order <- 4
Sec_derivative <- int2Lfd(max(0, n_order - 2))

wvar_list <- list()

for (station in colnames(Radiation_Data_imp_numeric)) {
  cat(paste("Calculating metrics for station (weighted):", station, "\n"))
  
  y <- Radiation_Data_imp_numeric[, station]
  
  wvar <- rollapply(y, width = 505, FUN = var, align = "center", fill = NA)
  wvar[is.na(wvar)] <- var(y, na.rm = TRUE)
  
  weight_vector <- 1 / wvar
  wvar_list[[station]] <- weight_vector
  
  gcvsave_weighted <- numeric(length(lambda_range_weighted))
  
  for (i in seq_along(lambda_range_weighted)) {
    lambda <- lambda_range_weighted[i]
    fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = lambda)
    
    smooth_result <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj, wtvec = weight_vector)
    
    gcvsave_weighted[i] <- smooth_result$gcv
  }
  
  min_gcv_index <- which.min(gcvsave_weighted)
  optimal_lambda_weighted <- lambda_range_weighted[min_gcv_index]
  min_gcv_weighted <- gcvsave_weighted[min_gcv_index]
  
  optimal_lambdas_radiation_weighted[optimal_lambdas_radiation_weighted$Station == station, ] <- list(
    Station = station,
    Optimal_lambda_weighted = optimal_lambda_weighted,
    Min_GCV_weighted = min_gcv_weighted
  )
  
  plot(lambda_range_weighted, gcvsave_weighted, type = "b", col = "blue", lwd = 2, pch = 19,
       xlab = "Lambda", ylab = "GCV", 
       main = paste("GCV Curve for station", station))
  abline(v = optimal_lambda_weighted, col = "red", lty = 2, lwd = 1.5)
}

print(optimal_lambdas_radiation_weighted)
cat("Descriptive statistics for Optimal Weighted Lambdas (Rolling-Variance Weights):\n")
print(summary(optimal_lambdas_radiation_weighted$Optimal_lambda_weighted))
cat("Standard Deviation:", sd(optimal_lambdas_radiation_weighted$Optimal_lambda_weighted), "\n\n")

cat("Descriptive statistics for Minimum Weighted GCVs (Rolling-Variance Weights):\n")
print(summary(optimal_lambdas_radiation_weighted$Min_GCV_weighted))
cat("Standard Deviation:", sd(optimal_lambdas_radiation_weighted$Min_GCV_weighted), "\n\n")

n_stations <- nrow(optimal_lambdas_radiation_weighted)
cat("Total number of stations:", n_stations, "\n")

fitted_values_matrix_radiation_weight <- matrix(NA, nrow = n_time, ncol = n_stations)
smoothed_stations_weighted <- list()

for (station in 1:n_stations) {
  station_name <- colnames(Radiation_Data_imp_numeric)[station]
  cat(paste("Smoothing station:", station_name, "\n"))
  
  optimal_lambda <- optimal_lambdas_radiation_weighted$Optimal_lambda_weighted[station]
  
  y <- Radiation_Data_imp_numeric[, station]
  
  wvar <- rollapply(y, width = 505, FUN = var, align = "center", fill = NA)
  wvar[is.na(wvar)] <- var(y, na.rm = TRUE)
  
  weight_vector <- 1 / wvar
  
  fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = optimal_lambda)
  
  smooth_result_weight <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj, wtvec = weight_vector)
  
  smoothed_stations_weighted[[station_name]] <- smooth_result_weight$fd
  fitted_values_matrix_radiation_weight[, station] <- eval.fd(time_points_radiation, smooth_result_weight$fd)
}

fitted_values_df <- as.data.frame(fitted_values_matrix_radiation_weight)
colnames(fitted_values_df) <- colnames(Radiation_Data_imp_numeric)
View(fitted_values_df)

combined_fd_radiation <- fd(
  coef = do.call(cbind, lapply(smoothed_stations_weighted, function(fd) fd$coefs)),
  basisobj = basis_fourier
)

#Plotting all smoothed functions
plot(
  combined_fd_radiation,
  xlab = "Days", ylab = "Smoothed Value",
  main = "Smoothed Functions for all stations"
)

#Performing Functional Principal Component Analysis
weighted_fpca_results_radiation <- pca.fd(fdobj = combined_fd_radiation, nharm = 3)
plot.pca.fd(weighted_fpca_results_radiation)
weighted_fpca_results_radiation$varprop

#Plotting the first and second harmonics 
par(mfrow = c(1, 1))
plot(weighted_fpca_results_radiation$harmonics[1], main = "First Harmonic (PC1)", xlab = "Days", ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180))
plot(weighted_fpca_results_radiation$harmonics[2], main = "Second Harmonic (PC2)", xlab = "Days", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180))
par(mfrow = c(1, 1))

fpc_scores <- weighted_fpca_results_radiation$scores

#Creating a scatter plot of the first two principal component scores
plot(
  fpc_scores[, 1],   
  fpc_scores[, 2],   
  xlab = "PC1 Score",
  ylab = "PC2 Score",
  main = "FPCA Scores Scatter Plot",
  pch = 19,          
  col = "blue"       
)

grid()

text(
  fpc_scores[, 1],
  fpc_scores[, 2],
  labels = colnames(Radiation_Data_imp_numeric),  
  pos = 3,      
  cex = 0.8,    
  col = "black"
)

