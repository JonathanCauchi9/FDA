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

Radiation_Data <- read_excel("C:/Users/jonat/OneDrive/Desktop/Thesis/Data Files/Global Radiation.xlsx")
View(Radiation_Data)

Radiation_Data <- Radiation_Data %>%
  pivot_wider(
    names_from = Station,
    values_from = `Global radiation (in J/cm2)`,
    id_cols = Date
  )

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

Radiation_Data_imp <- Radiation_Data %>%
  mutate(
    across(
      where(is.numeric),              
      ~ na_interpolation(.x, option = "spline")  
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

PHI_fourier <- eval.basis(time_points_radiation, basis_fourier)
dim(PHI_fourier)

Radiation_Data_imp_numeric<- as.matrix(Radiation_Data_imp)
n_stations<-ncol(Radiation_Data_imp_numeric)
n_stations
n_time <- nrow(Radiation_Data_imp_numeric)            
n_time

n_order <- 4                             
Sec_derivative <- int2Lfd(max(0, n_order - 2))
residuals_matrix_radiation <- matrix(NA, nrow = n_time, ncol = n_stations)
fitted_values_matrix_radiation <- matrix(NA, nrow = n_time, ncol = n_stations)


for (station in 1:n_stations) {
  cat(paste("Performing OLS for station", colnames(Radiation_Data_imp_numeric)[station], "\n"))
  
  y <- Radiation_Data_imp_numeric[, station]
  
  Phi <- eval.basis(time_points_radiation, basis_fourier)  
  
  c_hat <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% y
  
  fitted_values <- Phi %*% c_hat
  
  residuals <- y - fitted_values
  
  residuals_matrix_radiation[, station] <- residuals  
  fitted_values_matrix_radiation[, station] <- fitted_values  
}

View(residuals_matrix_radiation)  
View(fitted_values_matrix_radiation)
View(Radiation_Data_imp_numeric)

library(forecast)

n_stations <- ncol(residuals_matrix_radiation)

par(
  mfrow    = c(3, 2),        
  mar      = c(4, 4, 5, 1),  
  cex.main = 0.9,           
  ask      = FALSE           
)

for (station in seq_len(n_stations)) {
  p    <- pacf(residuals_matrix_radiation[, station], plot = FALSE, lag.max = 50)
  lags <- p$lag
  acfs <- p$acf
  idx  <- which.max(abs(acfs))
  
  st_idx  <- station
  st_name <- colnames(Radiation_Data_imp_numeric)[station]
  max_lag <- lags[idx]
  max_acf <- round(acfs[idx], 3)
  
  cat(sprintf("Station %d (%s) — Highest PACF at lag %d = %.3f\n", 
              st_idx, st_name, max_lag, max_acf))
  
  pacf(
    residuals_matrix_radiation[, station],
    main    = sprintf(" %s\nMax lag %d = %.3f", st_name, max_lag, max_acf),
    lag.max = 50
  )
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

#Ljung Box Test
lag_k1 <- 1

station_names       <- if (!is.null(colnames(residuals_matrix_radiation))) {
  colnames(residuals_matrix_radiation)
} else {
  paste0("Station", seq_len(n_stations))
}
lb_pvalues_radiation <- numeric(n_stations)
names(lb_pvalues_radiation) <- station_names

for (station in seq_len(n_stations)) {
  lb_test <- Box.test(
    residuals_matrix_radiation[, station],
    lag  = lag_k1,
    type = "Ljung-Box"
  )
  
  lb_pvalues_radiation[station] <- lb_test$p.value
  
  cat(sprintf(
    "Ljung-Box p-value for %s = %.4g\n",
    station_names[station],
    lb_test$p.value
  ))
}

print(lb_pvalues_radiation)
#P values are smaller than 0.05 for all stations so we accept H1: The auto-correlations up to lag 1 are not zero


par(mfrow = c(1, 1))

n_stations <- ncol(residuals_matrix_radiation)
n_time     <- nrow(residuals_matrix_radiation)

par(
  mfrow    = c(3, 2),        
  mar      = c(4, 4, 5, 1),  
  cex.main = 0.9,            
  ask      = FALSE           
)

for (station in seq_len(n_stations)) {
  st_name <- colnames(Radiation_Data_imp_numeric)[station]
  
  cat(sprintf("Plotting scatter for Station %d: %s\n", station, st_name))
  
  plot(
    residuals_matrix_radiation[1:(n_time - 1), station],
    residuals_matrix_radiation[2:n_time,     station],
    xlab = "Residual(t)",
    ylab = "Residual(t+1)",
    main = sprintf( st_name)
  )
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

#Fitting the AR(1) model and estimating the variance covariance matrix of the residuals

library(forecast)

n_stations <- ncol(residuals_matrix_radiation)
n_time <- nrow(residuals_matrix_radiation)
n_stations
rho_values <- numeric(n_stations)   
sigma2_values <- numeric(n_stations) 

cov_matrix <- vector("list", n_stations)  

#Fitting AR(1) model to each station's residuals and estimate parameters
wtvec_list <- list()  

for (station in 1:n_stations) {
  tryCatch({
    ar1_model <- arima(residuals_matrix_radiation[, station], order = c(1, 0, 0))  # AR(1) model
    
    rho_values[station] <- ar1_model$coef[1]  
    sigma2_values[station] <- ar1_model$sigma2  
    
    cat(paste("station", station, "AR(1) Model Parameters:\n"))
    cat("Autocorrelation parameter (rho):", rho_values[station], "\n")
    cat("Variance of noise (sigma^2):", sigma2_values[station], "\n")
    
    station_cov_matrix <- matrix(NA, nrow = n_time, ncol = n_time)  
    
    for (t1 in 1:n_time) {
      for (t2 in 1:n_time) {
        h <- abs(t1 - t2)
        
        station_cov_matrix[t1, t2] <- sigma2_values[station] * (rho_values[station]^h)
      }
    }
    
    cov_matrix[[station]] <- station_cov_matrix
    
    wtvec_list[[station]] <- solve(station_cov_matrix)
    
    cat(paste("Covariance Matrix for station", station, ":\n"))
    print(station_cov_matrix)
  }, error = function(e) {
    cat(paste("Error fitting AR(1) model for station", station, ": ", e$message, "\n"))
  })
}

parameter_estimates <- data.frame(
  station = 1:n_stations,
  Rho = rho_values,
  Sigma2 = sigma2_values
)
print(parameter_estimates)

lambda_range_weighted <- seq(1e-8, 11000, length.out = 50) 

optimal_lambdas_radiation_weighted <- data.frame(
  Station = colnames(Radiation_Data_imp_numeric),
  Optimal_lambda_weighted = NA,
  Min_GCV_weighted = NA
)

for (station in colnames(Radiation_Data_imp_numeric)) {
  cat(paste("Calculating metrics for station (weighted):", station, "\n"))
  
  y <- Radiation_Data_imp_numeric[, station]
  wtvec <- wtvec_list[[station]]
  
  gcvsave_weighted <- numeric(length(lambda_range_weighted))
  
  for (i in seq_along(lambda_range_weighted)) {
    lambda <- lambda_range_weighted[i]  
    fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = lambda)  
    
    smooth_result <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj, wtvec = wtvec)
    
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

fitted_values_matrix_radiation_weight <- matrix(NA, nrow = n_time, ncol = n_stations)


#Smoothing all stations with their respective optimal lambdas
smoothed_stations_weighted <- list()

for (station in 1:n_stations) {
  station_name <- colnames(Radiation_Data_imp_numeric)[station]
  cat(paste("Smoothing station:", station_name, "\n"))
  
  optimal_lambda <- optimal_lambdas_radiation_weighted$Optimal_lambda[station]
  
  y <- Radiation_Data_imp_numeric[, station]
  
  wtvec <- wtvec_list[[station_name]]
  
  fdParobj <- fdPar(basis_fourier, Sec_derivative, lambda = optimal_lambda)
  
  smooth_result_weight <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdParobj, wtvec = wtvec)
  
  smoothed_stations_weighted[[station_name]] <- smooth_result_weight$fd
  
  fitted_values_matrix_radiation_weight[, station] <- eval.fd(time_points_radiation, smooth_result_weight$fd)
}

fitted_values_df <- as.data.frame(fitted_values_matrix_radiation_weight)
colnames(fitted_values_df) <- colnames(Radiation_Data_imp_numeric)

View(fitted_values_df)

par(mfrow = c(1, 1))


#Combining smoothed functions into a single functional data object for plotting
combined_fd_radiation <- fd(
  coef = do.call(cbind, lapply(smoothed_stations_weighted, function(fd) fd$coefs)),
  basisobj = basis_fourier
)

#Plotting all smoothed functions
plot(
  combined_fd_radiation,
  xlab = "Time", ylab = "Smoothed Value",
  main = "Smoothed Functions for all stations"
)

#Performing Functional Principal Component Analysis
weighted_fpca_results_radiation <- pca.fd(fdobj = combined_fd_radiation, nharm = 2)
weighted_fpca_results_radiation$harmonics 
weighted_fpca_results_radiation$scores 
weighted_fpca_results_radiation$values 
weighted_fpca_results_radiation$varprop 
weighted_fpca_results_radiation$meanfd 
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

#Compute the first derivatives of the smoothed functions
combined_fd_deriv_rad <- deriv.fd(combined_fd_radiation, Lfdobj = 1)

#Plot the first derivatives for all assets in one plot
plot(combined_fd_deriv_rad, 
     xlab = "Time (Days)", 
     ylab = "First Derivative", 
     main = "First Derivative of Smoothed Asset Functions")


