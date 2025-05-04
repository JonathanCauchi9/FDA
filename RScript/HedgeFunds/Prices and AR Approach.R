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

spx_data <- pricesdata %>% 
  dplyr::select(Date, SPX.L)

View(spx_data)

Prices<-pricesdata
View(Prices)
dev.off()

#Data Preparation 
Prices$Date <- as.Date(Prices$Date)
time_points <- as.numeric(Prices$Date - min(Prices$Date)) 
time_points
length(time_points)
max(time_points)

Prices_data <- Prices %>% dplyr::select(-Date)
View(Prices_data)



#Checking for Periodicity
set.seed(16)  
sample_size <- 6  
sampled_assets <- sample(colnames(Prices_data), sample_size)

par(mfrow = c(ceiling(sample_size / 2), 2))  

for (asset in sampled_assets) {
  main_title <- paste("ACF of", asset, "Prices")
  
  acf(Prices_data[[asset]], main = main_title, lag.max = 50)
}

par(mfrow = c(1, 1))



#Step 2: Creating the B-spline Basis
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
PHI  

par(mfrow = c(1, 1))
#Plotting the basis functions and placement of the knots 
matplot(time_points,PHI,type='l',lwd=1,lty=1, xlab='days',ylab='B-spline function values',cex.lab=1,cex.axis=1)
for (i in 1:n_knots)
{
  abline(v=knots[i], lty=2, lwd=1)
}

#Finding the optimal number of basis functions for the fourier basis
K_values <- seq(4, 100, by = 2)  

results_df <- data.frame(
  Asset = character(),
  K     = numeric(),
  s2    = numeric(),
  stringsAsFactors = FALSE
)

n <- length(time_points)
assets <- colnames(Prices_data)

for (asset_name in assets) {
  
  Y <- as.numeric(Prices_data[[asset_name]])
  
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
  labs(title = paste("Residual variance vs K for asset:", assets[1]))

#By visualizing the elbow point in the plots, we deduce that 22 basis functions is adequate 

#Creating the Fourier Basis:
n_basis_fourier <- 23
basis_fourier <- create.fourier.basis(
  rangeval = c(min_time, max_time),
  nbasis   = n_basis_fourier
)

PHI_fourier <- eval.basis(time_points, basis_fourier)
dim(PHI_fourier)

mse_df <- data.frame(
  Asset         = colnames(Prices_data),
  MSE_Bspline   = NA_real_,
  MSE_Fourier   = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(colnames(Prices_data))) {
  
  Y <- as.numeric(Prices_data[[i]])
  
  coeff_bspline <- solve(t(PHI) %*% PHI, t(PHI) %*% Y)
  
  Yhat_bspline <- PHI %*% coeff_bspline
  
  mse_bspline <- mean((Y - Yhat_bspline)^2)
  
  coeff_fourier <- solve(t(PHI_fourier) %*% PHI_fourier, t(PHI_fourier) %*% Y)
  Yhat_fourier <- PHI_fourier %*% coeff_fourier
  
  mse_fourier <- mean((Y - Yhat_fourier)^2)
  
  mse_df$MSE_Bspline[i] <- mse_bspline
  mse_df$MSE_Fourier[i] <- mse_fourier
}

print(mse_df)
mse_df$Winner <- ifelse(mse_df$MSE_Fourier < mse_df$MSE_Bspline, 
                        "Fourier", 
                        ifelse(mse_df$MSE_Fourier > mse_df$MSE_Bspline,
                               "B-spline","Tie"))
table(mse_df$Winner)

# Data preparation
Prices_data_numeric <- as.matrix(Prices_data) 
View(Prices_data_numeric)
n_assets <- ncol(Prices_data_numeric)          
n_assets
n_time <- nrow(Prices_data_numeric)            
n_time

Sec_derivative <- int2Lfd(max(0, n_order - 2))

residuals_matrix <- matrix(NA, nrow = n_time, ncol = n_assets) 
fitted_values_matrix <- matrix(NA, nrow = n_time, ncol = n_assets) 

for (asset in 1:n_assets) {
  cat(paste("Performing OLS for asset", asset, "\n"))
  
  y <- Prices_data_numeric[, asset]
  
  Phi <- eval.basis(time_points, basis)  
  
  #Perform OLS estimation using the normal equation: c_hat = (Phi'Phi)^(-1) Phi'Y
  c_hat <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% y
  
  #Compute the fitted values: Y_hat = Phi * c_hat
  fitted_values <- Phi %*% c_hat
  
  #Calculate residuals: Observed - Fitted
  residuals <- y - fitted_values
  
  #Store results in matrices
  residuals_matrix[, asset] <- residuals  
  fitted_values_matrix[, asset] <- fitted_values  
}

View(residuals_matrix)  
View(fitted_values_matrix)
View(Prices_data_numeric)


#Checking whether there actually is auto-correlation between the residuals through numerous checks

#PACF Checks
library(forecast)

n_assets <- ncol(residuals_matrix)

par(
  mfrow = c(3, 2),
  mar   = c(4, 4, 5, 1),   
  cex.main = 0.9,          
  ask  = FALSE
)

for (i in seq_len(n_assets)) {
  asset_name <- colnames(Prices_data_numeric)[i]
  
  pacf(
    residuals_matrix[, i],
    main    = asset_name,  
    lag.max = 50
  )
}

par(mfrow = c(1,1), mar = c(5,4,4,2) + 0.1)

#The plots suggest fitting an AR(1) model as suspected

#Conducting the Ljung-Box Test to check the autocorrelation in the residuals of all assets
#Checking the null hypothesis of no autocorrelation up to a certain lag


# Number of lags to check for autocorrelation
lag_k <- 1

lb_pvalues <- numeric(n_assets)
names(lb_pvalues) <- colnames(Prices_data_numeric)

for (asset in seq_len(n_assets)) {
  asset_name <- colnames(Prices_data_numeric)[asset]
  
  #Performing the Ljung-Box test on the residuals of the current asset
  lb_test <- Box.test(
    residuals_matrix[, asset], 
    lag  = lag_k, 
    type = "Ljung-Box"
  )
  
  lb_pvalues[asset_name] <- lb_test$p.value
  
  cat(sprintf("Ljung-Box Test for Asset: %s — p-value = %.4f\n", 
              asset_name, lb_test$p.value))
}

print(lb_pvalues)
#P values are smaller than 0.05 for all assets so we accept H1: The auto-correlations up to lag 1 are not zero

n_assets <- ncol(residuals_matrix)
n_time   <- nrow(residuals_matrix)

par(
  mfrow    = c(3, 2),        
  mar      = c(4, 4, 5, 1),  
  cex.main = 0.9,            
  ask      = FALSE           
)

for (i in seq_len(n_assets)) {
  asset_name <- colnames(Prices_data_numeric)[i]
  
  plot(
    residuals_matrix[1:(n_time - 1), i],
    residuals_matrix[2:n_time,     i],
    xlab = "Residual(t)",
    ylab = "Residual(t+1)",
    main = paste("Scatter of residuals —", asset_name)
  )
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)


#So fit an AR(1) Model to the residuals 

#Fitting the AR(1) model and estimating the variance covariance matrix of the residuals

library(forecast)

n_assets <- ncol(residuals_matrix)
n_time <- nrow(residuals_matrix)

rho_values <- numeric(n_assets)    
sigma2_values <- numeric(n_assets) 

cov_matrix <- vector("list", n_assets)  

wtvec_list <- list()  

for (asset in 1:n_assets) {
  tryCatch({
    ar1_model <- arima(residuals_matrix[, asset], order = c(1, 0, 0))  
    
    rho_values[asset] <- ar1_model$coef[1]  
    sigma2_values[asset] <- ar1_model$sigma2  
    
    cat(paste("Asset", asset, "AR(1) Model Parameters:\n"))
    cat("Autocorrelation parameter (rho):", rho_values[asset], "\n")
    cat("Variance of noise (sigma^2):", sigma2_values[asset], "\n")
    
    asset_cov_matrix <- matrix(NA, nrow = n_time, ncol = n_time)  
    
    for (t1 in 1:n_time) {
      for (t2 in 1:n_time) {
        h <- abs(t1 - t2)
        
        asset_cov_matrix[t1, t2] <- sigma2_values[asset] * (rho_values[asset]^h)
      }
    }
    
    cov_matrix[[asset]] <- asset_cov_matrix
    
    wtvec_list[[asset]] <- solve(asset_cov_matrix)
    
    cat(paste("Covariance Matrix for Asset", asset, ":\n"))
    print(asset_cov_matrix)
  }, error = function(e) {
    cat(paste("Error fitting AR(1) model for Asset", asset, ": ", e$message, "\n"))
  })
}

parameter_estimates <- data.frame(
  Asset = 1:n_assets,
  Rho = rho_values,
  Sigma2 = sigma2_values
)

print(parameter_estimates)

lambda_range <- seq(1e-8, 1e4, length.out = 50)

optimal_lambdas_gcv <- data.frame(
  Asset = colnames(Prices_data_numeric),
  Optimal_lambda_GCV = NA,
  Min_GCV = NA,
  stringsAsFactors = FALSE
)

for (asset in colnames(Prices_data_numeric)) {
  cat(paste("Calculating GCV for asset:", asset, "\n"))
  
  y <- Prices_data_numeric[, asset]
  wtvec <- wtvec_list[[asset]]
  
  gcvsave <- numeric(length(lambda_range))
  
  for (i in seq_along(lambda_range)) {
    lambda <- lambda_range[i]  
    fdParobj <- fdPar(basis, Sec_derivative, lambda = lambda)  
    
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj,wtvec = wtvec)
    
    gcvsave[i] <- smooth_result$gcv  
  }
  
  min_gcv_index <- which.min(gcvsave)
  optimal_lambda_gcv <- lambda_range[min_gcv_index]
  min_gcv <- gcvsave[min_gcv_index]
  
  optimal_lambdas_gcv[optimal_lambdas_gcv$Asset == asset, ] <- list(
    Asset = asset,
    Optimal_lambda_GCV = optimal_lambda_gcv,
    Min_GCV = min_gcv
  )
  
  plot(lambda_range, gcvsave, type = "b", col = "blue", lwd = 2, pch = 19,
       xlab = "Lambda", ylab = "GCV", 
       main = paste("GCV Curve for Asset", asset))
  
  abline(v = optimal_lambda_gcv, col = "red", lty = 2, lwd = 1.5)
}

print(optimal_lambdas_gcv)


smoothed_assets <- list()

#Looping through each asset to apply smoothing using optimal lambdas
for (asset in optimal_lambdas_gcv$Asset) {
  cat(paste("Smoothing asset:", asset, "\n"))
  
  tryCatch({
    optimal_lambda <- optimal_lambdas_gcv[optimal_lambdas_gcv$Asset == asset, "Optimal_lambda_GCV"]
    
    y <- Prices_data_numeric[, asset]
    wtvec <- wtvec_list[[asset]]
    
    fdParobj <- fdPar(basis, Sec_derivative, lambda = optimal_lambda)
    
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj, wtvec = wtvec)
    
    smoothed_assets[[asset]] <- smooth_result$fd
  }, error = function(e) {
    cat(paste("Error smoothing asset:", asset, "- Skipping.\n"))
  })
}

coef_matrix <- sapply(smoothed_assets, function(fd_obj) fd_obj$coefs)
combined_fd <- fd(coef = coef_matrix, basisobj = basis)

#Plot all smoothed functions using plot.fd
plot(
  combined_fd,
  xlab = "Date", ylab = "Smoothed Value",
  main = "Smoothed Functions for All Assets"
)
#Conducting FPCA
fpca_results_AR <- pca.fd(fdobj = combined_fd, nharm = 2)
fpca_results_AR$harmonics 
fpca_results_AR$scores 
fpca_results_AR$values
fpca_results_AR$varprop 
fpca_results_AR$meanfd 
plot.pca.fd(fpca_results_AR)


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

#Convert FPCA scores to a data frame
fpc_scores <- as.data.frame(fpca_results_AR$scores)

colnames(fpc_scores) <- c("PC1", "PC2")

fpc_scores$Asset <- colnames(Prices_data_numeric)

fpc_scores <- fpc_scores[, c("Asset", "PC1", "PC2")]

print(fpc_scores)

View(fpc_scores)

#Plot the FPCA scores in a scatter plot
library(ggplot2)
ggplot(fpc_scores, aes(x = PC1, y = PC2, label = Asset)) +
  geom_point(color = "blue", size = 3) +             
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +    
  labs(title = "FPCA Scores Scatter Plot",
       x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()

n_assets <- ncol(Prices_data_numeric)

#Derivative Estimation
#Compute the first derivatives of the smoothed functions
combined_fd_deriv <- deriv.fd(combined_fd, Lfdobj = 1)
par(
  mfrow    = c(3, 2),        
  mar      = c(4, 4, 5, 1),  
  cex.main = 0.9,            
  ask      = FALSE           
)

for (i in seq_len(n_assets)) {
  asset_name <- colnames(Prices_data_numeric)[i]
  cat(sprintf("Plotting derivative for asset %d: %s\n", i, asset_name))
  
  plot(
    combined_fd_deriv[i],
    xlab = "Time (Days)",
    ylab = "First Derivative",
    main = asset_name
  )
}


