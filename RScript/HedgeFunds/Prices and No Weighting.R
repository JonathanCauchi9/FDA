install.packages("refund")
install.packages("tidyr")
install.packages("tseries")
install.packages("urca")
install.packages("FinTS")
install.packages("rugarch")
install.packages("tseries")
install.packages("seastests")
install.packages("gridExtra")
library(rugarch)
library(gridExtra)
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
Prices$Date <- as.Date(Prices$Date)
time_points <- as.numeric(Prices$Date - min(Prices$Date)) 
time_points
length(time_points)
max(time_points)

Prices_data <- Prices %>% dplyr::select(-Date)
View(Prices_data)



# Checking for Periodicity
assets <- colnames(Prices_data)        

op <- par(no.readonly = TRUE)          
on.exit(par(op), add = TRUE)           
par(ask = FALSE)                       

for (i in seq(1, length(assets), by = 6)) {
  
  grp <- assets[i : min(i + 5, length(assets))]   
  
  plot.new()                                      
  par(mfrow = c(3, 2), mar = c(4.5, 4, 3, 1))     
  
  for (asset in grp) {
    acf(Prices_data[[asset]],
        main    = paste("ACF –", asset),
        xlab    = "Lag",
        ylab    = "ACF",
        lag.max = 50)
  }
  
  if (length(grp) < 6)
    for (k in seq_len(6 - length(grp))) plot.new()
}


seasonality_results <- data.frame(
  Asset = character(),
  IsSeasonal = logical(),
  stringsAsFactors = FALSE
)

#Loop through each asset in Prices_data to test for seasonality
for (asset in colnames(Prices_data)) {
  ts_data <- ts(Prices_data[[asset]], frequency = 365)
  
  # Check for seasonality using isSeasonal
  seasonal_check <- isSeasonal(ts_data)
  
  if (seasonal_check) {
    cat("Seasonality detected for asset:", asset, "\n")
  } else {
    cat("No significant seasonality detected for asset:", asset, "\n")
  }
  
  seasonality_results <- rbind(
    seasonality_results,
    data.frame(Asset = asset, IsSeasonal = seasonal_check, stringsAsFactors = FALSE)
  )
}

print(seasonality_results)

library(tseries)
library(seastests)

assets <- colnames(Prices_data)
p_values <- numeric(length(assets))
is_stationary <- logical(length(assets))

Prices_data_adjusted <- Prices_data

for (i in seq_along(assets)) {
  ts_data <- ts(Prices_data[[assets[i]]], frequency = 365)
  
  adf_test <- adf.test(ts_data, alternative = "stationary", k = 0)
  
  p_values[i] <- adf_test$p.value
  is_stationary[i] <- p_values[i] < 0.05
  
  if (!is_stationary[i]) {
    differenced_data <- diff(ts_data, differences = 1)
    
    differenced_data <- c(differenced_data, tail(differenced_data, 1))
    
    Prices_data_adjusted[[assets[i]]] <- differenced_data
  }
}

stationarity_results <- data.frame(Asset = assets, P_Value = p_values, IsStationary = is_stationary)

print(stationarity_results)


p_values_adjusted <- numeric(length(assets))
is_stationary_adjusted <- logical(length(assets))

for (i in seq_along(assets)) {
  ts_data_adjusted <- ts(Prices_data_adjusted[[assets[i]]], frequency = 365)
  adf_test_adjusted <- adf.test(ts_data_adjusted, alternative = "stationary", k = 0)
  p_values_adjusted[i] <- adf_test_adjusted$p.value
  is_stationary_adjusted[i] <- p_values_adjusted[i] < 0.05
}

stationarity_results_adjusted <- data.frame(Asset = assets, 
                                            P_Value = p_values_adjusted, 
                                            IsStationary = is_stationary_adjusted)

print(stationarity_results_adjusted)


is_seasonal <- logical(length(assets))

for (i in seq_along(assets)) {
  ts_data_seasonal <- ts(Prices_data_adjusted[[assets[i]]], frequency = 365)
  is_seasonal[i] <- isSeasonal(ts_data_seasonal)
}

seasonality_results <- data.frame(Asset = assets, IsSeasonal = is_seasonal)

print(seasonality_results)

#Creating the B-spline Basis
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
    abline(v=knots[i], lty=2, lwd=1)}

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
assets <- unique(results_df$Asset)     

for (i in seq(1, length(assets), by = 6)) {
  
  grp <- assets[i : min(i + 5, length(assets))]      
  
  page_plots <- lapply(grp, function(sym) {
    ggplot(subset(results_df, Asset == sym),
           aes(x = K, y = s2)) +
      geom_line() +
      geom_point(size = 1.4) +
      labs(title = paste("s² vs K — Asset", sym),
           x = "K (Fourier basis size)",
           y = expression(sigma^2)) +
      theme_bw(base_size = 9)
  })
  
  if (length(page_plots) < 6)
    page_plots <- c(page_plots,
                    replicate(6 - length(page_plots), nullGrob(),
                              simplify = FALSE))
  
  grid.arrange(grobs = page_plots, ncol = 3)
}

#Creating the Fourier Basis:
n_basis_fourier <- 23
basis_fourier <- create.fourier.basis(
  rangeval = c(min_time, max_time),
  nbasis   = n_basis_fourier
)

# Evaluate Fourier basis functions at the time points
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




Prices_data_numeric <- as.matrix(Prices_data) 

#Define lambda range
lambda_range <- seq(1e-8, 1e4, length.out = 50)  
Sec_derivative <- int2Lfd(max(0, n_order - 2))


optimal_lambdas <- data.frame(
  Asset = colnames(Prices_data_numeric),
  Optimal_lambda = NA,
  Min_GCV = NA
)
n_assets <- ncol(Prices_data_numeric)
par(mfrow = c(3, 2))
for (asset in 1:n_assets) {
  asset_name <- colnames(Prices_data_numeric)[asset]
  cat(paste("Processing Asset:", asset_name, "\n"))
  
  y <- Prices_data_numeric[, asset]
  
  gcvsave <- numeric(length(lambda_range))
  
  for (i in seq_along(lambda_range)) {
    lambda <- lambda_range[i]
    fdParobj <- fdPar(basis, Sec_derivative, lambda)
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj)
    gcvsave[i] <- smooth_result$gcv  
  }
  
  min_gcv_index <- which.min(gcvsave)
  optimal_lambda <- lambda_range[min_gcv_index]
  min_gcv <- gcvsave[min_gcv_index]
  
  optimal_lambdas[asset, ] <- list(
    Asset = asset_name,
    Optimal_lambda = optimal_lambda,
    Min_GCV = min_gcv
  )
  
  plot(
    lambda_range, gcvsave, type = "b", col = "blue", lwd = 2, pch = 19,
    xlab = "Lambda", ylab = "GCV",
    main = paste("GCV Curve for Asset", asset_name)
  )
  abline(v = optimal_lambda, col = "red", lty = 2, lwd = 1.5)  
}

print(optimal_lambdas)
lowest_gcv_row <- optimal_lambdas[which.min(optimal_lambdas$Min_GCV), ]

cat("Descriptive statistics for Optimal Lambdas (Smoothing Analysis):\n")
print(summary(optimal_lambdas$Optimal_lambda))
cat("Standard Deviation:", sd(optimal_lambdas$Optimal_lambda), "\n\n")

cat("Descriptive statistics for Minimum GCVs (Smoothing Analysis):\n")
print(summary(optimal_lambdas$Min_GCV))
cat("Standard Deviation:", sd(optimal_lambdas$Min_GCV), "\n")

n_time <- length(time_points)  
cat("Total number of time points:", n_time, "\n")


#Smoothing all assets with their respective optimal lambdas
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

combined_fd <- fd(
  coef = do.call(cbind, lapply(smoothed_assets, function(fd) fd$coefs)),
  basisobj = basis
)

#Plotting all smoothed functions
plot(
  combined_fd,
  xlab = "Days", ylab = "Smoothed Functional Value"
)

#Performing Functional Principal Component Analysis
fpca_results <- pca.fd(fdobj = combined_fd, nharm = 2)
fpca_results$harmonics 
fpca_results$scores 
fpca_results$values 
fpca_results$varprop 
fpca_results$meanfd 
plot.pca.fd(fpca_results)


#Plotting the first harmonic
plot(fpca_results$harmonics[1], xlab = "Days", ylab = "PC1 Score", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180))

#Plotting the second harmonic
plot(fpca_results$harmonics[2], xlab = "Days", ylab = "PC2 Score", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180))

#Plotting the FPCA scores in a scatter plot
library(ggplot2)
ggplot(fpc_scores, aes(x = PC1, y = PC2, label = Asset)) +
  geom_point(color = "blue", size = 3) +             
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +    
  labs(title = "FPCA Scores Scatter Plot",
       x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()


