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

Prices$Date <- as.Date(Prices$Date)
time_points <- as.numeric(Prices$Date - min(Prices$Date)) 
time_points
length(time_points)
max(time_points)

Prices_data <- Prices %>% dplyr::select(-Date)
View(Prices_data)



# Checking for Periodicity
set.seed(16)  
sample_size <- 6  
sampled_assets <- sample(colnames(Prices_data), sample_size)

par(mfrow = c(ceiling(sample_size / 2), 2))  

#ACF Plots
for (asset in sampled_assets) {
  main_title <- paste("ACF of", asset, "Prices")
  
  acf(Prices_data[[asset]], main = main_title, lag.max = 50)
}

par(mfrow = c(1, 1))



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
  abline(v=knots[i], lty=2, lwd=1)
}

# GARCH Method
library(tseries)

#ADF Test
adf_results <- data.frame(
  Asset = colnames(Prices_data),
  ADF_p_value = NA,
  ADF_Test_Statistic = NA
)

differenced_data <- list()

for (asset in colnames(Prices_data)) {
  cat(paste("Performing ADF test for asset:", asset, "\n"))
  
  asset_data <- Prices_data[[asset]]
  
  adf_test <- adf.test(asset_data, alternative = "stationary")
  
  adf_results[adf_results$Asset == asset, ] <- data.frame(
    Asset = asset,
    ADF_p_value = adf_test$p.value,
    ADF_Test_Statistic = adf_test$statistic
  )
  
  differenced_data[[asset]] <- diff(asset_data)
}
adf_results
differenced_data <- do.call(cbind, differenced_data)

differenced_data <- as.data.frame(differenced_data)
View(differenced_data)
print(adf_results)

head(differenced_data)


stationary_assets <- adf_results %>% filter(ADF_p_value < 0.05)
stationary_assets
adf_results
View(differenced_data)

#Dataframe for adf on differenced data
adf_diff_results <- data.frame(
  Asset = colnames(differenced_data), 
  ADF_p_value = NA,
  ADF_Test_Statistic = NA
)

for (asset in colnames(differenced_data)) { 
  cat(paste("Performing ADF test on differenced data for asset:", asset, "\n"))
  
  diff_data <- differenced_data[[asset]]
  
  adf_test <- adf.test(diff_data, alternative = "stationary")
  
  adf_diff_results[adf_diff_results$Asset == asset, ] <- data.frame(
    Asset = asset,
    ADF_p_value = adf_test$p.value,
    ADF_Test_Statistic = adf_test$statistic
  )
}
adf_diff_results
#All p-values smaller than 0.05 so the data is no longer non-stationary so now we fit ARMA on the differenced data. 

#Fitting the ARMA model on each of the differenced assets
library(forecast)

arma_models <- list()
arma_residuals <- list()

#Fitting ARMA model for each differenced column (asset)
for (asset in colnames(differenced_data)) {  
  cat(paste("Fitting ARMA model for asset:", asset, "\n"))
  
  diff_data <- differenced_data[[asset]]
  
  arma_model <- auto.arima(diff_data, d = 0, seasonal = FALSE)
  
  arma_models[[asset]] <- arma_model
  arma_residuals[[asset]] <- residuals(arma_model)
  
  print(summary(arma_model))
}

arma_model_summaries <- data.frame(
 AR_Order = sapply(arma_models, function(model) model$arma[1]),
  MA_Order = sapply(arma_models, function(model) model$arma[2]),
  AIC = sapply(arma_models, AIC),
  BIC = sapply(arma_models, BIC)
)
arma_model_summaries

#Obtaining the residuals in the form of a matrix
n_assets <- length(arma_models)
n_assets
n_time <- nrow(differenced_data)
n_time

residuals_matrix <- matrix(NA, nrow = n_time, ncol = n_assets,
                           dimnames = list(NULL, names(arma_models)))

#Looping through each fitted ARMA model to extract residuals
for (asset in names(arma_models)) {
  cat(paste("Extracting residuals for asset:", asset, "\n"))
  
  residuals <- arma_residuals[[asset]]
  
  residuals_matrix[1:length(residuals), asset] <- residuals
}

residuals_df <- as.data.frame(residuals_matrix)
View(residuals_df)

#Testing for ARCH effects
library(tseries)

#Looping through each asset's residuals
garch_effects_results <- data.frame(
  Asset = names(arma_models),
  Ljung_Box_p_value = NA,
  ARCH_LM_p_value = NA
)

for (asset in names(arma_models)) {
  cat(paste("Testing for GARCH effects on asset:", asset, "\n"))
  
  residuals <- arma_residuals[[asset]]
  
  squared_residuals <- residuals^2
  
  #Performing Ljung-Box test on squared residuals
  ljung_box_test <- Box.test(squared_residuals, lag = 8, type = "Ljung-Box")
  
  #Performing ARCH LM test
  arch_lm_test <- ArchTest(residuals, lags = 8)
  
  garch_effects_results[garch_effects_results$Asset == asset, ] <- list(
    Asset = asset,
    Ljung_Box_p_value = ljung_box_test$p.value,
    ARCH_LM_p_value = arch_lm_test$p.value
  )
}
garch_effects_results
non_significant_garch <- garch_effects_results %>% 
  filter(Ljung_Box_p_value > 0.05 | ARCH_LM_p_value > 0.05)

significant_garch <- garch_effects_results %>%
  filter(Ljung_Box_p_value <= 0.05 & ARCH_LM_p_value <= 0.05)

print(non_significant_garch)

garch_models <- list()

#Looping through each asset and fit GARCH model for those with significant effects
for (asset in significant_garch$Asset) {
  cat(paste("Fitting GARCH(1,1) model for asset:", asset, "\n"))
  
  residuals <- arma_residuals[[asset]]
  
  #Defining the GARCH(1,1) specification
  garch_spec <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"  #Assuming normal distribution for residuals
  )
  
  garch_fit <- ugarchfit(spec = garch_spec, data = residuals)
  
  garch_models[[asset]] <- garch_fit
  
  print(garch_fit)
}

#Extracting the conditional variances for each asset
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

#Assigning equal variance (1) for assets with non-significant GARCH effects
for (asset in non_significant_garch$Asset) {
  cat(paste("Assigning equal weights for asset with no significant GARCH effects:", asset, "\n"))
  garch_variances[[asset]] <- rep(1, length(arma_residuals[[asset]]))
}

garch_variances_df <- as.data.frame(garch_variances)
View(garch_variances_df)

#Obtaining the reciprocal of the conditional variances for weighting
reciprocal_variances <- list()

for (asset in names(garch_variances)) {
  cat(paste("Calculating reciprocal of conditional variances for asset:", asset, "\n"))
  
  conditional_variances <- garch_variances[[asset]]
  
  reciprocal_variances[[asset]] <- ifelse(conditional_variances > 0, 
                                          1 / conditional_variances, 
                                          NA)
}

reciprocal_variances_df <- as.data.frame(reciprocal_variances)
View(reciprocal_variances_df)

new_row <- reciprocal_variances_df[nrow(reciprocal_variances_df), ]  
reciprocal_variances_df <- rbind(reciprocal_variances_df, new_row) 
View(reciprocal_variances_df)

Prices_data_numeric <- as.matrix(Prices_data)
View(Prices_data_numeric)

# List of assets to update
assets_to_update <- c("NZYM-B.CO", "ALIV-SDB.ST", "ASSA-B.ST", "SKF-B.ST", "VOLV-B.ST")

colnames(Prices_data_numeric) <- sapply(colnames(Prices_data_numeric), function(col) {
  if (col %in% assets_to_update) {
    sub("-", ".", col)  
  } else {
    col  
  }
})

#Obtaining the GCV across lambda values for the assets using the GARCH approach
Sec_derivative <- int2Lfd(max(0, n_order - 2))
n_assets <- ncol(reciprocal_variances_df)  
n_assets
lambda_range <- seq(1e-8, 1e4, length.out = 50)  

preliminary_results <- list()  
optimal_lambdas <- data.frame(Asset = colnames(reciprocal_variances_df), 
                              Optimal_lambda_GCV = NA, 
                              Min_GCV = NA)

par(mfrow = c(3, 2))
for (asset in colnames(reciprocal_variances_df)) {
  cat(paste("Calculating metrics for asset:", asset, "\n"))
  
  #Extracting the data and weight vector for the current asset
  y <- Prices_data_numeric[, asset]
  wtvec <- reciprocal_variances_df[[asset]]
  
  gcvsave <- numeric(length(lambda_range))
  
  for (i in seq_along(lambda_range)) {
    lambda <- lambda_range[i]  
    fdParobj <- fdPar(basis, Sec_derivative, lambda = lambda)  
    
    smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj, wtvec = wtvec)
    
    gcvsave[i] <- smooth_result$gcv  
  }
  
  min_gcv_index <- which.min(gcvsave)
  optimal_lambda_gcv <- lambda_range[min_gcv_index]
  min_gcv <- gcvsave[min_gcv_index]
  
  optimal_lambdas[optimal_lambdas$Asset == asset, ] <- list(
    Asset = asset,
    Optimal_lambda_GCV = optimal_lambda_gcv,
    Min_GCV = min_gcv
  )
  
  preliminary_results[[asset]] <- list(
    lambda_range = lambda_range,
    gcvsave = gcvsave
  )
  
  plot(lambda_range, gcvsave, type = "b", col = "blue", lwd = 2, pch = 19,
       xlab = "Lambda", ylab = "GCV", 
       main = paste("GCV Curve for Asset", asset))
  abline(v = optimal_lambda_gcv, col = "red", lty = 2, lwd = 1.5)  
  
}

print(optimal_lambdas)

#Optimal lambda statistics
cat("Descriptive statistics for Optimal Lambdas (Weighted Analysis):\n")
print(summary(optimal_lambdas$Optimal_lambda_GCV))
cat("Standard Deviation:", sd(optimal_lambdas$Optimal_lambda_GCV), "\n\n")

#Minimum GCV statistics
cat("Descriptive statistics for Minimum GCVs (Weighted Analysis):\n")
print(summary(optimal_lambdas$Min_GCV))
cat("Standard Deviation:", sd(optimal_lambdas$Min_GCV), "\n\n")

n_assets <- nrow(optimal_lambdas)
cat("Total number of assets:", n_assets, "\n")


#Smoothing all assets with optimal lambdas
smoothed_assets <- list()

for (asset in optimal_lambdas$Asset) {
  cat(paste("Smoothing asset:", asset, "\n"))
  
  optimal_lambda <- optimal_lambdas[optimal_lambdas$Asset == asset, "Optimal_lambda_GCV"]
  
  y <- Prices_data_numeric[, asset]
  wtvec <- reciprocal_variances_df[[asset]]
  
  fdParobj <- fdPar(basis, Sec_derivative, lambda = optimal_lambda)
  
  smooth_result <- smooth.basis(argvals = time_points, y = y, fdParobj = fdParobj, wtvec = wtvec)
  
  smoothed_assets[[asset]] <- smooth_result$fd
}
smoothed_assets
length(smoothed_assets)

#Combining smoothed functions into a single functional data object
combined_fd_garch <- fd(
  coef = do.call(cbind, lapply(smoothed_assets, function(fd) fd$coefs)), 
  basisobj = basis
)
combined_fd_garch$coefs[,1]
#Plotting the functions
plot(combined_fd_garch, xlab = "Days", ylab = "Value", main = "Smoothed Functions")
evaluated_values <- eval.fd(time_points, combined_fd_garch)
View(evaluated_values)
is.fd(combined_fd_garch)
#Conducting FPCA
fpca_results_garch <- pca.fd(fdobj = combined_fd_garch, nharm = 2)
fpca_results_garch$harmonics 
fpca_results_garch$scores 
fpca_results_garch$values 
fpca_results_garch$varprop 
fpca_results_garch$meanfd 
plot.pca.fd(fpca_results_garch)

#Plotting the harmonics individually (functional principal components)
#Plotting the first harmonic
plot(fpca_results_garch$harmonics[1], main = "First Harmonic (PC1)", xlab="Days",ylab = "PC1", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180))

#Plotting the second harmonic
plot(fpca_results_garch$harmonics[2], main = "Second Harmonic (PC2)",xlab="Days", ylab = "PC2", xaxt = 'n')
axis(1, at = seq(min(time_points), max(time_points), by = 180))
#Converting FPCA scores to a data frame
fpc_scores_garch <- as.data.frame(fpca_results_garch$scores)

colnames(fpc_scores_garch) <- c("PC1", "PC2")

fpc_scores_garch$Asset <- colnames(Prices_data_numeric)

fpc_scores_garch <- fpc_scores_garch[, c("Asset", "PC1", "PC2")]

print(fpc_scores_garch)

View(fpc_scores_garch)
library(ggplot2)
ggplot(fpc_scores_garch, aes(x = PC1, y = PC2, label = Asset)) +
  geom_point(color = "blue", size = 3) +             
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +    
  labs(title = "FPCA Scores Scatter Plot",
       x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()























