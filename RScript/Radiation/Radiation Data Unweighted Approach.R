install.packages("readxl")
install.packages("dplyr")
install.packages("tidyr")
install.packages("fda")
install.packages("imputeTS")
install.packages("ggplot2")
install.packages("forecast")
install.packages("seastests")
install.packages("tseries")
install.packages("gridExtra")
install.packages("GeneCycle")


library(gridExtra)   
library(GeneCycle)
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

Radiation_Data <- read_excel("C:/Users/jonat/OneDrive/Desktop/Thesis/Data Files/Global Radiation.xlsx")
View(Radiation_Data)

Radiation_Data <- Radiation_Data %>%
  pivot_wider(
    names_from = Station,
    values_from = `Global radiation (in J/cm2)`,
    id_cols = Date
  )

radts <- ts(Radiation_Data[, 23])

plot.ts(
  radts,
  xlab = "Days",
  ylab = "Global radiation (in J/cm2)",
  main = "Time Series of Global Radiation for station 330"
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

stations <- colnames(Radiation_Data_imp)

op <- par(no.readonly = TRUE)           
on.exit(par(op), add = TRUE)            

par(ask = FALSE)                       

for (i in seq(1, length(stations), by = 4)) {
  
  grp <- stations[i : min(i + 3, length(stations))]  
  
  plot.new()                        
  par(mfrow = c(2, 2), mar = c(4.5, 4, 4, 1))
  
  for (st in grp) {
    acf(Radiation_Data_imp[[st]],
        main     = paste("ACF – Station", st),
        xlab     = "Lag",
        ylab     = "ACF",
        lag.max  = 20,
        cex.main = 1.1)
  }
  
  if (length(grp) < 4) for (k in seq_len(4 - length(grp))) plot.new()
}

seasonality_results <- data.frame(
  Station = character(),
  IsSeasonal = logical(),
  stringsAsFactors = FALSE
)

for (station in colnames(Radiation_Data_imp)) {
  ts_data <- ts(Radiation_Data_imp[[station]], frequency = 365)
  
  seasonal_check <- isSeasonal(ts_data)
  
  if (seasonal_check) {
    cat("Seasonality detected for station:", station, "\n")
  } else {
    cat("No significant seasonality detected for station:", station, "\n")
  }
  
  seasonality_results <- rbind(
    seasonality_results,
    data.frame(Station = station, IsSeasonal = seasonal_check, stringsAsFactors = FALSE)
  )
}
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
matplot(time_points_radiation,PHI,type='l',lwd=1,lty=1, xlab='Days',ylab='B-spline function values',cex.lab=1,cex.axis=1)
for (i in 1:n_knots)
{
  abline(v=knots[i], lty=2, lwd=1)}

#Finding the optimal number of basis functions for the fourier basis
K_values <- seq(4, 100, by = 2)  

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
stations <- unique(results_df$Station)   

for (i in seq(1, length(stations), by = 4)) {
  
  grp  <- stations[i : min(i + 3, length(stations))]   
  pg   <- lapply(grp, function(st) {
    ggplot(subset(results_df, Station == st),
           aes(x = K, y = s2)) +
      geom_line() +
      geom_point(size = 1.5) +
      labs(x = "K",
           y = expression(S^2)) +
      theme_bw(base_size = 10)
  })
  
  if (length(pg) < 4) pg <- c(pg, replicate(4 - length(pg),
                                            nullGrob(), simplify = FALSE))
  
  grid.arrange(grobs = pg, ncol = 2)     
}

#Creating the Fourier Basis:
n_basis_fourier <- 7
basis_fourier <- create.fourier.basis(
  rangeval = c(min_time, max_time),
  nbasis   = n_basis_fourier
)

#Evaluate Fourier basis functions at the time points
PHI_fourier <- eval.basis(time_points_radiation, basis_fourier)
dim(PHI_fourier)

mse_df <- data.frame(
  Station         = colnames(Radiation_Data_imp),
  MSE_Bspline   = NA_real_,
  MSE_Fourier   = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(colnames(Radiation_Data_imp))) {
  
  Y <- as.numeric(Radiation_Data_imp[[i]])
  
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

#So in this case, a fourier basis is more optimal - Hinting at the periodicity of the data. 

Radiation_Data_imp_numeric<- as.matrix(Radiation_Data_imp)
n_stations<-ncol(Radiation_Data_imp_numeric)
n_stations
#Define lambda range
#Unweighted Analysis
lambda_range_unweighted <- seq(1e-8, 11000, length.out = 50)  # Lambda values to test
Sec_derivative <- int2Lfd(max(0, n_order - 2))

optimal_lambdas_radiation_unweighted <- data.frame(
  Station = colnames(Radiation_Data_imp_numeric),
  Optimal_lambda_unweighted = NA,
  Min_GCV_unweighted = NA
)

par(mfrow = c(1, 1))

optimal_lambdas_radiation_unweighted <- data.frame(
  Station                 = character(n_stations),
  Optimal_lambda          = numeric(n_stations),
  Min_GCV                 = numeric(n_stations),
  stringsAsFactors        = FALSE
)

par(
  mfrow    = c(3, 2),
  mar      = c(4, 4, 5, 1),
  cex.main = 0.9,
  ask      = FALSE
)

for (i in seq_len(n_stations)) {
  station_name <- colnames(Radiation_Data_imp_numeric)[i]
  cat(sprintf("Processing station %d: %s\n", i, station_name))
  
  y <- Radiation_Data_imp_numeric[, i]
  
  gcv_vals <- numeric(length(lambda_range_unweighted))
  for (j in seq_along(lambda_range_unweighted)) {
    lam     <- lambda_range_unweighted[j]
    fdP     <- fdPar(basis_fourier, Sec_derivative, lam)
    smoothr <- smooth.basis(argvals = time_points_radiation, y = y, fdParobj = fdP)
    gcv_vals[j] <- smoothr$gcv
  }
  
  idx_opt <- which.min(gcv_vals)
  lam_opt <- lambda_range_unweighted[idx_opt]
  gcv_opt <- gcv_vals[idx_opt]
  
  optimal_lambdas_radiation_unweighted[i, ] <- list(
    Station        = station_name,
    Optimal_lambda = lam_opt,
    Min_GCV        = gcv_opt
  )
  
  plot(
    lambda_range_unweighted, gcv_vals,
    type  = "b",      
    lwd   = 2,        
    pch   = 19,       
    xlab  = "Lambda",
    ylab  = "GCV",
    main  = paste(strwrap(station_name, width = 30), collapse = "\n")
  )
  abline(v = lam_opt, col = "red", lty = 2, lwd = 1.5)
}

par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)

print(optimal_lambdas_radiation_unweighted)

print(optimal_lambdas_radiation_unweighted)
optimal_lambda_values <- as.numeric(optimal_lambdas_radiation_unweighted$Optimal_lambda_unweighted)
min_gcv_values <- as.numeric(optimal_lambdas_radiation_unweighted$Min_GCV_unweighted)



fitted_values_matrix_radiation_unweight <- matrix(NA, nrow = n_time, ncol = n_stations)

#Smooth all stations with their respective optimal lambdas
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

fitted_values_df_un <- as.data.frame(fitted_values_matrix_radiation_unweight)
colnames(fitted_values_df_un) <- colnames(Radiation_Data_imp_numeric)

View(fitted_values_df_un)

#Combine smoothed functions into a single functional data object for plotting
combined_fd_radiation <- fd(
  coef = do.call(cbind, lapply(smoothed_stations_unweighted, function(fd) fd$coefs)),
  basisobj = basis_fourier
)

#Plot all smoothed functions
plot(
  combined_fd_radiation,
  xlab = "Days", ylab = "Smoothed Functional Value"
)

#Performing Functional Principal Component Analysis 
unweighted_fpca_results_radiation <- pca.fd(fdobj = combined_fd_radiation, nharm = 3)
unweighted_fpca_results_radiation$harmonics 
unweighted_fpca_results_radiation$scores 
unweighted_fpca_results_radiation$values 
unweighted_fpca_results_radiation$varprop 
unweighted_fpca_results_radiation$meanfd 
plot.pca.fd(unweighted_fpca_results_radiation)
unweighted_fpca_results_radiation$varprop

par(mfrow = c(1, 1))
#Plotting the first harmonic
plot(unweighted_fpca_results_radiation$harmonics[1], xlab = "Days", ylab = "PC1 Score", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

#Plotting the second harmonic
plot(unweighted_fpca_results_radiation$harmonics[2], xlab = "Days", ylab = "PC2 Score", xaxt = 'n')
axis(1, at = seq(min(time_points_radiation), max(time_points_radiation), by = 180), 
)

fpc_scores <- as.data.frame(unweighted_fpca_results_radiation$scores)

colnames(fpc_scores) <- c("PC1", "PC2")

fpc_scores$Station <- colnames(Radiation_Data_imp_numeric)

fpc_scores <- fpc_scores[, c("Station", "PC1", "PC2")]

print(fpc_scores)

View(fpc_scores)

ggplot(fpc_scores, aes(x = PC1, y = PC2, label = Station)) +
  geom_point(color = "blue", size = 3) +         # Plot the points
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) + # Add station labels near points
  labs(x = "PC1 Score",
       y = "PC2 Score") +
  theme_minimal()


combined_fd_radiation_deriv <- deriv.fd(combined_fd_radiation, Lfdobj = 1)

#Plotting the first derivatives for all stations in one plot
plot(combined_fd_radiation_deriv, 
     xlab = "Time (Days)", 
     ylab = "First Derivative", 
     main = "First Derivative of Smoothed Radiation Functions")

