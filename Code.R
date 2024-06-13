---
title: "Thesis_Code"
author: "Ivo Bonfanti"
date: "2024-02-05"
---

  
Chunk 0: Libraries
```{r}  
library(zoo)
library(ggplot2)
library(tidyr)
library(bvartools)
library(LaplacesDemon)
library(dplyr)
library(textshape)
library(gridExtra)
setwd("XXXXXX")
```

Chunk 1: Functions
```{r}
# Extract submatrices (row)-------------------------------------------------
extract_submatrix <- function(i, matrix) {
  pattern <- paste0("0", i)
  #all row names including "0i"
  rows_to_extract <- grep(pattern, rownames(matrix))
  return(matrix[rows_to_extract, , drop = FALSE])
}

# Extract submatrices (col)-------------------------------------------------
extract_submatrix1 <- function(i, matrix) {
  pattern <- paste0("0", i)
  #all col names including "0i"
  col_to_extract <- grep(pattern, colnames(matrix))
  return(matrix[, col_to_extract, drop = FALSE])
}

#compute the mean per lag---------------------------------------------------
lagmean <- function(vect, lag){
  vect = as.double(vect)
  meanlist <- c()
  lim = (length(vect)/lag) -1
  for(i in 0: lim){
    mean = 0
    pos = i*lag
    for(l in 1: lag){
      mean = mean + vect[pos+l]*(1 / 2^(l))
    }
    meanlist[i+1] <- mean
  }
  return(meanlist)
}

# Function to plot a random row from a df ---------------------------------
plot_random_row <- function(df) {
  # Select a random row index
  random_index <- sample(nrow(df), 1)
  # Select the random row
  random_row <- df[random_index, ]
  # Create the plot
  ggplot(data = data.frame(Value = as.vector(t(random_row))), 
         aes(x = 1:length(random_row), y = Value)) +
         geom_line() +
         labs(title = paste("Coefficient:", random_index),
              x = "Column Index",
              y = "Value") +
         theme_minimal()
}


# Signal adaptive variable selector function ------------------------------
SAVS <- function (matrix1, matrix2) {
  # Initialize an empty matrix for the sparse coefficients
  sparse_matrix<-list()
  for (i in 1:lag) {
    subcoeff<-extract_submatrix(i, matrix1)
    subexog<-extract_submatrix1(i, matrix2)
    r<-nrow(subcoeff)
    c<-ncol(subcoeff)
    # Iterate along the subcoeff entries
    for (s in 1:r) {
      for (j in 1:c) {
        beta<-subcoeff[s, j]
        sql2norm<-t(subexog[,s])%*%subexog[,s]
        # If own lags (diagonal)
        if (s==j) {
          mu<-(i-1)^2*(beta)^(-2)
          mu <- ifelse(is.nan(mu), 0, mu)
          if (abs(beta)*sql2norm <= mu) {
            subcoeff[s,j]<-0
          }else {
            subcoeff[s,j]<-sign(beta)*(abs(beta)*sql2norm - mu)/sql2norm
          }
        } else {
          mu<-i^2*beta^(-2)
          if (abs(beta)*sql2norm <= mu) {
            subcoeff[s,j]<-0
          } else {
            subcoeff[s,j]<-sign(beta)*(abs(beta)*sql2norm - mu)/sql2norm
          }
        }
      }
    }
    sparse_matrix[[i]] <- subcoeff
    # Sparse matrix of coefficients
    final_matrix <- do.call(rbind, sparse_matrix)
  }
  return(final_matrix)
}
```

Chunk 2: Preprocessing and plots
```{r}
# Upload the data
series<-read.csv("sovereign.csv")
series<-read.csv("banking.csv")

# Linear interpolation and transformation --------------------------------
# Check missing values
print(na_count <- colSums(is.na(series))) 
# Liner interpolation
series[, -1]<-lapply(series[, -1], na.approx)
print(colSums(is.na(series)))

# From absolute to % to basis points
series[, -1]<-series[, -1]*100*100
# Compute log first differences and multiply by 100
log_diff<-apply(series[, -1], 2, function(x) diff(log(x))*100)

# Plot sovereign series ---------------------------------------------------
# Greece excluded for a matter of visualization
sov_plot <- series[, -which(names(series) == "Hellenic.Rep")]

# Converts to long format
sov_plot <- pivot_longer(sov_plot, cols = -Date, names_to = "Countries", 
                         values_to = "Spread") 

sov_plot$Date <- as.Date(sov_plot$Date, format = "%Y-%m-%d") 

sov_graph<- ggplot(sov_plot, aes(x = Date, y = Spread, color = Countries)) +
            geom_line() +
            labs(x = "Date", y = "Spread", color = "Countries") +
            ggtitle("Sovereign CDS") +theme_minimal() +
            geom_line(linewidth = 0.5) +
            theme(
              axis.title.x = element_text(family = "Arial", size = 12),
              axis.title.y = element_text(family = "Arial", size = 12),
              axis.text.x = element_text(angle = 45, family = "Arial"),
              axis.text.y = element_text(family = "Arial")) +
            scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m-%d") +
            ylim(0, 400)+
            guides(color = "none") #remove legend

ggsave("sov_graph.png", width = 15, height = 7.5, dpi=700, bg="white")


# Plot banking series -----------------------------------------------------
# MPS excluded for a matter of visualization
ban_plot <- series[, -which(names(series) ==
                              "Bca.Monte.dei.Paschi.di.Siena.S.p.A")]

ban_plot <- pivot_longer(ban_plot, cols = -Date, names_to = "Banks", 
                         values_to = "Spread")

ban_plot$Date <- as.Date(ban_plot$Date, format = "%Y-%m-%d") 

ban_graph<- ggplot(ban_plot, aes(x = Date, y = Spread, color = Banks)) +
            geom_line() +
            labs(x = "Date", y = "Spread", color = "Banks") +
            ggtitle("Banking CDS") +theme_minimal() +
            geom_line(linewidth = 0.5) +
            theme(
              axis.title.x = element_text(family = "Arial", size = 12),
              axis.title.y = element_text(family = "Arial", size = 12),
              axis.text.x = element_text(angle = 45, family = "Arial"),
              axis.text.y = element_text(family = "Arial")) +
            scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m-%d") +
            ylim(0, 300)+
            guides(color = "none")

ggsave("ban_graph.png", width = 15, height = 7.5, dpi=700, bg="white")


# Plot sovereign first differences ----------------------------------------
# Add dates back
sov_first <- cbind(series[-1, 1, drop = FALSE], log_diff) 

sov_plot1 <- pivot_longer(sov_first, cols = -Date, names_to = "Countries", 
                          values_to = "LogDiff") 

sov_plot1$Date <- as.Date(sov_plot1$Date, format = "%Y-%m-%d") 

sov_graph1<- ggplot(sov_plot1, aes(x = Date, y = LogDiff, color = Countries)) +
             geom_line() +
             labs(x = "Date", y = "LogDiff", color = "Countries") +
             ggtitle("Sovereign CDS log-differences") +theme_minimal() +
             geom_line(linewidth = 0.5) +
             theme(
              axis.title.x = element_text(family = "Arial", size = 12),
              axis.title.y = element_text(family = "Arial", size = 12),
              axis.text.x = element_text(angle = 45, family = "Arial"),
              axis.text.y = element_text(family = "Arial")) +
             scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m-%d") +
             ylim(-50,50)+
             guides(color = "none") 

ggsave("sov_graph1.png", width = 15, height = 7.5, dpi=700, bg="white")

# Plot banking first differences ------------------------------------------
# Add dates back
ban_first <- cbind(series[-1, 1, drop = FALSE], log_diff) 

ban_plot1 <- pivot_longer(ban_first, cols = -Date, names_to = "Banks", 
                          values_to = "LogDiff") 

ban_plot1$Date <- as.Date(ban_plot1$Date, format = "%Y-%m-%d") 

ban_graph1<- ggplot(ban_plot1, aes(x = Date, y = LogDiff, color = Banks)) +
             geom_line() +
             labs(x = "Date", y = "LogDiff", color = "Banks") +
             ggtitle("Banking CDS log-differences") +theme_minimal() +
             geom_line(linewidth = 0.5) +
             theme(
              axis.title.x = element_text(family = "Arial", size = 12),
              axis.title.y = element_text(family = "Arial", size = 12),
              axis.text.x = element_text(angle = 45, family = "Arial"),
              axis.text.y = element_text(family = "Arial")) +
             scale_x_date(date_breaks = "1 year", date_labels = "%Y-%m-%d") +
             ylim(-1, 1)+
             guides(color = "none")

ggsave("ban_graph1.png", width = 15, height = 7.5, dpi=700, bg="white")
```

Chunk 3: Full BVAR to be run on Colab
```{r}
# Model specification ------------------------------------------------------
set.seed(19)

# Data 
data<-as.ts(log_diff)

#Define the number of lags
lag=3 

# Generate the model
model <- gen_var(data, p=lag)

# Extract the data
# Endogenous variables
y <- t(model$data$Y)
# Exogenous variables
x <- t(model$data$Z)

# Minnesota prior defined
Litterman_prior <- minnesota_prior(
  model,
  kappa0 = 0.04,
  kappa1 = 0.25,
  kappa2 = NULL,
  kappa3 = 5,
  max_var = NULL,
  coint_var = TRUE,
  sigma = "VAR"
)

# Gibbs sampler specification ---------------------------------------------
# Number of iterations of the Gibbs sampler
draws <- 10000 
# Number of burn-in draws
burnin <- 1000
store = draws-burnin

# Number of observations 
tt <- ncol(y)
# Number of endogenous variables
k <- nrow(y) 
# Number of estimated coefficients for each VAR
m <- k * nrow(x)

# Litterman prior parameters 
# Vector of prior parameter means
beta_mu_prior <- Litterman_prior$mu 
# Inverse of the prior covariance matrix
beta_v_i_prior <- Litterman_prior$v_i 

# Inverse Wishart parameters
# Degrees of freedom
u_nu_prior <- k+1 
# Scale matrix
zeta_scale_prior<- diag(1, k)
# Posterior degrees of freedom
u_nu_post <- tt + u_nu_prior
# Sigma initialization (ones could take 10*I or use OLS estimator instead)
u_sigma <- diag(1, k) 

# Gibbs sampler -----------------------------------------------------------
# Data containers for posterior draws
# For each iteration m coefficients
draws_a <- matrix(NA, m, store)
# For each iteration k^2 coefficients --> VAR-COV matrices
draws_sigma <- matrix(NA, k*k, store)

# Define progress bar
pb <- txtProgressBar(min = 0, max = draws, style = 3, width = 50, char = "=")

# Start Gibbs sampler
for (draw in 1:draws) {
  setTxtProgressBar(pb, draw) #progress bar
  
  # Draw conditional mean parameters
  beta <- post_normal(y, x, u_sigma, beta_mu_prior, beta_v_i_prior)
  
  # Draw variance-covariance matrix
  u <- y - matrix(beta, k) %*% x # Obtain residuals
  u_zeta_post <-zeta_scale_prior + tcrossprod(u) # Obtain posterior scale matrix
  u_sigma <-rinvwishart(u_nu_post, u_zeta_post) # Estimate Sigma with iW
  
  # Store Beta and Sigma
  if (draw > burnin) {
    draws_a[, draw - burnin] <- beta 
    draws_sigma[, draw - burnin] <- u_sigma 
  }
}

# Gibbs convergence check -------------------------------------------------
# Upload stored and burn-in coefficients computed through Colab
setwd("XXXXXX")
draws_a<-read.csv("model7.csv")
draws_a1000<-read.csv("sovdrawsa1000.csv")

# Merge the two dfs and select one column every 20
gibbs <- cbind(draws_a, draws_a1000)
index <- seq(1, ncol(gibbs), by = 20)
gibbs <- gibbs[,index]

# Create a list to store plots
plots_list <- lapply(1:6, function(i) plot_random_row(gibbs))

# Arrange plots in a grid
grid.arrange(grobs = plots_list, ncol = 3)

# Coefficient matrix ------------------------------------------------------
# Obtain means for every row
A <- rowMeans(draws_a)
# Transform mean vector into a matrix
A <- matrix(A, k) 
# Round values
A <- round(A, 3) 
# Rename matrix dimensions
dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]])
#Convert the matrix to a dataframe
A <- data.frame(A)
#Order the columns alphabetically
order = sort(colnames(A))
A <- A[, order]
#remove the constant
A <- A %>% select(-c('const'))

# Signal adaptive variable selector ---------------------------------------
# Prepare the design matrix
exog<-t(x[-nrow(x), ])
# Coefficients matrix
coeff<-t(A)

# Get the order of rows and columns in alphabetical order
sorted_row_indices <- order(rownames(coeff))
sorted_col_indices <- order(colnames(coeff))
# Rearrange the matrix alphabetically
coeff <- coeff[sorted_row_indices, sorted_col_indices]

#Do the same with the exogenous variables matrix
sorted_col_indices <- order(colnames(exog))
exog <-exog[, sorted_col_indices]

# SAVS function to obtain a sparse matrix of coefficients
final_matrix<-SAVS(coeff, exog)

# Create the adjacency matrix ---------------------------------------------
# Reorder the matrix alphabetically
sorted_row <- order(rownames(final_matrix))
final_matrix<- final_matrix[sorted_row, ]
final_matrix<-t(final_matrix)

# Create a dataframe
col_names <- row.names(final_matrix)
df = data.frame(matrix(nrow = 0, ncol = (length(col_names)) ))
colnames(df) <- col_names

# Derive the weighted mean for the lags of the coefficients
for(j in 1: k){
  lmean = c()
  row = final_matrix[j, ]
  lmean = lagmean(row, lag)
  df[nrow(df)+1, ] <- lmean
} 

# Adjust column names
df['col_names'] <- col_names
row.names(df) <- df$col_names
df <- subset(df, select = -col_names)

# From df to matrix
ad_mat <- as.matrix(df)

# set to zero the elements of the main diagonal 
diag(ad_mat) <- 0

# Retrieve the mean of the absolute values of the matrix non zero entries
threshold = mean(abs(ad_mat[ad_mat !=0])) 

# Set to 0 entries below the threshold and round the results
ad_mat[abs(ad_mat) < threshold] <- 0
ad_mat=round(ad_mat,3)

# Output directory
output_directory <- "XXXXXX"

# Filename for the CSV file
output_filename <- "sovadj.csv"
#output_filename <- "banadj.csv"

# Full path for the output file
output_path <- file.path(output_directory, output_filename)

# Write the data frame to the CSV file
write.csv(ad_mat, file = output_path, row.names = FALSE)
```

Chunk 4: BVAR with rolling window to be run on Colab
```{r}
# Rolling window set up ---------------------------------------------------
# 13 windows (tot 1961 observations)
# About 3 years
windows_size<-791
# About 3 months
step_size<-90
num_obs <- nrow(log_diff)
num_vars <- ncol(log_diff)
num_models<-(num_obs-windows_size)/step_size

# Compute a BVAR for each rolling window ----------------------------------
set.seed(19)
lag=1 

# Loop to run the model for each time period
for (i in 0:(num_models)) {
  # Model and Gibbs specification -----------------------------------------
  # Subset the dataset
  data<-as.ts(log_diff[((i*step_size)+1):((i*step_size)+windows_size),])
  # Retrieve the subset starting and ending dates
  dates<-series[c((i*step_size)+1, (i*step_size)+windows_size), 1]
  
  # Generate the model
  model <- gen_var(data, p=lag)
  # Endogenous variables
  y <- t(model$data$Y)
  # Exogenous variables
  x <- t(model$data$Z)
  
  # Minnesota prior defined
  Litterman_prior <- minnesota_prior(
    model,
    kappa0 = 0.04,
    kappa1 = 0.25,
    kappa2 = NULL,
    kappa3 = 5,
    max_var = NULL,
    coint_var = TRUE,
    sigma = "VAR"
  )
  # Litterman prior parameters 
  # Vector of prior parameter means
  beta_mu_prior <- Litterman_prior$mu 
  # Inverse of the prior covariance matrix
  beta_v_i_prior <- Litterman_prior$v_i 
  
  # Number of observations 
  tt <- ncol(y)
  # Number of endogenous variables
  k <- nrow(y) 
  # Number of estimated coefficients for each VAR
  m <- k * nrow(x)
  
  # Inverse Wishart parameters
  # Degrees of freedom
  u_nu_prior <- k+1 
  # Scale matrix
  zeta_scale_prior<- diag(1, k)
  # Posterior degrees of freedom
  u_nu_post <- tt + u_nu_prior
  # Sigma initialization 
  u_sigma <- diag(1, k) 
  
  # Gibbs sampler ---------------------------------------------------------
  # Data containers for posterior draws
  # For each iteration m coefficients
  draws_a <- matrix(NA, m, store)
  # For each iteration k^2 coefficients --> VAR-COV matrices
  draws_sigma <- matrix(NA, k*k, store)
  # Define progress bar
  pb <- txtProgressBar(min = 0, max = draws, style = 3, width = 50, char = "=")
  # Start Gibbs sampler
  for (draw in 1:draws) {
    setTxtProgressBar(pb, draw) #progress bar
    
    # Draw conditional mean parameters
    beta <- post_normal(y, x, u_sigma, beta_mu_prior, beta_v_i_prior)
    
    # Draw variance-covariance matrix
    u <- y - matrix(beta, k) %*% x # Obtain residuals
    u_zeta_post <-zeta_scale_prior + tcrossprod(u) # Obtain posterior scale matrix
    u_sigma <-rinvwishart(u_nu_post, u_zeta_post) # Estimate Sigma with iW
    
    # Store Beta and Sigma
    if (draw > burnin) {
      draws_a[, draw - burnin] <- beta 
      draws_sigma[, draw - burnin] <- u_sigma 
    }
  }
  
  # Coefficient matrix ------------------------------------------------------
  # Obtain means for every row
  A <- rowMeans(draws_a)
  # Transform mean vector into a matrix
  A <- matrix(A, k) 
  # Round values
  A <- round(A, 3) 
  # Rename matrix dimensions
  dimnames(A) <- list(dimnames(y)[[1]], dimnames(x)[[1]])
  #Convert the matrix to a dataframe
  A <- data.frame(A)
  #Order the columns alphabetically
  order = sort(colnames(A))
  A <- A[, order]
  #remove the constant
  A <- A %>% select(-c('const'))
  
  # Signal adaptive variable selector ---------------------------------------
  # Prepare the design matrix
  exog<-t(x[-nrow(x), ])
  # Coefficients matrix
  coeff<-t(A)
  
  # Get the order of rows and columns in alphabetical order
  sorted_row_indices <- order(rownames(coeff))
  sorted_col_indices <- order(colnames(coeff))
  # Rearrange the matrix alphabetically
  coeff <- coeff[sorted_row_indices, sorted_col_indices]
  
  #Do the same with the exogenous variables matrix
  sorted_col_indices <- order(colnames(exog))
  exog <-exog[, sorted_col_indices]
  
  # SAVS function to obtain a sparse matrix of coefficients
  final_matrix<-SAVS(coeff, exog)
  
  # Create the adjacency matrix ---------------------------------------------
  # Reorder the matrix alphabetically
  sorted_row <- order(rownames(final_matrix))
  final_matrix<- final_matrix[sorted_row, ]
  final_matrix<-t(final_matrix)
  
  # Create a dataframe
  col_names <- row.names(final_matrix)
  df = data.frame(matrix(nrow = 0, ncol = (length(col_names)) ))
  colnames(df) <- col_names
  
  # Derive the weighted mean for the lags of the coefficients
  for(j in 1: k){
    lmean = c()
    row = final_matrix[j, ]
    lmean = lagmean(row, lag)
    df[nrow(df)+1, ] <- lmean
  } 
  
  # Adjust column names
  df['col_names'] <- col_names
  row.names(df) <- df$col_names
  df <- subset(df, select = -col_names)
  
  # From df to matrix
  ad_mat <- as.matrix(df)
  
  # set to zero the elements of the main diagonal 
  diag(ad_mat) <- 0
  
  # Retrieve the mean of the absolute values of the matrix non zero entries
  threshold = mean(abs(ad_mat[ad_mat !=0])) 
  
  # Set to 0 entries below the threshold and round the results
  ad_mat[abs(ad_mat) < threshold] <- 0
  ad_mat=round(ad_mat,3)
  
  # Output directory
  output_directory <- "XXXXXX"
  #output_directory <- "XXXXXX"
  
  # Filename for the CSV file
  output_filename <- paste(dates[1], dates[2], sep="_")
  
  # Full path for the output file
  output_path <- file.path(output_directory, output_filename)
  
  # Write the data frame to the CSV file
  write.csv(ad_mat, file = output_path, row.names = FALSE)
}
```