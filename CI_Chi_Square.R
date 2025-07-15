#------------------------------------------------------------
# Simulation study to estimate the empirical coverage probability
# of the chi-square theory confidence interval for the between-unit
# variance component (sigma^2_u) in a one-way random effects model
#
# For each (a, r) design (number of units, number of replicates),
# we generate repeated datasets, estimate variance components,
# construct confidence intervals, and compute coverage probability.

# Note: This CI does not incorporate sigma2_e in confidence interval.
#------------------------------------------------------------

rm(list = ls())
set.seed(2025-07-14)

# Design grid: combinations of (a, r)
ar_matrix <- data.frame(
  # a = c(6, 8, 12, 24, 32, 48),      # Number of units
  # r = c(16, 12, 8, 4, 3, 2)         # Number of replicates per unit
  a = c(6),
  r = c(16)
)

# True variance component
sigma2_u <- 0.5
sigma2_e <- 1

# Confidence interval parameters
alpha <- 0.05

# Simulation parameters
n_designs <- 5e5                               # Number of simulated datasets per design

#------------------------------------------------------------
# Loop over each (a, r) design
#------------------------------------------------------------

for (ar_pair in 1:nrow(ar_matrix)) {
  a <- ar_matrix[ar_pair, "a"]
  r <- ar_matrix[ar_pair, "r"]
  
  cat("Processing a =", a, ", r =", r, "\n")      # Progress message
  
  covered_chi <- 0          # Counter for how many times the true sigma2_u is covered
  sum_length_chi <- 0       # sum of CI lengths 
  count_zero <- 0           # count how frequently boundary estimates occur
  
  # Confidence interval parameters using chi-square approximation
  chi_U <- qchisq(1-alpha/2,a-1,lower.tail = TRUE)
  chi_L <- qchisq(alpha/2,a-1,lower.tail = TRUE)
  
  #------------------------------------------------------------
  # Inner loop: simulate many datasets under the current (a, r) design
  #-----------------------------------------------------------
  for (design in 1:n_designs){
    
    #------------------------------------------------------------
    # Step 1: Generate a dataset with 'a' units and 'r' replications each
    # Each unit has a random effect from N(0, sigma2_u)
    # Each measurement includes an independent error from N(0, sigma2_e)
    #------------------------------------------------------------
    
    data_profile = matrix(nrow = a, ncol =r)
    for (i in 1:a){
      Y <- rnorm(n=1, mean = 0, sd=sqrt(sigma2_u)) + rnorm(n=r, mean = 0, sd=sqrt(sigma2_e))
      data_profile[i,] <- Y 
    }
    
    #------------------------------------------------------------
    # Step 2: Estimate variance components using MLE method
    #------------------------------------------------------------
    
    group_mean <- rowMeans(data_profile)       # Mean for each unit
    overal_mean <- mean(group_mean)            # Overal mean across all data
    
    SS_u <- r*sum((group_mean - overal_mean)^2)    # Sum of squares between units
    SS_e <- sum((data_profile - matrix(rep(group_mean, each = r), nrow = a, byrow = TRUE)) ^ 2)   
                                                   # Sum of squares within units
    SS_t <- SS_u + SS_e                            # Total sum of squares
    
    MS_u <- SS_u/(a-1)
    MS_e <- SS_e/(a*(r-1))
    beta <- a/(a-1)
    
    # MLE estimates (constrained to be within plausible range)
    
    sigma2_u_MLE = pmax(0, (beta^{-1}*MS_u - MS_e)/r)
    sigma2_e_MLE = pmin(SS_t/(a*r), MS_e)
    count_zero <- count_zero + as.numeric(sigma2_u_MLE == 0) 
    
    #------------------------------------------------------------
    # Step 3: Construct chi-square CI and check if true sigma2_u is inside
    #------------------------------------------------------------

    # Check if the CI contains true value
    covered_chi <- covered_chi + as.numeric((a*sigma2_u_MLE/chi_U <= sigma2_u) & (sigma2_u <= a*sigma2_u_MLE/chi_L))
    
    # Track total CI length for this design
    sum_length_chi <- sum_length_chi + (a*sigma2_u_MLE/chi_L - a*sigma2_u_MLE/chi_U)
    
  }
  #------------------------------------------------------------
  # Attach coverage probabilities and average lengths to design matrix and display
  #------------------------------------------------------------
  
  # Empirical coverage probability for this (a, r) design
  ar_matrix$Pr_chi[ar_pair] <- covered_chi/n_designs
  
  # Average length of confidence interval for this (a, r) design
  ar_matrix$mean_length_chi[ar_pair] <- sum_length_chi / n_designs
  
  # Count the frequency of boundry estimates
  ar_matrix$prop_zero[ar_pair] <- count_zero / n_designs
  
}

# View results
print(ar_matrix)

