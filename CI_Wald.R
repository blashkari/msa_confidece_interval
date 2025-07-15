#------------------------------------------------------------
# Simulation study to estimate the empirical coverage probability
# of the normal-theory confidence interval for the between-unit
# variance component (sigma^2_u) in a one-way random effects model
#
# For each (a, r) design (number of units, number of replicates),
# we generate repeated datasets, estimate variance components,
# construct confidence intervals, and compute coverage probability.
#------------------------------------------------------------

rm(list = ls())
set.seed(2025-07-14)

# Design grid: combinations of (a, r)
ar_matrix <- data.frame(
   a = c(6, 8, 12, 24, 32, 48),      # Number of units
   r = c(16, 12, 8, 4, 3, 2)         # Number of replicates per unit
  # a = c(6),
  # r = c(16)
)

# True variance component
sigma2_u <- 0.5
sigma2_e <- 1


# Confidence interval parameters
alpha <- 0.05
z_alpha_over_2 <- qnorm(1-alpha/2)     # Critical value from standard normal


# Simulation parameters
n_designs <- 5e5                       # Number of simulated datasets per design

#------------------------------------------------------------
# Loop over each (a, r) design
#------------------------------------------------------------

for (ar_pair in 1:nrow(ar_matrix)) {
  a <- ar_matrix[ar_pair, "a"]
  r <- ar_matrix[ar_pair, "r"]

  cat("Processing a =", a, ", r =", r, "\n")      # Progress message
  
  covered_wald <- 0          # Counter for how many times the true sigma2_u is covered
  sum_length_norm <- 0        # sum of CI lengths 
  count_zero <- 0             # count how frequently boundary estimates occur
  
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
    # Step 3: Construct normal-theory CI and check if true sigma2_u is inside
    #------------------------------------------------------------
    
    var_sigma2_u_MLE <- 2*((sigma2_u_MLE + sigma2_e_MLE/r)^2 + sigma2_e_MLE^2/(r^2*(r-1)))
    
    # Compute z-score for distance from true value
    z_score <- abs(sigma2_u_MLE - sigma2_u)/sqrt(var_sigma2_u_MLE/a)
    
    # If the z-score is less than the critical value, the CI contains true value
    covered_wald <- covered_wald + as.numeric(z_score < z_alpha_over_2)
    
    # Track total CI length for this design
    sum_length_norm <- sum_length_norm + (
      sigma2_u_MLE + z_alpha_over_2 * sqrt(var_sigma2_u_MLE / a) -
        pmax(0, sigma2_u_MLE - z_alpha_over_2 * sqrt(var_sigma2_u_MLE / a))
    )
    
  }
  #------------------------------------------------------------
  # Attach coverage probabilities and average lengths to design matrix and display
  #------------------------------------------------------------
  
  # Empirical coverage probability for this (a, r) design
  ar_matrix$Pr_wald[ar_pair] <- covered_wald/n_designs
  
  # Average length of confidence interval for this (a, r) design
  ar_matrix$mean_length_wald <- sum_length_norm / n_designs
  
  # Count the frequency of boundry estimates
  ar_matrix$prop_zero[ar_pair] <- count_zero / n_designs
  
}


# View results
print(ar_matrix)

