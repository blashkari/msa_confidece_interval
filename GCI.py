import numpy as np
import pandas as pd
from scipy.stats import norm, chi2

# Set random seed
np.random.seed(2025_07_14)

# GCI function for sigma^2_u
def gci_sigma2_u(MSU, MSE, a, r, n_sim=10000):
    gpq_vals = np.empty(n_sim)
    
    for i in range(n_sim):
        W1 = chi2.rvs(df=a - 1)
        W2 = chi2.rvs(df=a * (r - 1))
        
        MSE_star = MSE * a * (r - 1) / W2
        MSU_star = MSU * (a - 1) / W1
        
        sigma2_u_star = max(0, (MSU_star - MSE_star) / r)
        gpq_vals[i] = sigma2_u_star
        
    return np.quantile(gpq_vals, [0.025, 0.975])

#--------------------------------------------------
# Design grid: combinations of (a, r)
designs = pd.DataFrame({
    #'a': [6],
    #'r': [16]
    'a': [6, 8, 12, 24, 32, 48],     # Number of units
    'r': [16, 12, 8, 4, 3, 2]        # Number of replicates per unit
})
#--------------------------------------------------

# True variance components
sigma2_u = 0.5
sigma2_e = 1

# Confidence level
alpha = 0.05
z_alpha_over_2 = norm.ppf(1 - alpha / 2)

# Simulation parameters
n_designs = 500000

#------------------------------------------------------------
# Loop over each (a, r) design
#------------------------------------------------------------

for idx, row in designs.iterrows():

    a = int(row['a'])
    r = int(row['r'])
    print(f"Processing a = {a}, r = {r}")

    covered_gci = 0
    CI_lengths_gci = np.empty(n_designs)

    #------------------------------------------------------------
    # Inner loop: simulate many datasets under the current (a, r) design
    #-----------------------------------------------------------
    for d in range(n_designs):
        # Generate data
        data = np.zeros((a, r))
        for i in range(a):
            unit_effect = np.random.normal(0, np.sqrt(sigma2_u))
            errors = np.random.normal(0, np.sqrt(sigma2_e), size=r)
            data[i, :] = unit_effect + errors
        
        # Estimate variance components
        group_means = data.mean(axis=1)
        overall_mean = group_means.mean()

        SS_u = r * np.sum((group_means - overall_mean) ** 2)
        SS_e = np.sum((data - group_means[:, np.newaxis]) ** 2)
        SS_t = SS_u + SS_e

        MS_u = SS_u / (a - 1)
        MS_e = SS_e / (a * (r - 1))
        beta = a / (a - 1)

        sigma2_u_MLE = max(0, (MS_u / beta - MS_e) / r)
        sigma2_e_MLE = min(SS_t / (a * r), MS_e)

        # Compute GCI and check coverage
        gci = gci_sigma2_u(MS_u, MS_e, a, r)
        lower, upper = gci

        if lower <= sigma2_u <= upper:
            covered_gci += 1
        
        CI_lengths_gci[d] = upper - lower

    # Store results
    designs.loc[idx, 'Pr_gci'] = covered_gci / n_designs
    designs.loc[idx, 'mean_length_gci'] = np.mean(CI_lengths_gci)
    designs.loc[idx, 'median_length_gci'] = np.median(CI_lengths_gci)

print(designs)