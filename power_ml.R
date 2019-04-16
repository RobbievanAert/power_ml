################################################################################
##### STATISTICAL POWER FOR TESTING NULL HYPOTHESIS OF NO EFFECT AND NO    #####
##### BETWEEN-STUDY VARIANCE IN A META-ANALYSIS                            #####
##### Author: Robbie C.M. van Aert                                         #####
################################################################################

rm(list = ls())

#################
### FUNCTIONS ###
#################

### Function for computing estimate of tau2 with Paule-Mandel estimator
PM <- function(tau2, yi, vi) 
{
  df <- length(yi) - 1 # Degrees of freedom of chi-square distribution
  wi <- 1/(vi+tau2) # Weights in meta-analysis
  theta <- sum(yi*wi)/sum(wi) # Meta-analytic effect size
  Q <- sum(wi*(yi-theta)^2) # Q-statistic
  Q - df # Q-statistic minus degrees of freedom
}

### Function to get power in meta-analysis for SMD and correlation
get_power <- function(es, n1i, n2i, ni, rhos, mus, alpha_mu = 0.05, tail = "two", 
                      alpha_tau2 = 0.05, reps = 10000, tau2_max = 5, report = TRUE)
{
  # es = effect size measure used (SMD or correlation)
  # n1i = sample size group 1 (SMD)
  # n2i = sample size group 2 (SMD)
  # ni = sample size (correlation)
  # rhos = intra-class corelations
  # mus = true effect sizes
  # alpha_mu = alpha-level for testing null hypothesis of no effect
  # tail = whether one or two-tailed tests is conducted in primary studies
  # alpha_tau2 = alpha-level for Q_test
  # reps = number of replications simulation results are based upon
  # tau2_max = upper bound for root-finding algorithm
  # report = whether you want to get a HTML report of the results
  
  if (es == "SMD")
  { # Standardized mean difference (Hedges' g)
    k <- length(n1i) # Number of studies in meta-analysis
    df <- n1i+n2i-2 # Degrees of freedom
    
  } else if (es == "COR")
  {
    k <- length(ni) # Number of studies in meta-analysis
  }
  
  ### Empty object for storing results
  pow_mu <- pow_tau2 <- mean_I2 <- matrix(NA, nrow = length(mus), ncol = length(rhos), 
                                          dimnames = list(as.character(mus), as.character(rhos)))
  
  ### Create progress bar
  pb <- txtProgressBar(min = 0, max = length(rhos)*length(mus), style = 3)
  m <- 0
  Sys.sleep(0.1)
  
  for (rho in rhos)
  {
    for (mu in mus)
    {
      
      null_sim <- Q_sim <- I2 <- numeric(reps) # Empty objects for storing results    
      
      tau2 <- -rho/(rho-1) # Compute tau^2 based on ICC and sigma2 = 1
      
      for (i in 1:reps)
      {
        if (es == "SMD")
        { # Standardized mean difference (Hedges' g)
          mdiffi <- rnorm(k, mean = mu, sd = sqrt(1/n1i+1/n2i+tau2)) # Observed mean difference
          s21i <- 1/(n1i-1) * rchisq(k, df = n1i-1)         # Observed variance group 1
          s22i <- 1/(n2i-1) * rchisq(k, df = n2i-1)         # Observed variance group 2
          pool <- ((n1i-1)*s21i + (n2i-1)*s22i)/(n1i+n2i-2) # Observed pooled variances of mean difference
          di <- mdiffi/sqrt(pool) # Cohen's d
          J <- exp(lgamma(df/2) - log(sqrt(df/2)) - lgamma((df-1)/2)) # Hedges' g correction factor
          yi <- J * di # Compute Hedges' g
          ### Unbiased estimator of variance (Viechtbauer, 2007, equation 23)
          vi <- 1/n1i+1/n2i+(1-(n1i+n2i-2-2)/((n1i+n2i-2)*J^2))*yi^2 
          
        } else if (es == "COR")
        { # Correlation coefficient (after Fisher-z transformation)
          yi <- rnorm(k, mean = mu, sd = sqrt(1/(ni-3)+tau2))
          vi <- 1/(ni-3)
        }
        
        ### Estimate tau^2 with Paule-Mandel estimator
        if (PM(tau2 = 0, yi = yi, vi = vi) < 0)
        { # If estimate is smaller than zero, set it equal to zero
          tau2_PM <- 0
        } else
        { # Estimate tau2 with PM if estimate is larger than zero
          tau2_PM <- uniroot(PM, interval = c(0, tau2_max), yi = yi, vi = vi)$root
        }
        
        wi_star <- 1/(vi+tau2_PM) # Weights RE model
        est <- sum(wi_star*yi)/sum(wi_star) # Estimate RE model
        se_est <- sqrt(1/sum(wi_star)) # SE of estimate RE model
        
        wi <- 1/vi # Weights EE model
        s2 <- (k-1)*sum(wi)/(sum(wi)^2-sum(wi^2)) # Typical within-study variance
        I2[i] <- tau2_PM/(s2+tau2_PM)*100 # I2-statistic
        
        #########################
        ### TEST OF NO EFFECT ###
        #########################
        
        if (tail == "two")
        { # Compute two-tailed p-value
          pval <- ifelse(est > 0, 2*pnorm(est/se_est, lower.tail = FALSE), 2*pnorm(est/se_est))
        } else if (tail == "one")
        { # Compute one-tailed p-value
          pval <- pnorm(est/se_est, lower.tail = FALSE)
        }
        
        null_sim[i] <- pval < alpha_mu # Check whether p-value is smaller than alpha_mu
        
        ############################################
        ### TEST OF HOMOGENEOUS TRUE EFFECT SIZE ###
        ############################################
        
        est <- sum(wi*yi)/sum(wi) # Estimate EE model
        Qstat <- sum(wi*(yi-est)^2) # Q-statistic
        pval_Q <- pchisq(Qstat, df = k-1, lower.tail = FALSE) # P-value Q-statistic
        Q_sim[i] <- pval_Q < alpha_tau2 # Check whether p-value is smaller than alpha_tau2
        
      }
      
      ### Compute statitical power across replications
      pow_mu[as.character(mu),as.character(rho)] <- mean(null_sim)
      pow_tau2[as.character(mu),as.character(rho)] <- mean(Q_sim)
      
      ### Mean I2-statistic across replications
      mean_I2[as.character(mu),as.character(rho)] <- mean(I2)
      
      ### Update progress bar
      m <- m + 1
      setTxtProgressBar(pb, m)
      
    }
  }
  
  close(pb) # Close progress bar
  
  if (report == TRUE)
  { # If the user want to see the report
    res <- list(pow_mu = pow_mu, pow_tau2 = pow_tau2, mean_I2 = mean_I2)
    save(res, file = "res.RData") # Save results to working directory
    rmarkdown::render("report_power_ml.Rmd") # Create report
    browseURL(file.path("report_power_ml.html")) # Open report
  }
  
  return(list(pow_mu = pow_mu, pow_tau2 = pow_tau2, mean_I2 = mean_I2))
}

################################################################################
################################################################################
################################################################################

### THE USER HAS TO SPECIFY THE FOLLOWING INFORMATION FOR APPLYING THE FUNCTION:
# es = effect size measure used --> standardized mean difference ("SMD") or correlation ("COR")
# n1i = vector of sample sizes group 1 (only for SMD)
# n2i = vector of sample sizes group 2 (only for SMD)
# ni = vector of sample sizes (only for COR)
# rhos = vector of intra-class corelations
# mus = vector of true effect sizes
# alpha_mu = alpha-level for testing null hypothesis of no effect (default = 0.05)
# tail = whether null-hypothesis of no effect is tested one- ("one") or two-tailed 
  # ("two") (default = "two")
# alpha_tau2 = alpha-level for Q_test (default = 0.05)
# reps = number of replications for simulations (default = 10000)
# tau2_max = upper bound for root-finding algorithm for estimating tau2 (default = 5)
# report = whether you want to get a HTML report of the results (in order to create 
  # the report two files will be saved to the working directory of your computer, 
  # default = TRUE)

### Example standardized mean difference
rhos <- c(0, 0.1, 0.25) # Intra-class correlations
mus <- c(0, 0.5, 1) # True effect size
n1i <- n2i <- c(15, 20, 30, 40, 50, 60, 70, 80, 90, 100) # Sample sizes group 1 and 2

get_power(es = "SMD", n1i = n1i, n2i = n2i, rhos = rhos, mus = mus)

### Example correlation
rhos <- c(0, 0.1, 0.25) # Intra-class correlations
mus <- c(0, 0.5, 1) # True effect size
ni <- c(15, 20, 30, 40, 50, 60, 70, 80, 90, 100) # Sample sizes 

get_power(es = "COR", ni = ni, rhos = rhos, mus = mus)