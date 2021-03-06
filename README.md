Function to assess the statistical power for testing the null-hypothesis of no effect and 
the null-hypothesis of homogeneous true effect size of a Many Labs project. Monte-Carlo 
simulations are used for computing the statistical power.

The function can be used for standardized mean differences and raw correlation coefficients 
as effect size measures. Different values for the intra-class correlation (rho) and true 
effect size (mu) have to be specified together with information on the sample sizes in the 
primary studies.

By setting the argument "report" to TRUE, a HTML report is created. Note that a prerequisite 
for this is that the files "power_ml.R" and "report_power_ml.Rmd" are saved in the working 
directory of the R session. If "report" is set to FALSE, the output is only shown in the 
R console.

The function has the following arguments:
- es = effect size measure used --> standardized mean difference ("SMD") or correlation ("COR")
- n1i = vector of sample sizes group 1 (only for SMD)
- n2i = vector of sample sizes group 2 (only for SMD)
- ni = vector of sample sizes (only for COR)
- rhos = vector of intra-class corelations
- mus = vector of true effect sizes
- alpha_mu = alpha-level for testing null hypothesis of no effect (default = 0.05)
- tail = whether null-hypothesis of no effect is tested one- ("one") or two-tailed ("two") (default = "two")
- alpha_tau2 = alpha-level for Q_test (default = 0.05)
- reps = number of replications for simulations (default = 10000)
- tau2_max = upper bound for root-finding algorithm for estimating tau2 (default = 5)
- report = whether you want to get a HTML report of the results (in order to create the report two files will be saved to the working directory of your computer, default = TRUE)

``` r
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
```

If you suspect a bug or have a question, please report a report [here](https://github.com/RobbievanAert/power_ml/issues).