Function to assess the statistical power for testing the null-hypothesis of no effect and 
the null-hypothesis of homogeneous true effect size of a Many Labs project. Monte-Carlo 
simulations are used for computing the statistical power.

The function can be used for standardized mean differences and raw correlation coefficients 
as effect size measures. Different values for the intra-class correlation (rho) and true 
effect size (mu) have to be specified together with information on the sample sizes in the 
primary studies.

By setting the argument "report" to TRUE, a HTML report is created. Note that a prerequisite 
for this is that the files "power_ml.R" and "report_power_ml.Rmd" should be save in the same 
directory of your computer. If "report" is set to FALSE, the output is only shown in the 
R console.