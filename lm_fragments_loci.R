## load data from clipboard
sim.data <- read.table(file = "clipboard", header=TRUE)
emp.data <- read.table(file = "clipboard", header=TRUE)

## reduce data to working samples
sim.red <- subset(emp.data, taxon != "idahoensis" & final_loci_all > 5000)

## linear regression models on fragments
fit.sim.gf0 <- lm(div_from_mela ~ genome_framents_0, data = sim.data)
summary(fit.sim.gf0)

# Call:
# lm(formula = div_from_mela ~ genome_framents_0, data = sim.data)

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -2.6437 -0.5449  0.1226  1.0603  1.8852 

# Coefficients:
                    # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)       13.5444974  3.5681853   3.796 0.000694 ***
# genome_framents_0 -0.0005166  0.0001602  -3.225 0.003117 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.407 on 29 degrees of freedom
# Multiple R-squared:  0.2639,	Adjusted R-squared:  0.2385 
# F-statistic:  10.4 on 1 and 29 DF,  p-value: 0.003117

##linear regression on fragments without novomexicana
fit.sim.gf0.no31 <- lm(div_from_mela ~ genome_framents_0, data = sim.data, subset = (1:length(div_from_mela) != 31))
summary(fit.sim.gf0.no31)

# Call:
# lm(formula = div_from_mela ~ genome_framents_0, data = sim.data, 
    # subset = (1:length(div_from_mela) != 31))

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -2.5597 -0.6539  0.1411  1.1379  1.8729 

# Coefficients:
                    # Estimate Std. Error t value Pr(>|t|)   
# (Intercept)       12.5310324  3.8826359   3.227  0.00318 **
# genome_framents_0 -0.0004726  0.0001736  -2.723  0.01102 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.42 on 28 degrees of freedom
# Multiple R-squared:  0.2093,	Adjusted R-squared:  0.1811 
# F-statistic: 7.413 on 1 and 28 DF,  p-value: 0.01102

##linear regression model on loci
fit.sim.fl <- lm(div_from_mela ~ final_loci_all, data = sim.data)
summary(fit.sim.fl)

# Call:
# lm(formula = div_from_mela ~ final_loci_all, data = sim.data)

# Residuals:
     # Min       1Q   Median       3Q      Max 
# -2.44886 -0.65461  0.05341  1.15409  1.79491 

# Coefficients:
                 # Estimate Std. Error t value Pr(>|t|)    
# (Intercept)     7.055e+00  1.624e+00   4.344 0.000156 ***
# final_loci_all -1.870e-04  6.016e-05  -3.109 0.004179 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.421 on 29 degrees of freedom
# Multiple R-squared:   0.25,	Adjusted R-squared:  0.2242 
# F-statistic: 9.668 on 1 and 29 DF,  p-value: 0.004179

##linear regression on loci without novomexicana
fit.sim.fl.no31 <- lm(div_from_mela ~ final_loci_all, data = sim.data, subset = (1:length(div_from_mela) != 31))
summary(fit.sim.fl.no31)

# Call:
# lm(formula = div_from_mela ~ final_loci_all, data = sim.data, 
    # subset = (1:length(div_from_mela) != 31))

# Residuals:
    # Min      1Q  Median      3Q     Max 
# -2.4423 -0.6680  0.1896  1.1616  1.7927 

# Coefficients:
                 # Estimate Std. Error t value Pr(>|t|)   
# (Intercept)     7.001e+00  2.038e+00   3.436  0.00186 **
# final_loci_all -1.851e-04  7.454e-05  -2.483  0.01927 * 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 1.446 on 28 degrees of freedom
# Multiple R-squared:  0.1805,	Adjusted R-squared:  0.1512 
# F-statistic: 6.167 on 1 and 28 DF,  p-value: 0.01927

##linear regression on loci of empirical data
> fit.emp.fl <- lm(div_from_mela ~ final_loci_all, data = emp.red)
> summary(fit.emp.fl)

# Call:
# lm(formula = div_from_mela ~ final_loci_all, data = emp.red)

# Residuals:
   # Min     1Q Median     3Q    Max 
# -2.304 -1.378 -1.074  1.880  3.171 

# Coefficients:
                # Estimate Std. Error t value Pr(>|t|)
# (Intercept)    5.939e-01  1.428e+00   0.416    0.684
# final_loci_all 7.656e-05  1.034e-04   0.740    0.471

# Residual standard error: 1.951 on 14 degrees of freedom
# Multiple R-squared:  0.03768,	Adjusted R-squared:  -0.03105 
# F-statistic: 0.5483 on 1 and 14 DF,  p-value: 0.4713


