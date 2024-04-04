# OptimInterim
The source code is implemented for OptimInterim, which can determine the optimal timing for interim analyses in group sequential trials. The preprint is available at .

## Get Optimal Information Rates
### Description 
Return the optimal information rates (minimising the expected sample size under the alternative hypothesis) for testing means or rates in one or two samples.

### Usage 
```{r}
getOptimalInformationRates(
    ...,
    groups = 2L,
    allocationRatioPlanned = NA_real_,
    kMax = NA_integer_,
    alpha = NA_real_,
    beta = NA_real_,
    sided = 1L,
    futilityBounds = NA_real_,
    typeOfDesign = c("OF", "P", "WT", "PT", "HP", "WToptimum", "asP", "asOF", "asKD", "asHSD", "asUser", "noEarlyEfficacy"),
    deltaWT = NA_real_,
    deltaPT1 = NA_real_,
    deltaPT0 = NA_real_,
    optimizationCriterion = c("ASNH1", "ASNIFH1", "ASNsum"),
    gammaA = NA_real_,
    typeBetaSpending = c("none", "bsP", "bsOF", "bsKD", "bsHSD", "bsUser"),
    userAlphaSpending = NA_real_,
    userBetaSpending = NA_real_,
    gammaB = NA_real_,
    bindingFutility = NA,
    betaAdjustment = NA,
    constantBoundsHP = 3,
    twoSidedPower = NA,
    delayedInformation = NA_real_,
    tolerance = 1e-08,
    seed = NA_real_
) 
```

### Arguments 
All arguments can be found in the [getDesignGroupSequential](https://rdrr.io/cran/rpact/man/getDesignGroupSequential.html) and [getSampleSizeMeans](https://rdrr.io/cran/rpact/man/getSampleSizeMeans.html) (or [getSampleSizeRates](https://rdrr.io/cran/rpact/man/getSampleSizeRates.html)) functions from the [rpact](https://rdrr.io/cran/rpact) package.

### Details 
The getOptimalInformationRates function

### Value 
Return the optimal information rates.

### Author(s)
The implementation of OptimInterim was written by Zhangyi He working for Laurent Billot and Suzie Cro at the George Institute for Global Health and the Imperial Clinical Trial Unit, Imperial College London.

### References
He, Z., Cro, S., & Billot, L. (2024). Optimal timing for interim analyses in group sequential trials. 

### Examples
```{r}
# Calculate optimal information rates in a four-stage group sequential design 
# with an equal allocation and the O'Brien & Fleming type alpha- and beta-spending functions;
# one-sided significance alpha = 0.025, power 1 - beta = 90%:
getOptimalInformationRatesRates(
  groups = 2,
  allocationRatioPlanned = 1,
  kMax = 4,
  alpha = 0.025,
  beta = 0.1,
  sided = 1,
  typeOfDesign = "asOF",
  typeBetaSpending = "bsOF")
```
See more examples from [OptimInterim_Eg](https://github.com/zhangyi-he/GSD_OptimInterim/blob/main/OptimInterim_Eg.R)).
