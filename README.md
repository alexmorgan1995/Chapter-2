# Chapter 2 Model - Modelling the effects of livestock antibiotic usage on human food-borne disease

## Foodborne AMR Model

Listed in this GitHub repository are the models (some are deprecated...) for my 2nd PhD chapter. This chapter looks at creating a simple modelling framework to explore the impact of livestock antibiotic usage (τ) on human health, specifically overall levels of foodborne disease and the proportion of antibiotic-resistant foodborne disease in humans.

**As of [07/11/19] the final script for the analysis is named: ```Chapter 2 - Finalised Model - Campylobacter in Poultry```**

This script is currently being altered to make the model analysis I carried out more streamlined. Especially with regards to testing different parameter combinations and plotting these parameter sets. 

## Model Details

All code was written and run in R-Studio. The packages required to run the code are as follows:

```library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2")```

deSolve is the primary package used to solve the ODEs used in the modelling approach. A function was used in this model to help tidy up the prevalence values obtained from numerically solving the set of ODEs - ```rounding```. This is defined as:

```
rounding <- function(x) {
  if(as.numeric(x) < 1e-10) {x <- 0
  } else{signif(as.numeric(x), digits = 6)}
}
```

This was set at a lower limit of 1x10^-10 before the prevalence is calcualted as 0. This is to prevent any rounding errors with the integrator from negatively altering the model output. 

There are currently two types of model structure used in these scripts. 
- With scaling parameters (λ - scaling parameter for antibiotic-mediated recovery) and (ζ - scaling parameter for fitness effects of antibiotic-resistance on transmission)
- Without scaling parameters 

### Code for Model w/ Scaling Parameters 

```
amr <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSa = ua + ra*(Ia + Ira) + zeta*tau*Ia - (betaAA*Ia*Sa) - (betaAH*Ih*Sa) - lambda*(betaAH*Irh*Sa) - lambda*(betaAA*Ira*Sa) - ua*Sa  
    dIa = betaAA*Ia*Sa + betaAH*Ih*Sa + phi*Ira - zeta*tau*Ia - tau*theta*Ia - ra*Ia - ua*Ia
    dIra = lambda*betaAH*Irh*Sa + lambda*betaAA*Ira*Sa + tau*theta*Ia - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ih+Irh) - (betaHH*Ih*Sh) - lambda*(betaHH*Irh*Sh) - (betaHA*Ia*Sh) - lambda*(betaHA*Ira*Sh) - uh*Sh 
    dIh = betaHH*Ih*Sh + betaHA*Ia*Sh - rh*Ih - uh*Ih 
    dIrh = lambda*(betaHH*Irh*Sh) + lambda*(betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIa,dIra,dSh,dIh,dIrh)))
  })
}
```

### Code For Model w/o Scaling Parameters 

```
amrold <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dSa = ua + ra*(Ia + Ira) + tau*Ia - (betaAA*Ia*Sa) - (betaAH*Ih*Sa) - (betaAH*Irh*Sa) - (betaAA*Ira*Sa) - ua*Sa  
    dIa = betaAA*Ia*Sa + betaAH*Ih*Sa + phi*Ira - tau*Ia - tau*theta*Ia - ra*Ia - ua*Ia
    dIra = betaAH*Irh*Sa + betaAA*Ira*Sa + tau*theta*Ia - phi*Ira - ra*Ira - ua*Ira
    
    dSh = uh + rh*(Ih+Irh) - (betaHH*Ih*Sh) - (betaHH*Irh*Sh) - (betaHA*Ia*Sh) - (betaHA*Ira*Sh) - uh*Sh 
    dIh = betaHH*Ih*Sh + betaHA*Ia*Sh - rh*Ih - uh*Ih 
    dIrh = (betaHH*Irh*Sh) + (betaHA*Ira*Sh) - rh*Irh - uh*Irh 
    return(list(c(dSa,dIa,dIra,dSh,dIh,dIrh)))
  })
}
```
