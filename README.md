# Chapter 2 Model - Modelling the effects of livestock antibiotic usage on human food-borne disease

## Foodborne AMR Model

Listed in this GitHub repository are the models (some are deprecated...) for my 2nd PhD chapter. This chapter looks at creating a simple modelling framework to explore the impact of livestock antibiotic usage (τ) on human health, specifically overall levels of foodborne disease and the proportion of antibiotic-resistant foodborne disease in humans.

As of [07/11/19] the final script for the analysis is named: ```Chapter 2 - Finalised Model - Campylobacter in Poultry```

This script is currently being altered to make the model analysis I carried out more streamlined. Especially with regards to testing different parameter combinations and plotting these parameter sets. 

## Model Details

All code was written and run in R-Studio. The packages required to run the code are as follows:

```library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("tidyr")```

deSolve is the primary package used to solve the ODEs used in the modelling approach 


This platform will provide a basic model in which further model iterations can be based off. 
