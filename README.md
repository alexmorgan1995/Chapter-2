# Chapter 2 Model - Modelling the effects of livestock antibiotic usage on human food-borne disease

## Foodborne AMR Model

Listed in this GitHub repository are the models (some are deprecated...) for my 2nd PhD chapter. This chapter looks at creating a simple modelling framework to explore the impact of livestock antibiotic usage (τ) on human health, specifically overall levels of foodborne disease and the proportion of antibiotic-resistant foodborne disease in humans.

**As of [03/21] the final version of these analyses can be found in the ```/NewFits_041021``` folder**

## Folder Structure

We currently have two main folders, a ```/Deprecated``` and a ```/NewFits_041021``` folder. The main set of analyses can be found in the latter folder. This folder contains both the raw extracted data for the main analyses, including human and livestock antibiotic resistance and usage datasets.

This was used for the ABC-SMC model fitting procedure in the main ```/models``` folder for four main case-studies. Within this folder there are a number of sub-folders which correspond to either descriptive analyses to explore the data (```/models/Desc_Anal```) or to conduct the modelling analyses in the chapter (```/models/Analyses```). There is also a folder for deprecated model analyses files (```/models/Deprecated```).

## Model Details

All code was written and run in R-Studio. The packages required to run the code are as follows:

```library("deSolve"); library("fast"); library("sensitivity"); library("ggplot2"); library("plotly"); library("reshape2")```

deSolve is the primary package used to solve the ODEs used in the modelling approach.

