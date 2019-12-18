setwd("//csce.datastore.ed.ac.uk/csce/biology/users/s1678248/PhD/Mathematical Models")

rm(list=ls())

library("fitR"); library("MASS"); data(SEITL_stoch); data(FluTdC1971)
plotTraj(data = FluTdC1971)

traj <- FluTdC1971

#The next part is to think of a sumamry statistic that you might be able to use for the model. Examples for this could include...
#1) Maximum Epidemic Size at Peak 
#2) Overall Number of Cases 
#3) Number of Cases at time = 7 
#4) Length of Epidemic Time 
#5) Number of Cases after time 25 

#Summary Statistics 
summarystatmax <- function(traj) {
  return(max(traj$obs, na.rm = T))
}

summarystatlength <- function(traj) { 
  return(min(traj[which(traj$obs == max(traj$obs)),]$time, na.rm = T))
}

summarystattwenty <- function(traj) {
  return(sum(traj$obs[traj$time >= 20], na.rm = T))
}

#Simualting Dataset based on Parameters 

theta <- c(R0 = 2, D_lat = 2, D_inf = 2, alpha = 0.9, D_imm = 13, rho = 0.85)

init.state <- c(S = 250, E = 0, I = 4, T = 0, L = 30, Inc = 0)

#The replicated function allows you to replicate atomic coordinates using boundary conditions - repeat a data.frame or function multiple times?
#So the code detailed below essentially repeats the model run 100 times 

hist(replicate(100, summarystatmax(rTrajObs(SEITL_stoch, theta, init.state, FluTdC1971$time))))
hist(replicate(100, summarystatlength(rTrajObs(SEITL_stoch, theta, init.state, FluTdC1971$time))))

#Time to calculate the difference between the summary statistics identified from data and the model statistics 

simu <- rTrajObs(SEITL_stoch, theta, init.state, FluTdC1971$time)

funclist <- list(summarystatmax, summarystatlength, summarystattwenty)

#Interestingly the x arguement is essentially a function fro the sum stats input

my_distance <- function(sum.stats, data.obs, model.obs) {
  sumdata <- sapply(sum.stats, function(x) {
    abs(x(data.obs) - x(model.obs))
  })
  return(sumdata)
}

my_distance(funclist, 
            FluTdC1971, 
            simu)

#it can take these extremely basic inputs without any need for further subsetting (i.e giving the function just time)
#Because I have previously specified the subsetting in the previous function equations (traj$obs etc.)
#This function allows me to determine the difference between my model data and the actual observed data 

#The next step is to try and determine the distribution of distances - but first it may be beneficial to create one function which can take #
#All of the inputs into the model - rather than taking a dataframe calcualted from a pervious model run - Integrate into one function

computeDistanceABC(sum.stats = list(summarystatmax, summarystatlength, summarystattwenty),
                   distanceABC = my_distance,
                   fitmodel = SEITL_stoch,
                   theta = theta,
                   init.state = init.state,
                   data = FluTdC1971)

#It would next be interesting to see this replicated many times so you can see a distribution of distances 
#You would probably only be able to run this function with a single summary statistic each time - so cut this down

hist(replicate(100, computeDistanceABC(sum.stats = list(summarystatmax),
                                  distanceABC = my_distance,
                                  fitmodel = SEITL_stoch,
                                  theta = theta,
                                  init.state = init.state,
                                  data = FluTdC1971)))

hist(replicate(100, computeDistanceABC(sum.stats = list(summarystatlength),
                                       distanceABC = my_distance,
                                       fitmodel = SEITL_stoch,
                                       theta = theta,
                                       init.state = init.state,
                                       data = FluTdC1971)))

hist(replicate(100, computeDistanceABC(sum.stats = list(summarystattwenty),
                                       distanceABC = my_distance,
                                       fitmodel = SEITL_stoch,
                                       theta = theta,
                                       init.state = init.state,
                                       data = FluTdC1971)))

#We next need a function that returns the posterior distribution for the theta parameter using the ABC rejection algorithm - 1 generation
#We can use the gamma distribution - this means that the parameter values we use for the model need to be positive 
#We can estimate two parameters Dlat and Dinf - DInf Duration of Infection, Dinf - Duration of latent period 
#The first thing to do is choose informative priors for my parameter -  we choose the shape as 16 and the rate as 8 - from this we calculate the mean and variance
#This is nice as we actually have a distribution which is informative with a gamma distribution 

# a function that takes 7 arguments:
# - N :the number of desired samples
# - epsilon: a vector (if the distance function returns a vector) or
#            single number, the tolerance for ABC
# - sum.stats: list of summary statistic functions
# - distanceABC: ABC distance function
# - fitmodel: model to compare to the data
# - init.state: initial state for the simulation
# - data: data to compare the model to

ABC_algorithm <- function(N, epsilon, sum.stats, distanceABC, fitmodel, init.state, data) {
  dump1 <- matrix(nrow = 0 , ncol = 6)
  i <- 0
  while(i < N) {
    dlat <- rgamma(1, shape = 16, rate = 8)
    dinf <- rgamma(1, shape = 16, rate = 8)
    theta <- c(R0 = 2, D_lat = dlat, D_inf = dinf, alpha = 0.9, D_imm = 13, rho = 0.85)
    dist <- computeDistanceABC(sum.stats = sum.stats,
                       distanceABC = distanceABC,
                       fitmodel = fitmodel,
                       theta = theta,
                       init.state = init.state,
                       data = data)
    if(dist[1] <= epsilon[1] | dist[2] <= epsilon[2]) {
      dump1 <- rbind(dump1, theta)
    }
    i <- dim(dump1)[1]
    print(i)
  }
  return(dump1)
}

postdist <- data.frame(ABC_algorithm(N = 1000, epsilon = c(5,5), sum.stats = list(summarystatlength, summarystatmax),
              distanceABC = my_distance, fitmodel = SEITL_stoch, init.state = init.state, data = FluTdC1971))

hist(postdist$D_lat)
hist(postdist$D_inf)
hist(rgamma(1000, shape = 16, rate = 8))

#### ABC-SMC #### 

# use the ABC rejection algorithm to find population 1 in the ABC-SMC algorithm
ABC_algorithm <- function(N, epsilon, sum.stats, distanceABC, fitmodel, init.state, data) {
  dump1 <- matrix(nrow = 0 , ncol = 6)
  i <- 0
  while(i < N) {
    dlat <- rgamma(1, shape = 16, rate = 8)
    dinf <- rgamma(1, shape = 16, rate = 8)
    theta <- c(R0 = 2, D_lat = dlat, D_inf = dinf, alpha = 0.9, D_imm = 13, rho = 0.85)
    dist <- computeDistanceABC(sum.stats = sum.stats,
                               distanceABC = distanceABC,
                               fitmodel = fitmodel,
                               theta = theta,
                               init.state = init.state,
                               data = data)
    if(dist[1] <= epsilon[1]) {
      dump1 <- rbind(dump1, theta)
    }
    i <- dim(dump1)[1]
    print(i)
  }
  return(dump1)
}

gen1 <- data.frame(ABC_algorithm(N = 1000, epsilon = 5, sum.stats = list(summarystatlength),
                                     distanceABC = my_distance, fitmodel = SEITL_stoch, init.state = init.state, data = FluTdC1971))

ABC_SMC_algorithm <- function(N, epsilon_2, sum.stats, distanceABC, fitmodel, init.state, data) {
  sigma <- matrix(c(0.5,0,0.5,0), 2, 2)
  results <- matrix(nrow = 0, ncol = 6)
  i <- 0 
  while(i < N) {
    row_no <- sample(1000, 1)
    sample1 <- as.numeric(gen1[row_no,c("D_lat","D_inf")]) # sample used in the mvrnorm function as mu - which is the vector giving the means of the parameters
    #of the variables - which is d_lat and d_inf
    gaus_dist <- mvrnorm(1, mu = sample1, sigma) # with this we already have a parameter which has been peturbed by the gaussian peturbatio kernel
    theta <- c(R0 = 2, D_lat = gaus_dist[[1]], D_inf = gaus_dist[[2]], alpha = 0.9, D_imm = 13, rho = 0.85)
    print(c(theta[2], theta[3]))
    dist <- computeDistanceABC(sum.stats = sum.stats,
                             distanceABC = distanceABC,
                             fitmodel = fitmodel,
                             theta = theta,
                             init.state = init.state,
                             data = data)
    if(dist <= epsilon_2) {
      results <- rbind(results, theta)
      }
    i <- dim(results)
  }
  return(results)
}

gen2 <- data.frame(ABC_SMC_algorithm(N = 1000, epsilon_2 = 5, sum.stats = list(summarystatlength),
                                 distanceABC = my_distance, fitmodel = SEITL_stoch, init.state = init.state, data = FluTdC1971))

# Basic framework is here - but the model must account for the increases in "generation" not just particles 
# I also need to calculate weights for every new generation after the first one
# I also need to do the vector for the thresholds (epsilon which should come in a vector most likely - linked to the no. of generations)

sigma <- matrix(c(0.5,0,0.5,0), 2, 2)
test1 <- as.numeric(gen1[sample(1000, 1), c("D_lat", "D_inf")])
gaus_dist <- mvrnorm(2, mu = test1, sigma)
theta <- c(R0 = 2, D_lat = gaus_dist[[1]], D_inf = gaus_dist[[2]], alpha = 0.9, D_imm = 13, rho = 0.85)




