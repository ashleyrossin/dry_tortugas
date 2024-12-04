#we're going to run a bayes model and we are not going to cry

##lets bring some data in
library(readr)
library(tidyverse)
library(patchwork)
library(Hmisc)
library(lme4)
library(MCMCpack)
library(jagsUI)
library(dplyr)
MOTE2_measures <- read_csv("MOTE2_measures .csv")
View(MOTE2_measures)

##ok this step is VERY important or you get a bunch of NAs, definitely something to remember

MOTE2_measures$value <- as.numeric(MOTE2_measures$value)
MOTE2_measures$measurement <- as.factor(MOTE2_measures$measurement)

##get one value each
df1Tall <- MOTE2_measures %>% group_by(individual, species, sample_health_state, measurement) %>%
  summarise(max = max(value), min = min(value), mean = mean(value), std = sd(value))

##making it wider and usable for the bayes model
df4 <- df1Tall %>% pivot_wider(id_cols = c(individual, sample_health_state, species), 
                               names_from = measurement, values_from = c(max, min, mean, std)) 

df4$ID <- paste(df4$individual, df4$sample_health_state, sep = "_")

View(df4)



# Group condition
## lets also remove any possible NAs before we start
## we keep health state private from us as we run these quantifications so they get a random letter normally
##in this case, healthy was "J" and diseased was "M" which we learn the day work in R starts
df4$state <- ifelse(df4$sample_health_state %in% c("J"), 0, 1)
View(df4)
df5 <- df4[,c(1,2, 3, 4, 11, 22, 23, 26, 37)]
na.omit(df5)
df6 <- df5 %>% drop_na()
complete.cases(df6)


# Data prep for JAGS model
#centering becoming standardizing, smaller range of values so theyre centered on zero and won't have huge outlierts
#this could also be done by the function scale: as.numeric(scale(lenght)) but this also centers columns and rows
#would no longer be its own vector
##zooxarea has been the name of vacuolization since the onset of this study and while it makes no sense, we keep it for nostalgia
# Exopercent = x1
df6$x1 <- (df6$mean_prop_exo - mean(df6$mean_prop_exo))/sd(df6$mean_prop_exo)

# Gastrosep = x2
df6$x2 <- (df6$mean_gastro - mean(df6$mean_gastro))/sd(df6$mean_gastro)

# Zooxarea = x3
df6$x3 <- (df6$max_symb_vac - mean(df6$max_symb_vac))/sd(df6$max_symb_vac)

# symbiont size = x4
df6$x4 <- (df6$max_avg_symb - mean(df6$max_avg_symb))/sd(df6$max_avg_symb)

# degraded symbionts = x5
df6$x5 <- (df6$mean_degr_symb - mean(df6$mean_degr_symb))/sd(df6$mean_degr_symb)


##what a great time to see if things are correlated
#values from -1 to 1, 0 means no linear relationship
cor(df6$mean_prop_exo, df6$max_symb_vac, use='complete.obs') #-0.056
cor(df6$mean_prop_exo, df6$mean_gastro, use='complete.obs') #0.237
cor(df6$mean_prop_exo, df6$max_avg_symb, use='complete.obs') #-0.0418
cor(df6$mean_prop_exo, df6$mean_degr_symb, use='complete.obs') #NA
cor(df6$mean_gastro, df6$max_symb_vac, use='complete.obs') #-0.0866
cor(df6$mean_gastro, df6$max_avg_symb, use='complete.obs') #0.011
cor(df6$mean_gastro, df6$mean_degr_symb, use='complete.obs') #NA
cor(df6$max_symb_vac, df6$max_avg_symb, use='complete.obs') #0.339
cor(df6$max_symb_vac, df6$mean_degr_symb, use='complete.obs') #NA
cor(df6$mean_degr_symb, df6$max_avg_symb, use='complete.obs') #NA

##nothing is correlated!!!


##now we get into making things for the bayesian model itself
# Enumerate species for indexing as random effect
df6$sp <- as.numeric(as.factor(df6$species))

# Number of species
J <- length(unique(df6$species))

# Number of parameters
##this is one more than params because of intercept
K <- 6

# Create identity matrix for Wishart dist'n
# (Number of parameters to estimate (K))
W <- diag(K)

# Load data
data <- list(y = df6$state, 
             group = df6$sp, 
             n = dim(df6)[1], 
             J = J,
             x1 = df6$x1,
             x2 = df6$x2,
             x3 = df6$x3, 
             x4 = df6$x4,
             x5 = df6$x5,
             K = K, 
             W = W)

## model time!
# Define the model in the BUGS language and write a text file
##
sink("model.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    y[i] ~ dbin(p[i],1)				# distributional assumption
    p[i] <- exp(lp[i])/(1 + exp(lp[i])) # logit link function
    lp[i] <- alpha[group[i]] + beta1[group[i]] * x1[i] + beta2[group[i]] * x2[i] + beta3[group[i]] * x3[i] + beta4[group[i]] * x4[i]	+ beta5[group[i]] * x5[i]# linear predictor    
    } 
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta1[j] <- BB[j,2] # Exo
    beta2[j] <- BB[j,3] # Gastro
    beta3[j] <- BB[j,4] # Zoox
    beta4[j] <- BB[j,5] # symb size
    beta5[j] <- BB[j,6] # degraded symb
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) # bivariate normal
    
    BB.hat[j,1] <- mu.alpha
	  BB.hat[j,2] <- mu.beta1
	  BB.hat[j,3] <- mu.beta2
	  BB.hat[j,4] <- mu.beta3
	  BB.hat[j,5] <- mu.beta4
	  BB.hat[j,6] <- mu.beta5
	  # Use below for 2nd-level of model
    #BB.hat[j,1] <- mu.alpha + gamma0b * z1[j] + gamma0b2 * z2[j] + gamma0b3 * z3[j] 
    #BB.hat[j,2] <- mu.beta + gamma1b * z1[j] + gamma1b2 * z2[j] + gamma1b3 * z3[j] 
    }
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta1 ~ dnorm(0, 0.0001)
    mu.beta2 ~ dnorm(0, 0.0001)
    mu.beta3 ~ dnorm(0, 0.0001)
    mu.beta4 ~ dnorm(0, 0.0001)
    mu.beta5 ~ dnorm(0, 0.0001)
    #gamma0b ~ dnorm(0, 0.0001)
    #gamma0b2 ~ dnorm(0, 0.0001)
    #gamma0b3 ~ dnorm(0, 0.0001)

    #gamma1b ~ dnorm(0, 0.0001)
    #gamma1b2 ~ dnorm(0, 0.0001)
    #gamma1b3 ~ dnorm(0, 0.0001)
 
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill = TRUE)
sink()

# Initial values
inits <- function(){
  list(mu.alpha = rnorm(1), 
       mu.beta1 = rnorm(1), 
       mu.beta2 = rnorm(1), 
       mu.beta3 = rnorm(1), 
       mu.beta4 = rnorm(1), 
       mu.beta5 = rnorm(1),
       BB = matrix(rnorm(J*K),nrow=J,ncol=K),
       Tau.B = rwish(K+1,diag(K))
       #gamma0b = rnorm(1), 
       #gamma1b=rnorm(1),
       #gamma0b2=rnorm(1),
       #gamma0b3=rnorm(1),
       #gamma1b2=rnorm(1),
       #gamma1b3=rnorm(1) 
  )
}

# Parameters monitored
parameters <- c("mu.alpha","mu.beta1","mu.beta2","mu.beta3", "mu.beta4",
                "mu.beta5", "BB", "Sigma.B"
                #"gamma0b","gamma1b","gamma0b2",
                #"gamma0b3","gamma1b2","gamma1b3"
)


# MCMC settings

## These can be changed as needed, up to what you need

ni <- 40000
nt <- 2
nb <- 5000
nc <- 3

out1 <- jags(data, inits, parameters, "model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

print(out1)
## if it is false, the mean doesn't cross 0, so it's what we would normally call significant, but here is "strongly correlated"


##make a logistic regression figure for all
## this has been the best way I've found to visualize the different parameters all together
##it can be really overwhelming at first, but just walk through it slowly and the patterns are clear
View(df6$species)
table(df6$species, df6$sample_health_state)


f1 <- ggplot(df6, aes(x = mean_prop_exo, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = .8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.0,0.7)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "OFRA", "MCAV")), ncol = 5)


f1

f2 <- ggplot(df6, aes(x = mean_gastro, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm", se=TRUE, fullrange=TRUE, 
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(0.6,0.9)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "OFRA", "MCAV")), ncol = 5)
f2

f3 <- ggplot(df6, aes(x = max_symb_vac, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.4, 1)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "OFRA", "MCAV")), ncol = 5)

f3

f4 <- ggplot(df6, aes(x = max_avg_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(20,40)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "OFRA", "MCAV")), ncol = 5)

f4

f5 <- ggplot(df6, aes(x = mean_degr_symb, y = state)) +
  geom_point(alpha = 0.5) +
  stat_smooth(method="glm",
              method.args = list(family=binomial), level = 0.8) +
  theme_classic(base_size = 10) + scale_x_continuous(limits=c(.7,1)) +
  facet_wrap(~factor(species, levels = c("CNAT", "OFAV", "OFRA", "MCAV")), ncol = 5)
f5

f1/f3/f4

f1/f2/f3/f4/f5


## strong correlation can be limiting, so we added an 80% CI to these studies as well
## when writing manuscripts we do strong correlation for 95% CI (above) and weak correlation for 80%

##by species here we go
##CNAT
quantile(out1$sims.list$BB[ ,1,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,1,6], probs = c(0.1, 0.9))

##OFAV
quantile(out1$sims.list$BB[ ,2,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,2,6], probs = c(0.1, 0.9))

##OFRA
quantile(out1$sims.list$BB[ ,3,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,3,6], probs = c(0.1, 0.9))

##PAST
quantile(out1$sims.list$BB[ ,4,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,4,6], probs = c(0.1, 0.9))

##PSTR
quantile(out1$sims.list$BB[ ,5,1], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,2], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,3], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,4], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,5], probs = c(0.1, 0.9))

quantile(out1$sims.list$BB[ ,5,6], probs = c(0.1, 0.9))



