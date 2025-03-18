# Preface ----------------------------------------------------------------------
# MOSRD Example: code is from Oliver Gimenez's "bayesRD" folder on gitHub 
# link: https://github.com/oliviergimenez/bayesRD#bayesian-implementation-of-capture-recapture-models-with-robust-design
# Tinkering with this code to build my knowledge in Bayesian stats and start learning about the robust design

# Initialization ---------------------------------------------------------------
library(R2jags) #to run JAGS
library(shinystan) #to run shiny stan
library(tidyverse) #to utilize pipe operators
library(RMark) #to load the test data set

# Create Data ------------------------------------------------------------------
#load data
data(robust)

#create time intervals
time.intervals <- c(0, 1,            # primary occasion 1, 2 secondary occasions
                    0, 1,            # primary occasion 2, 2 secondary occasions
                    0, 0, 0, 1,      # primary occasion 3, 4 secondary occasions
                    0, 0, 0, 0, 1,   # primary occasion 4, 5 secondary occasions
                    0)               # primary occasion 5, 2 secondary occasions

chx <- NULL
for (i in 1:nrow(robust)){
  chx <- rbind(chx, unlist(str_split(robust$ch[i],'')))
}
head(chx)

#Create encounter histories
encounter <- NULL 
for (i in 1:length(robust$freq)){
  if (robust$freq[i] == 1) encounter <- rbind(encounter, as.numeric(chx[i,]))
  if (robust$freq[i] > 1) encounter <- rbind(encounter, matrix(rep(as.numeric(chx[i,]), robust$freq[i]), nrow = robust$freq[i], byrow= T))
}
head(encounter)

#Calculate important quantities
n.ind <- nrow(encounter) # number of individuals
n.secondary <- c(2, 2, 4, 5, 2) # number of secondary occasions per primary occasion
n.primary <- length(n.secondary) # number of primary occasions
index <- list(1:2,
              3:4,
              5:8,
              9:13,
              14:15) # the secondary occasions

#Calculate the number of individuals caught in each primary occasion
caught <- rep(NA, n.primary)
for (i in 1:n.primary){
  tmp <- encounter[,index[[i]]]
  caught[i] <- nrow(tmp[rowSums(tmp)!=0,])
}
caught

#Observation array: number of individuals, number of primary occasions, number of secondary occasions
obs <- array(NA, dim = c(n.ind, n.primary, max(n.secondary)))
for (i in 1:n.primary){
  obs[,i,1:n.secondary[i]] <- encounter[,index[[i]]]
}
dim(obs)

# Format Data ------------------------------------------------------------------
#format data for the Bayesian robust design
ch <- matrix(NA, n.ind, n.primary)
for (i in 1:n.ind){
  for (t in 1:n.primary){
    ifelse(any(obs[i,t,1:n.secondary[t]] == 1), ch[i,t] <- 1, ch[i,t] <- 2)
  }
}

#summarize detections by primary and secondary occassions
test <- matrix(NA, n.ind, n.primary)
for (i in 1:nrow(test)){
  for (j in 1:ncol(test)){
    test[i,j] <- sum(obs[i,j,], na.rm = TRUE)
  }
}

seen <- array(NA, c(n.ind, n.primary, max(n.secondary)))
missed <- array(NA, c(n.ind, n.primary, max(n.secondary)))

for (i in 1:nrow(test)){
  for (t in 1:ncol(test)){
    for (j in 1:n.secondary[t]){
      if(test[i,t] > 1 & obs[i,t,j] == 1){seen[i,t,j] <- 1}
      if(test[i,t] >= 1 & obs[i,t,j] == 0){missed[i,t,j] <- 1}
    }
  }
}

yes <- matrix(NA, n.primary, max(n.secondary))
no <- matrix(NA, n.primary, max(n.secondary))

for (i in 1:nrow(yes)){
  for (j in 1:ncol(yes)){
    yes[i,j] <- sum(seen[,i,j], na.rm = TRUE)
    no[i,j] <- sum(missed[,i,j], na.rm = TRUE)
  }
}

total <- yes + no
total

#first occasion of capture
get.first <- function(x)min(which (x != 2))
first <- apply(ch,1,get.first); first[first == "Inf"] <- NA

#remove individuals released in last primary occasion
ch <- subset(ch, first != n.primary)
first <- subset(first, first != n.primary)

#define initial values for the latent states
z.init <- matrix(NA, nrow(ch), ncol(ch))
for (i in 1:nrow(ch)){
  if(first[i] < ncol(z.init)){
    z.init[i,(first[i] + 1):ncol(z.init)] <- 1
  }
}

#Constructing JAGS model
model <- function() { 
  
  # priors
  phi ~ dunif(0,1)     # survival 
  gP ~ dunif(0,1)      # gamma'
  gPP ~ dunif(0,1)     # gamma'' 
  gamP <- 1 - gP       # MARK parameterization
  gamPP <- 1 - gPP     # MARK parameterization
  mean.p ~ dunif(0,1)  # detection
  
  # secondary occasions p's
  for (t in 1:n.years){
    for (j in 1:max(n.sec[1:n.years])){
      p[t,j] <- mean.p
    }
  }   
  
  for (t in 1:n.years){
    for (j in 1:n.sec[t]){
      yes[t,j] ~ dbin(p[t,j], total[t,j])
    }
  }   
  
  # Primary occasions p's or pooled detection probability
  for (t in 1:n.years){
    pstar[t] <- 1 - prod(1 - p[t,1:n.sec[t]])
  }
  
  # state matrices
  s[1,1] <- phi * gPP
  s[1,2] <- phi * (1 - gPP)
  s[1,3] <- 1 - phi
  s[2,1] <- phi * gP
  s[2,2] <- phi * (1 - gP)
  s[2,3] <- 1 - phi
  s[3,1] <- 0
  s[3,2] <- 0
  s[3,3] <- 1
  
  # observation matrices
  for (t in 1:n.years){
    o[1,t,1] <- pstar[t]
    o[1,t,2] <- 1 - pstar[t]
    o[2,t,1] <- 0
    o[2,t,2] <- 1
    o[3,t,1] <- 0
    o[3,t,2] <- 1
  }
  
  # likelihood
  for (i in 1:n.ind){
    z[i,first[i]] <- ch[i,first[i]]
    for (t in (first[i]+1):n.years){
      z[i,t] ~ dcat(s[z[i,t-1], ])   # state equations
      ch[i,t] ~ dcat(o[z[i,t], t, ]) # obsevation equations
    } 
  }
  
}

# Model analysis ---------------------------------------------------------------
#create the list of data
dat <- list(first = first, 
            ch = ch, 
            n.sec = n.secondary, 
            n.years = ncol(ch), 
            n.ind = nrow(ch),
            yes = yes, 
            total = total)

#Set initial values
inits <- function(){list(z = z.init)}



#define parameters
pars <- c('pstar','mean.p','phi','gamP','gamPP')

#mcmc settings 
n.chains <- 1
n.iter <- 1000
n.burnin <- 500

#Call JAGS from R
res_markovian <- jags(data = dat, 
                      inits = inits, 
                      parameters.to.save = pars, 
                      model.file = model, 
                      n.chains = n.chains,
                      n.iter = n.iter, 
                      n.burnin = n.burnin)

print(res_markovian, digits = 3)

library(lattice)
jagsfit.mcmc <- as.mcmc(res_markovian)
densityplot(jagsfit.mcmc)

#k<-mcmcplots::as.mcmc.rjags(res_markovian)%>%as.shinystan()%>%launch_shinystan() #making it into a MCMC, each list element is a chain, then puts it through to shiny stan

