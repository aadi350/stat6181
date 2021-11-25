#!/home/aadi/miniconda3/envs/r_env/bin/Rscript

# establish fixed values
v <- 30
sigma_sq <- 15 
mu <- 35

# install.packages("extraDistr") to attain rinvchisq
# install.packages("extraDistr")
library(extraDistr)

set.seed(117)
N <- 50

# generating random values from the inverse chi-squared 
# rinvchisq(n, nu, tau)
# V is a vector of V_i
V <- rinvchisq(N, v, sigma_sq^0.5)
Y <- rnorm(N, mu, sqrt(V))

sample_indep_mh <- function(rg, dg, df, cur) {
    u <- runif(1)

    proposed <- rg()
    ratio <- (df(proposed)*dg(cur)/df(cur)*dg(proposed))
    accept = u < ratio
    
    if(accept) {
        return(proposed)
    }
    return(cur)
}

sample_rw_mh <- function(x) {

}

# defining full conditionals  
sample_V <- function(y, mu, v, sigma_sq, N) {
    
    return(rinvchisq(N, 
              v+1, 
              sqrt(((y-mu)^2 + sigma_sq*v)/(v+1))))
}

sample_sigma_sq <- function(v, n, V_i) {
    return(rgamma(1,
                 shape=0.5*n*v + 1,
                  rate=(v/2)*sum(1/V_i)
                 ))
}

# MH Sampler for \mu    
sample_mu <- function(Y, V_i, sigma_sq, v, mu) {
    
    mu_mchain <- rep(NA, 10000)
    mu_mchain[1] <- mu

    "
        MH ratio = f(x*)g(x|x*)/g(x*)f(x*|x)
        rg : proposal generation
        dg : proposal density 
        df : f(x) 
    "

    # proposal random distribution 
    rg <- function() {
        rnorm(1,mu, 1)
    }

    #  
    df <- function(x) {
        exp(-0.5 * sum((Y-x)^2/V_i))
    }
    
    dg <- function(x) {
        dnorm(x, mu, 1)

    }
    
    for (i in 2:10000) {
            mu_mchain[i] <- sample_indep_mh(rg, dg, df, mu_mchain[i-1]) 
    } 
    
    return(mean(mu_mchain))
}

# Gibbs Sampling
ITER_NUM <- 100
data <- Y

mchain <- matrix(NA,ITER_NUM, 2)
# set starting values
# column1=mean, column2=sigma_sq
mchain[1, 1] <- 30 
mchain[1, 2] <- 20


V_i <- matrix(NA, ITER_NUM, N)
V_i[1,] <- rep(1, N) 

v <- 30

for(i in 2:ITER_NUM) {
    cur_V_i <- V_i[i-1,]
    
    cur_mu <- mchain[i-1, 1]
    cur_sigma_sq <- mchain[i-1, 2]
    

    new_V_i <- sample_V(data, cur_mu, v, cur_sigma_sq, N)
    
    # MH in Gibbs step    
    new_mu <- sample_mu(data, cur_V_i, cur_sigma_sq, v, cur_mu)
    
    new_sigma_sq <- sample_sigma_sq(v, N, cur_V_i)  

    V_i[i, ] <- new_V_i
    mchain[i, 1] <- new_mu
    mchain[i, 2] <- new_sigma_sq 
  
    
}

print("Mean MU")
print(mean(mchain[,1]))

print("Mean Sigma sq")
print(mean(mchain[,2]))

print("Differennce in V_i")
print(colMeans(V_i) - V)