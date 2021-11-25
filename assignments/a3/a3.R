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

# defining full conditionals  
sample_V <- function(Y, mu, v, sigma_sq, N) {
    "return(rinvchisq(N, 
              v+1, 
              sqrt(((Y-mu)^2 + sigma_sq*v)/(v+1))))"
    sample <- (1/rgamma(N,
        0.5*(v+1),
        0.5*((Y-mu)^2+v*sigma_sq)
        ) # close rinvgamma
    ) # close return
   return(sample)
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
V_i[1,] <- rnorm(N, v, 1) 

v <- 30

for(i in 2:ITER_NUM) {
    cur_V_i <- V_i[i-1,]
    
    cur_mu <- mchain[i-1, 1]
    cur_sigma_sq <- mchain[i-1, 2]
    

    new_V_i <- sample_V(data, cur_mu, v, cur_sigma_sq, N)
    print(new_V_i) 
    # MH in Gibbs step    
    new_mu <- sample_mu(data, cur_V_i, cur_sigma_sq, v, cur_mu)
    
    new_sigma_sq <- sample_sigma_sq(v, N, cur_V_i)  

    V_i[i, ] <- new_V_i
    mchain[i, 1] <- new_mu
    mchain[i, 2] <- new_sigma_sq 
  
    
}

print("MU")
print(mean(mchain[,1]))



print("Sigma sq")
print(mean(mchain[,2]))

print("Differennce in V_i")
print(colMeans(V_i)-V)

png("sigma_sq_acf.png")
acf_sigma_sq <- acf(mchain[,2], plot=FALSE)
plot(acf_sigma_sq, main="Sigma-Squared ACF")
dev.off()

png("mu_acf.png")
acf_mu <- acf(mchain[,1], plot=FALSE)
plot(acf_mu, main="mu ACF")
dev.off()

png("mu_dist.png")
dist_mu <- density(mchain[,1], plot=FALSE)
plot(dist_mu, main="mu Density")
dev.off()

png("sigma_sq_dist.png")
dist_sigma_sq <- density(mchain[,2], plot=FALSE)
plot(dist_sigma_sq, main="Sigma Sq Density")
dev.off()