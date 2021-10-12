# QUESTION I
## Derivation of the EM algorithm

# Expectation Step
e_step <- function(x, mean_1, mean_2, pi) {
    return(
        (pi * dpois(x, mean_2, log=FALSE))/
        ((1-pi)*dpois(x, mean_1, log=FALSE) + pi*dpois(x, mean_2, log=FALSE))
    )
}
# Maximization Step
m_step <- function(x, gamma) {
    mean_1 = sum(x*(1-gamma))/sum(1-gamma)
    mean_2 = sum(x*gamma)/sum(gamma)

    return(c(mean_1, mean_2))
}

# Putting together E and M using function from lecture slides
e_m <- function(x, mean_1_0, mean_2_0, pi_0, tol=1e-10, max_iter=1000) {
    # specifies change of means at every iteration
    delta <- Inf
    # wrap both means into vector
    means <- c(mean_1_0, mean_2_0)
    pi <- pi_0
    iter <- 1
    
    
    while (delta > tol) {
        # update new gamma values in expectation step
        gamma_prev <- e_step(x, means[1], means[2], pi)
        # use updated gamma to maximize both means in  maximization step
        new_means <- m_step(x, gamma_prev)
        means <- new_means
        # update pi (proportion) 
        pi <- sum(gamma_prev)/length(x)
        
        # calculate change in old and new estimates 
        delta <- abs(new_means[1]-means[1] + new_means[2]-means[2])
        iter <- iter + 1
        
        
        # if gradients explode, return error
        if (is.na(delta)) {
            return("ERROR: EXPLOSION")
        }
        
        # if more than max_iter or change is below tolerance, stop and return
        if (iter >= max_iter | delta < tol) {
             return(c(new_means, pi))  
        }
 
                
    }
    return(c(new_means, pi))  
}
# place given X values into vector
x = c(1, 2, 3, 8 ,12)
# running single iteration of E-M
#  first find estimates for gamma
gamma <- e_step(x, 2, 8, 0.5)
# perform maximization step
m_step(x, gamma)

# Generating 100 values from Poisson with mean 5 and 900 with mean 10
# 	and running EM
set.seed(12020569)

# Set means for both distributions
means <- c(5, 10)

# let number of values generated be n
n <- c(100, 900)

# generate 1000 samples from distributions in one column
x <- c(rpois(n[1], means[1]), rpois(n[2], means[2]))

e_m(x, 2, 7, 0.5)

# QUESTION 2
set.seed(46692)

N = 100000
# generate N Uniform Samples 
uniform_X <- runif(N, min=0, max=1)
uniform_Y <- runif(N, min=0, max=1)

# create matrix to store resulting samples
joint <- matrix(nrow=N, ncol=2)


for (i in seq_len(N)) {
    u_X = uniform_X[i]
    u_Y = uniform_Y[i]
        # getting the X value, using the CDF to determine which X value to 'choose'
    # if the value from the uniform is less than 5/20, choose X=1
    if (u_X < 5/20) {
        joint[i, 1] = 1
    } 
    # if the value from the uniform is greater than 5/20 and less than 14/20, choose X=2
    else if (u_X >= 5/20 && u_X < 14/20) {
        joint[i, 1] = 2
    } 
    # if the value from the uniform is greather than 14/20, choose X=3
    else {
       joint[i, 1] = 3
    }
    
    # getting the Y value, using the CDF to determine which Y value to 'choose'
    # if the value from the uniform is less than 4/20, choose Y=1
    if (u_Y < 4/20) {
        joint[i, 2] = 1
    } 
    # if the value from the uniform is greater than 5/20, choose Y=2
    else {
        joint[i, 2] = 2
    }
}

# for Y
table(joint[,1])/N

# for X
table(joint[,2])/N

set.seed(46692)

P <- matrix(nrow=2, ncol=3)
P[1, ] <- c(1/20,2/20,1/20) 
P[2, ] <- c(4/20, 7/20, 5/20)


# using joint probabilities to find expectation
#  this is the sample mean 
expectation = 0
for (idx in seq_len(N)) {
    x <- joint[idx, 1] 
    y <- joint[idx, 2]
    P_xy = P[y,x]
    expectation = expectation + x * y
}
expectation = expectation/N
expectation

P[1,1]

expectation/100000
