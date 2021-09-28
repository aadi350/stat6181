install.packages("plotly")

library(plotly)
library(ggplot2)

# Specify range of x and y values
xs <- seq(from = -2.0, to = 2.0, by = 0.5)
ys <- seq(from = -2.0, to = 2.0, by = 0.5)

#define the function
fx <- function (x, y) {
    return (exp(-(1/3)*(x^2) + x - y^3))
}
# find matrix of values in both directions and apply function
z <- outer(xs, ys, fx)

# plot graph
fig = plot_ly(z = z, x=xs, y=ys, type="contour", colorscale='Bluered',width = 800, height = 500)
fig = fig %>%
  layout(title='Contour plot', plot_bgcolor = "#e5ecf6", xaxis = list(title = 'x'), 
         yaxis = list(title = 'y'))

fig

fig = plot_ly(z = z, x=xs, y=ys, colorscale='Bluered', width=800, height=800)
fig <- fig %>% add_surface()
fig = fig %>%
  layout(title='Surface plot', plot_bgcolor = "#e5ecf6", xaxis = list(title = 'x'), 
         yaxis = list(title = 'y'))


fig

# redefine function so that optim() works
fx <- function (var) {
  return(exp(-(1/3)*var[1]^2 + var[1] - var[2]^3))
}

optim(
    c(0,0), 
    fx, 
    method='L-BFGS-B',
    lower=c(-2,-2), 
    upper=c(2,2))

library(ggplot2)
library(plotly)

weibull <- function(x, shape=1, scale=5) {
    (shape/scale)*(x/scale)^(shape-1)*exp(-(x/scale)^(shape))
}

x <- seq(0, 2, by=0.01)
f_x <- weibull(x, shape=5, scale=1)

#generate plot
data = data.frame(x, f_x)
fig <- plot_ly(data, x=x, y=f_x, name='2-param Weibull', type='scatter', mode='lines', width = 500, height = 500)

fig <- fig %>% layout(title = '2-Parameter Weibull Distribution', plot_bgcolor = "#e5ecf6")

fig

set.seed(1123)
b <- 5
eta <- 1
n <- 1000
x <- rweibull(n, shape=b, scale=eta)
# write data to output file for submission
write.table(x, file="weibull.txt", row.names=FALSE, col.names=FALSE)


i = seq(1, n, 1)
data = data.frame(i, x)

fig <- plot_ly(data, x=i, y=x, name='2-param Weibull', type='scatter', mode='lines')

fig

fig <- plot_ly(x = x, type = "histogram")
fig

#log likelihood

f <- function(t)
{
    # scale
    eta <- t[1]; 
    # shape
    b <- t[2];
    if( (eta > 0) & (b > 0) )
    {
        return( sum( dweibull(x, shape=b, scale=eta, log=TRUE) ) )
    } else
    {
        return(-Inf)
    }
}

df <- function(t) {
    eta <- t[1]
    b <- t[2]
    
    score = rep(0,2)
    score[1] = n/b - n*log(eta) + sum(log(x)) - sum((x/eta)^(b) * log(x/eta))    
    score[2] = -(n*b)/eta + (eta/b)*sum((x/eta)^b )
    
    return(score)
}

df2 <- function(t) {
    eta <- t[1]
    b <- t[2]
        
    h <- matrix(0,2,2)    
    h[1,1] <- (b/eta^2)*(n-(b-1)*sum((x/eta)^b))
    h[2,2] <- n/b^2 - sum((x/eta)^(b) * (log(x/eta))^2)
    h[1,2] <- -(1/eta)*(n - sum((x/eta)^b) - b*sum(log(x/eta)*(x/eta)^b) )
    h[2,1] <- h[1,2]
    
    return(h)
}

# from lecture notes
newton <- function(x0, f, df, d2f, tol=1e-4, pr=FALSE) {
    k <- 1
    fval <- f(x0) 
    grad <- df(x0)
    hess <- d2f(x0)
    xk_1 <- x0 
    cond1 <- sqrt( sum(grad^2) ) 
    cond2 <- Inf
    if( (cond1 < tol) ) return(x0) 
    
    while( (cond1 > tol) & (cond2 > tol)) {
        L <- 1
        bool <- TRUE 
        while(bool == TRUE) {
        xk <- xk_1 - L * solve(hess) %*% grad
            if( f(xk) > fval ) {
            bool = FALSE 
            grad <- df(xk) 
            fval <- f(xk)
            hess <- d2f(xk) 
            } else {
                L = L/2
                if( abs(L) < 1e-20 ) {
                    return("Failed to find uphill step - try new start values")
                }
            } 
        }
        cond1 <- sqrt( sum(grad^2) ) 
        cond2 <- sqrt( sum( (xk-xk_1)^2 ))/(tol + sqrt(sum(xk^2)))
        k <- k + 1
        xk_1 <- xk 
    }
    if(pr == TRUE) print( sprintf("Took %i iterations", k) )
    return(xk)
}

newton(x0=c(2, 4), f=f, df=df, d2f=df2, pr=TRUE)

# from lecture notes
newton <- function(x0, f, df, d2f, tol=1e-4, pr=FALSE) {
    k <- 1
    fval <- f(x0) 
    grad <- df(x0)
    hess <- d2f(x0)
    xk_1 <- x0 
    cond1 <- sqrt( sum(grad^2) ) 
    cond2 <- Inf
    if( (cond1 < tol) ) return(x0) 
    
    while( (cond1 > tol) & (cond2 > tol)) {
        L <- 1
        bool <- TRUE 
        while(bool == TRUE) {
        xk <- xk_1 - L * solve(hess) %*% grad
            if( f(xk) > fval ) {
            bool = FALSE 
            grad <- df(xk) 
            fval <- f(xk)
            hess <- d2f(xk) 
            } else {
                L = L/2
                if( abs(L) < 1e-20) {
                    return("Failed to find uphill step - try new start values")
                }
            } 
        }
        cond1 <- sqrt( sum(grad^2) ) 
        cond2 <- sqrt( sum( (xk-xk_1)^2 ))/(tol + sqrt(sum(xk^2)))
        k <- k + 1
        xk_1 <- xk 
    }
    if(pr == TRUE) print( sprintf("Took %i iterations", k) )
    return(xk)
}

newton(x0=c(2, 4), f=f, df=df, d2f=df2, pr=TRUE, tol=1e-10)
