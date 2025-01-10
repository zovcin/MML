optim03 <- function(britn, x0, fn, gfn, hfn) {
  # Iterative method for deterministic optimization problem
  # 
  # Input:
  #   britn - maximum number of iterations
  #   x0 - starting point
  #   fn - objective function
  #   gfn - gradient of the objective function
  #   hfn - Hessian of the objective function
  #
  # Output:
  #   x1 - final value
  #   y1 - function value at the final point
  #   g1 - gradient at the final point
  #   termcode - termination reason
  #              0 = reached maximum number of iterations
  #              2 = gradient below tolerance
  #   itncount - number of iterations
  
  gradtol <- .Machine$double.eps^(3 / 4)  # Gradient tolerance
  
  # Compute initial gradient
  g0 <- gfn(x0)
  
  # Check for early exit
  if (sqrt(sum(g0^2)) <= gradtol) {
    termcode <- 2
    return(list(x1 = x0, y1 = fn(x0), g1 = g0, termcode = termcode, itncount = 0))
  }
  
  termcode <- 0  # Default termination code
  
  # Loop
  for (itncount in 1:britn) {
    # Compute Hessian
    H0 <- hfn(x0)
    
    # Compute descent direction
    d <- -solve(H0, g0)
    
    # Update point
    x1 <- x0 + d
    
    # Compute new gradient
    g1 <- gfn(x1)
    
    # Check for termination
    if (sqrt(sum(g1^2)) <= gradtol) {
      termcode <- 2
      break
    }
    
    # Prepare for next iteration
    x0 <- x1
    g0 <- g1
  }
  
  # Final function value
  y1 <- fn(x1)
  
  return(list(x1 = x1, y1 = y1, g1 = g1, termcode = termcode, itncount = itncount))
}

# Deterministic function
fn <- function(x) {
  n <- 4
  a <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in i:n) {
      a[i, j] <- 1
    }
  }
  a <- a / n
  ax <- a %*% x
  ax3 <- ax^3
  ax4 <- ax^4
  
  y <- t(x) %*% (t(a) %*% a) %*% x + sum(ax3) / 10 + sum(ax4) / 100
  return(as.numeric(y))
}

# Gradient of the function
gfn <- function(x) {
  g=c(
    (2*x[1]+2*x[2]+2*x[3]+2*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3),
    (2*x[1]+4*x[2]+4*x[3]+4*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(3/10)*(x[2]+x[3]+x[4])^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3+(1/25)*(x[2]+x[3]+x[4])^3),
    (2*x[1]+4*x[2]+6*x[3]+6*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(3/10)*(x[2]+x[3]+x[4])^2+(3/10)*(x[3]+x[4])^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3+(1/25)*(x[2]+x[3]+x[4])^3+(1/25)*(x[3]+x[4])^3),
    (2*x[1]+4*x[2]+6*x[3]+8*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(3/10)*(x[2]+x[3]+x[4])^2+(3/10)*(x[3]+x[4])^2+(3/10)*x[4]^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3+(1/25)*(x[2]+x[3]+x[4])^3+(1/25)*(x[3]+x[4])^3+(1/25)*x[4]^3)
  )
  return(matrix(g,ncol=1))
}

# Hessian of the function
hfn <- function(x) {
  g=c(2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2, 2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2, 2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2, 2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2,
    2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2, 4+(3/5)*x[1]+(6/5)*x[2]+(6/5)*x[3]+(6/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2, 4+(3/5)*x[1]+(6/5)*x[2]+(6/5)*x[3]+(6/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2, 4+(3/5)*x[1]+(6/5)*x[2]+(6/5)*x[3]+(6/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2,
    2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2, 4+(3/5)*x[1]+(6/5)*x[2]+(6/5)*x[3]+(6/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2, 6+(3/5)*x[1]+(6/5)*x[2]+(9/5)*x[3]+(9/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2+(3/25)*(x[3]+x[4])^2, 6+(3/5)*x[1]+(6/5)*x[2]+(9/5)*x[3]+(9/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2+(3/25)*(x[3]+x[4])^2,
    2+(3/5)*x[1]+(3/5)*x[2]+(3/5)*x[3]+(3/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2, 4+(3/5)*x[1]+(6/5)*x[2]+(6/5)*x[3]+(6/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2, 6+(3/5)*x[1]+(6/5)*x[2]+(9/5)*x[3]+(9/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2+(3/25)*(x[3]+x[4])^2, 8+(3/5)*x[1]+(6/5)*x[2]+(9/5)*x[3]+(12/5)*x[4]+(3/25)*(x[1]+x[2]+x[3]+x[4])^2+(3/25)*(x[2]+x[3]+x[4])^2+(3/25)*(x[3]+x[4])^2+(3/25)*x[4]^2)
  return(matrix(g,ncol=4))
}

n <- 4

result <- optim02(50000, matrix(1+0*(1:n),ncol = 1)/5, fn)
result <- optim03(britn = 50000, matrix(1+0*(1:n),ncol = 1)/5, fn, gfn, hfn)
print(result)
result$x1  # Final value of x
result$y1  # Final function value
result$g1  # Gradient at the final point
result$termcode  # Termination reason
result$itncount  # Number of iterations
