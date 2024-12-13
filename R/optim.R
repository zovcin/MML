optim <- function(britn, x0, fn, gfn ) {
  # Iterative method for deterministic optimization problem
  #
  # Usage:
  # result <- optim(30, c(5, 5, 5, 5), fn, gfn)
  #
  # Inputs:
  # britn: maximum number of iterations
  # x0: initial point
  # fn: function to optimize
  # gfn: gradient of the function (if not provided, finite differences are used)
  #
  # Outputs:
  # List containing:
  #   x0: obtained final value
  #   y0: function at ther obtained final value
  #   g0: gradient value at the final point
  #   termcode: termination reason
  #     0 = reached maximum number of iterations
  #     1 = reached maximum function evaluations
  #     2 = gradient norm below gradtol
  #     3 = step size issues (stuck)
  #   itncount: number of iterations

  gradtol <- .Machine$double.eps^(3/4)

  alpha <- 0.085

  termcode <- 0 # Default termination: max iterations

  g0 <- gfn(x0)

  # Main loop
  for (itncount in 1:britn) {
    # Compute new point in the chosen direction
    x1 <- x0 - alpha * t(g0)
    
    # Check if gradient norm is below tolerance
    if (sqrt(sum(g0^2)) <= gradtol) {
      termcode <- 2
      break
    }

    # Prepare for next iteration
    x0 <- x1
    g0 <- gfn(x0)
  }

  y0 = fn(x0)
  
  return(list(x0 = x0, y0 = y0, g0 = g0, termcode = termcode, itncount = itncount))
}

# Deterministic function
fn <- function(x) {
  a <- matrix(c(2,1,1,20), ncol = 2)
  b <- matrix(c(5,3), ncol = 1)
  y <- t(x) %*% a %*% x - t(b) %*% x
  return(y)
}

# Gradient of deterministic function
gfn <- function(x) {
  a <- matrix(c(2,1,1,20), ncol = 2)
  b <- matrix(c(5,3), ncol = 1)
  g <- t(x) %*% a - t(b)
  return(g)
}

result <- optim(britn = 200, x0 = matrix(c(-3,-1),ncol = 1), fn = fn, gfn = gfn)
print(result)
result$x0  # Final value of x
result$y0  # Final function value
result$g0  # Gradient at the final point
result$termcode  # Termination reason
result$itncount  # Number of iterations
