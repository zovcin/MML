optim02 <- function(britn, x0, fn, gfn = NULL) {
  # Iterative method for deterministic optimization problem
  #
  # Usage:
  # result <- optim02(30, c(5, 5, 5, 5), f4, gf4)
  #
  # Inputs:
  # britn: maximum number of iterations
  # x0: initial point
  # fn: function to optimize
  # gfn: gradient of the function (if not provided, finite differences are used)
  #
  # Outputs:
  # List containing:
  #   x1: obtained final value
  #   y1: value of the objective function at the final point
  #   g1: gradient value at the final point
  #   termcode: termination reason
  #     0 = reached maximum number of iterations
  #     1 = reached maximum function evaluations
  #     2 = gradient norm below gradtol
  #     3 = step size issues (stuck)
  #   itncount: number of iterations

  # Check gradient input
  analytic_grad <- !is.null(gfn)
  if (missing(fn)) {
    stop("Insufficient number of input arguments")
  }

  # Set initial parameters
  n <- length(x0)
  maxstep <- 1000 * max(sqrt(sum(x0^2)), 1)
  steptol <- .Machine$double.eps^(2/3)
  gradtol <- .Machine$double.eps^(3/4)
  eta <- .Machine$double.eps
  maxfcalcls <- 100

  # Save initial x0 and y0
  itncount <- 1
  y0 <- fn(x0)

  # Compute gradient
  if (analytic_grad) {
    g0 <- gfn(x0)
  } else {
    g0 <- fdgrad(x0, y0, fn, eta)
  }

  termcode <- 0 # Default termination: max iterations

  # Main loop
  for (itncount in 1:(britn - 1)) {
    # Search for descending direction
    d <- -g0

    # Compute new point in the chosen direction
    result <- ls(x0, y0, fn, g0, d, maxstep, steptol, maxfcalcls)
    x1 <- result$x1
    y1 <- result$y1
    retcode <- result$retcode
    lambda <- result$lambda
    maxtaken <- result$maxtaken
    d <- result$d

    # Check termination condition
    if (retcode == 1) { # Could not find step size (lambda)
      termcode <- 3
      break
    }

    # Compute gradient at the new point
    if (analytic_grad) {
      g1 <- gfn(x1)
    } else {
      g1 <- fdgrad(x1, y1, fn, eta)
    }

    # Check if gradient norm is below tolerance
    if (sqrt(sum(g1^2)) <= gradtol) {
      termcode <- 2
      break
    }

    # Prepare for next iteration
    x0 <- x1
    g0 <- g1
    y0 <- y1
  }

  itncount <- itncount + 1

  return(list(x1 = x1, y1 = y1, g1 = g1, termcode = termcode, itncount = itncount))
}

# Forward difference gradient approximation
fdgrad <- function(xc, fc, fn, eta) {
  n <- length(xc)
  g <- numeric(n)
  sqrteta <- sqrt(eta)

  for (j in 1:n) {
    if (xc[j] == 0) {
      stepsizej <- sqrteta
    } else {
      stepsizej <- sqrteta * max(abs(xc[j]), 1) * sign(xc[j])
    }
    tempj <- xc[j]
    xc[j] <- xc[j] + stepsizej
    stepsizej <- xc[j] - tempj
    fj <- fn(xc)
    g[j] <- (fj - fc) / stepsizej
    xc[j] <- tempj
  }

  return(g)
}

# Line search function
ls <- function(x0, y0, fn, g0, d, maxstep, steptol, maxfcalcls) {
  n <- length(x0)
  maxtaken <- FALSE
  retcode <- 3
  fcalcls <- 0
  alpha <- 0.0001

  Newtlen <- sqrt(sum(d^2))
  if (Newtlen > maxstep) {
    d <- d * (maxstep / Newtlen)
    Newtlen <- maxstep
  }

  initslope <- sum(g0 * d)

  rellength <- 0
  for (i in 1:n) {
    rellength <- max(rellength, abs(d[i]) / max(abs(x0[i]), 1))
  }

  minlambda <- steptol / rellength
  lambda <- 1

  if (maxfcalcls == 0) {
    x1 <- x0
    y1 <- y0
    retcode <- 2
    return(list(x1 = x1, y1 = y1, retcode = retcode, lambda = lambda, maxtaken = maxtaken, d = d, fcalcls = fcalcls))
  }

  while (retcode > 2) {
    x1 <- x0 + lambda * d
    y1 <- fn(x1)
    fcalcls <- fcalcls + 1

    if (y1 <= y0 + alpha * lambda * initslope) {
      retcode <- 0
      if ((lambda == 1) && (Newtlen > 0.99 * maxstep)) {
        maxtaken <- TRUE
      }
    } else if (lambda < minlambda) {
      retcode <- 1
      x1 <- x0
      y1 <- y0
    } else if (fcalcls >= maxfcalcls) {
      retcode <- 2
    } else {
      if (lambda == 1) {
        lambda_temp <- -initslope / (2 * (y1 - y0 - initslope))
      } else {
        tempm <- c(1 / (lambda - lambda_prev)) *
          matrix(c(1 / lambda^2, -1 / lambda_prev^2,
                   -lambda_prev / lambda^2, lambda / lambda_prev^2),
                 nrow = 2, byrow = T) %*%
          matrix(c(y1 - y0 - lambda * initslope, y_prev - y0 - lambda_prev * initslope), ncol = 1)
        a <- tempm[1]
        b <- tempm[2]
        disc <- b^2 - 3 * a * initslope

        if (a == 0) {
          lambda_temp <- -initslope / (2 * b)
        } else {
          lambda_temp <- (-b + sqrt(disc)) / (3 * a)
        }

        if (lambda_temp > 0.5 * lambda) {
          lambda_temp <- 0.5 * lambda
        }
      }

      lambda_prev <- lambda
      y_prev <- y1

      if (lambda_temp <= 0.1 * lambda) {
        lambda <- 0.1 * lambda
      } else {
        lambda <- lambda_temp
      }
    }
  }

  return(list(x1 = x1, y1 = y1, retcode = retcode, lambda = lambda, maxtaken = maxtaken, d = d, fcalcls = fcalcls))
}

# Deterministic function
fn <- function(x) {
  n <- length(x)
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
  return(c(y))
}

n <- 4
result <- optim02(britn = 500, x0 = matrix(1+0*(1:n),ncol = 1)/5, fn = fn)
print(result)
result$x1  # Final value of x
result$y1  # Final function value
result$g1  # Gradient at the final point
result$termcode  # Termination reason
result$itncount  # Number of iterations
