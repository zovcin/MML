GRAD_TOL <- .Machine$double.eps^(3/4)
STEP_TOL <- .Machine$double.eps^(2/3)

grad <- function(f, x, h=10^(-3)) {
  n = length(x)
  g = numeric(n)
  for (i in 1:n) {
    fx = f(x)
    bak = x[i]
    x[i] = x[i]+h
    fxh = f(x)
    x[i] = bak
    g[i] = (fxh-fx)/h
  }
  g
}

optimize <- function(f, x0, max_iter) {
  i = 0
  x = x0

  while (i < max_iter) {
    g = grad(f, x)
    if (sqrt(sum(g^2)) <= GRAD_TOL) {
      print("gradient brake")
      return(x)
    }
    
    ret = line_search(f, x, g, -g, 10^(-4), 0.1)
    x = x-ret$alpha*g;
    if (ret$code == 1) {
      #print("max iter reached")
    } else if (ret$code == 2) {
      # print("small alpha")
      return(x)
    } else if (ret$code != 0) {
      stop("unhandled case")
    } 
    i = i+1
  }
  return(x)
}

get_minalpha <- function(x0, p) {
  rellength <- 0
  for (i in 1:length(x0)) {
    rellength <- max(rellength, abs(p[i]) / max(abs(x0[i]), 1))
  }
  STEP_TOL / rellength
}

line_search <- function(f, x0, g0, p, c1, c2) {
  i = 0
  MAX_ITER = 100
  
  alpha_prev = 0
  alpha = 1
  min_alpha = get_minalpha(x0, p)
  
  y0 = f(x0)
  y0grad = g0
  slope0 = sum(y0grad * p)
  y_prev = NULL
  
  while (TRUE) {
    if (i == MAX_ITER) { return(list(code=1, alpha=alpha)) }
    if (min_alpha > alpha) { return(list(code=2, alpha=alpha)) }
    i = i + 1
    
    x = x0+alpha*p
    y = f(x)

    # sufficient decrease check
    if (!(y <= y0+c1*alpha*slope0)) {
      alpha_next = zoom(y0, y, y_prev, slope0, alpha_prev, alpha)
    } else {
      # curvature decrease check
      slope = sum(grad(f, x)*p)
      if (!(abs(slope) <= c2*abs(slope0))) {
        if (slope >= 0) {
          alpha = zoom(y0, y, y_prev, slope0, alpha_prev, alpha)
          break
        } else {
          alpha_next = alpha*1.01
        }  
      }
    }
    
    alpha_prev = alpha
    alpha = alpha_next
    y_prev = y
  }
  
  return(list(code=0, alpha=alpha))
}

zoom <- function(y0, y, y_prev, slope0, alpha_prev, alpha) {
  if (alpha == 1) {
    alpha_next = inter2_optim(y0, y, slope0)
  } else {
    alpha_next = inter3_optim(y0, y_prev, y, slope0, alpha_prev, alpha)
  }
  
  if (alpha_next <= 0.1*alpha) {
    alpha_next = 0.1*alpha; 
  }
  
  alpha_next
}

inter2_optim <- function(y0, y, slope0) {
  -slope0 / (2*(y-y0-slope0))
}

inter3_optim <- function(y0, y_prev, y, slope0, alpha_prev, alpha) {
  m0 = c(1/(alpha - alpha_prev))
  
  m1 = matrix(byrow=TRUE, nrow = 2, c(
    1/alpha^2, -1/alpha_prev^2,
    -alpha_prev/alpha^2, alpha / alpha_prev^2
  ))
  
  m2 = matrix(byrow=TRUE, nrow=2, c(
    y - y0 - slope0*alpha,
    y_prev - y0 - slope0*alpha_prev
  ))
  
  t = m0 * (m1 %*% m2)
  
  a = t[1]
  b = t[2]
  if (a == 0) {
    alpha_next = -slope0 / 2*b
  } else {
    disc = b^2 - 3*a*slope0
    alpha_next = (-b + sqrt(disc)) / 3*a
  }
  
  if (alpha_next > 0.5*alpha) {
    alpha_next = 0.5*alpha;
  }
  
  alpha_next
}

M = 5
C = c(0.1, 0.2, 0.2, 0.4, 0.6)
A = matrix(nrow=5, byrow=TRUE, c(
  4.0, 4.0, 4.0, 4.0,
  1.0, 1.0, 1.0, 1.0,
  8.0, 8.0, 8.0, 8.0,
  6.0, 6.0, 6.0, 6.0,
  3.0, 7.0, 3.0, 7.0
))

shekel <- function(x) {
  r = 0
  for (i in 1:M) {
    r = r + (C[i] + sum((x-A[i,])^2))^(-1)
  }
  -r
}

x = optimize(shekel, c(1,3,5,6), 250)
print(x)
print(shekel(x))