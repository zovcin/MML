theta_fn <- function(x1, x2) {
    if (abs(x1) < 1e-10) {  
        if (x2 >= 0) return(0.25)
        else return(0.75)
    }
    if (x1 > 0) {
        return((1/(2*pi)) * atan(x2/x1))
    } else {
        return(0.5 + (1/(2*pi)) * atan(x2/x1))
    }
}
f1 <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    x3 <- x[3]
    return(10 * (x3 - 10 * theta_fn(x1, x2)))
}

f2 <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    return(10 * (sqrt(x1^2 + x2^2) - 1))
}

f3 <- function(x) {
    return(x[3])
}

Fhelical <- function(x) {
    return(f1(x)^2 + f2(x)^2 + f3(x)^2)
}

grad_f1 <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    #x3 <- x[3]
    dx1 = 100*x2/(2*pi*(x1^2+x2^2))
    dx2 = -100*x1/(2*pi*(x1^2+x2^2))
    dx3 = 10
    return( c(dx1, dx2, dx3))
}

grad_f2 <- function(x) {
    x1 <- x[1]
    x2 <- x[2]
    dx1 = 10*x1/(sqrt(x1^2+x2^2))
    dx2 = 10*x2/(sqrt(x1^2+x2^2))
    dx3 = 0
    return(c(dx1, dx2, dx3))
}

grad_f3 <- function(x) {
    return(c(0,0,1))
}


grad_fhelical <- function(x) {
    f1x <- f1(x)
    f2x <- f2(x)
    f3x <- f3(x)
    gf1 <- grad_f1(x)
    gf2 <- grad_f2(x)
    gf3 <- grad_f3(x)
    return( 2*f1x*gf1 + 2*f2x*gf2 + 2*f3x*gf3)
}

line_search <- function(x0, f0, fn, g0, p,
                             alpha_min   = 1e-6,  
                             max_f_evals = 30,    
                             c           = 1e-4    
) {

  g0 <- grad_fhelical(x0)
  initslope <- sum(g0 * p)  
  retcode   <- 1            
  alpha     <- 1.450 
  f_evals   <- 0
  x1        <- x0
  f1        <- f0

  alpha_prev <- NA  
  f_prev     <- NA  

  while (retcode > 0) {
    
    x_trial <- x0 + alpha * p
    f_trial <- fn(x_trial)
    f_evals <- f_evals + 1
    #print(paste("Iteration:", f_evals, 
           # "Alpha:", alpha, 
           # "f_trial:", f_trial, 
           # "Target:", f0 + c * alpha * initslope))

    #Armijo
    if (f_trial <= f0 + c * alpha * initslope) {

      retcode <- 0
      x1 <- x_trial
      f1 <- f_trial
      #print(paste("retcode:", retcode))
      break

    } else if (alpha < alpha_min) {
      retcode <- 1
      x1 <- x0
      f1 <- f0
      #print(paste("retcode:", retcode))
      break

    } else if (f_evals >= max_f_evals) {
      retcode <- 2
      x1 <- x0
      f1 <- f0
      #print(paste("retcode:", retcode))
      break

    } else {

      if (is.na(alpha_prev)) {
        denom <- (f_trial - f0 - initslope)
        if (abs(denom) < 1e-14) {
          alpha_temp <- alpha * 0.5
          retcode <- 3
          #print(paste("retcode:", retcode))
        } else {
          alpha_temp <- -initslope / (2 * denom)
          retcode <- 4
          #print(paste("retcode:", retcode))
        }
      } else {
        d1 <- f_prev  - f0 - initslope * alpha_prev
        d2 <- f_trial - f0 - initslope * alpha

        denom <- (alpha - alpha_prev) * (alpha^2) * (alpha_prev^2)
        if (abs(denom) < 1e-14) {
          # fallback => half alpha
          alpha_temp <- alpha * 0.5
          retcode <- 5
          #print(paste("retcode:", retcode))
        } else {

          retcode <- 6
          #print(paste("retcode:", retcode))
          alphaC <- alpha_prev - ( (alpha^2 * d1) - (alpha_prev^2 * d2) ) / ( (alpha^2 - alpha_prev^2) * initslope * 2 )
          if (is.na(alphaC) || is.nan(alphaC) || alphaC <= 0) {
            alpha_temp <- 0.5 * alpha
            retcode <- 7
            #print(paste("retcode:", retcode))
          } else {
            alpha_temp <- alphaC
            retcode <- 8
            #print(paste("retcode:", retcode))
            #print(paste("Current alpha:", paste(alpha, collapse=", "), "Length:", length(alpha)))
          }
        }
      }

      if (alpha_temp > 0.5 * alpha) {
        alpha_temp <- 0.5 * alpha
        retcode <- 9
        #print(paste("retcode:", retcode))
        #print(paste("Current alpha:", paste(alpha, collapse=", "), "Length:", length(alpha)))
      } else if (alpha_temp <= 0) {
        alpha_temp <- 0.5 * alpha
        retcode <- 10

        #print(paste("retcode:", retcode))
      }

      alpha_prev <- alpha
      f_prev     <- f_trial

      # update alpha
      alpha <- alpha_temp
      #print(paste("alpha:", alpha))
    } 
  }  

  return(list(
    x1       = x1,
    f1       = f1,
    alpha    = alpha,
    retcode  = retcode,
    f_evals  = f_evals
  ))
}


optimize_helical <- function(f, x0, max_iter=600, tol=1e-6) {
    x <- x0
    iter <- 0
    f_best <- f(x)
    x_best <- x
    
    while (iter < max_iter) {
        g <- grad_fhelical(x)
        grad_norm <- sqrt(sum(g^2))
        
        if (!is.na(grad_norm) && !is.nan(grad_norm) && grad_norm < tol) {
            return(list(x=x_best, f=f_best, iter=iter, status="converged"))
        }
        
        p <- -g
        
        f_current <- f(x)
        result <- line_search(x, f_current, f, g, p)
        
        x_new <- result$x1
        f_new <- result$f1
        
        if (!is.na(f_new) && !is.nan(f_new) && f_new < f_best) {
            f_best <- f_new
            x_best <- x_new
        }
        
        step_size <- sqrt(sum((x_new - x)^2))
        if (!is.na(step_size) && !is.nan(step_size) && step_size < tol) {
            return(list(x=x_best, f=f_best, iter=iter, status="small_step"))
        }
        
        x <- x_new
        iter <- iter + 1
    }
    
    return(list(x=x_best, f=f_best, iter=iter, status="max_iter"))
}

x0 <- c(-1, 0, 0)    


result <- optimize_helical(Fhelical, x0)
distance <- sqrt(sum((result$x - c(1,0,0))^2))

print(result)
print(paste("Starting from:", paste(x0, collapse=" ")))
print(paste("Final point:", paste(result$x, collapse=" ")))
print(paste("Function value:", result$f))
print(paste("Distance to optimal:", distance))
print(paste("Status:", result$status))
    



    
