b = 1
x0 = c(0.5, 0.5)
phi_b = function(x, b_val = 1) {
  norm_sq = x[1]^2 + b*x[2]^2
  if(norm_sq == 0) return(0)
  u = 0.1 * sqrt(norm_sq)
  return(exp(u) + exp(-u) - 2)
}
c = 200/phi_b(c(10,10))

f=function(x) {
  return(c * phi_b(x))
}

grad_f=function(x) {
  x1 = x[1]
  x2 = x[2]
  if (x1 == 0 && x2 == 0) {
    return(c(0,0))
  }
  u = sqrt(x1^2 + x2^2)
  common=c*0.1*(exp(0.1*u) - exp(-0.1*u))
  dx1 = common * (x1/u)
  dx2 = common * (x2/u)
  return(c(dx1,dx2))
}

line_search=function(x,p,g,c1=1e-4,c2=0.9) {
  alpha = 1.0  
  alpha_lo = 0  
  alpha_hi = Inf  
  max_iter = 20
  
  f_x=f(x)
  slope0 = sum(g*p)  
  for (i in 1:max_iter) {
    x_new = x + alpha*p
    f_new = f(x_new)
    
    if((f_new > f_x+c1*alpha*slope0) || (i > 1 && f_new >= f_lo)) {
      return(zoom(x,p,g,alpha_lo,alpha))
    }
    g_new=grad_f(x_new)
    slope_new=sum(g_new*p)
    
    if (abs(slope_new) <= -c2*slope0) {
      return(alpha)
    }
    
    if (slope_new >= 0) {
      return(zoom(x,p,g,alpha,alpha_lo))
    }
    alpha_lo=alpha
    f_lo=f_new
    alpha=alpha*2 
  }
  
  return(alpha)
}

zoom = function(x,p,g,alpha_lo,alpha_hi,max_iter=20) {
  f_x = f(x)
  slope0 = sum(g*p)
  c1 = 1e-4
  c2 = 0.9
  
  for (i in 1:max_iter) {
    alpha = (alpha_lo + alpha_hi) / 2
    x_new = x + alpha*p
    f_new = f(x_new)
    
    if (f_new > f_x + c1*alpha*slope0 || f_new >= f(x+alpha_lo*p)) {
      alpha_hi=alpha 
    } else {
      g_new = grad_f(x_new)
      slope_new = sum(g_new*p)
      
      if (abs(slope_new) <= -c2*slope0) {
        return(alpha)  #uspeh
      }
      
      if (slope_new * (alpha_hi - alpha_lo) >= 0) {
        alpha_hi = alpha_lo
      }
      alpha_lo = alpha
    }
  }
  return(alpha)
}

quadappr=function(x0,max_iter=100,tolerance=1e-6) {
  x_current=x0
  final_ter = max_iter 
  
  for (i in 1:max_iter) {
    g_current = grad_f(x_current)
    grad_norm = sqrt(sum(g_current^2))
    
    if (grad_norm < tolerance) {
      final_iter = i
      break
    }
    
    p_current = -g_current
    
    alpha = line_search(x_current, p_current, g_current)
    
    x_current = x_current + alpha*p_current
  }
  
  return(list(
    optimal_x = x_current,
    iterations = final_iter,
    f_val = f(x_current)
  ))
}

result = quadappr(x0)
cat(sprintf("Konvergencija postignuta u %d iteraciji\n",result$iterations))
cat("Optimalna tačka (x1, x2):", sprintf("(%.6e, %.6e)", result$optimal_x[1], result$optimal_x[2]), "\n")
cat("Vrednost funkcije u optimumu f(x):", sprintf("%.6e", result$f_val), "\n")