library('matlib')
n=4
A=matrix(c(1,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1), ncol = n,byrow=TRUE);A
x0=c(0.2,0.2,0.2,0.2);x0
f <- function(x) {
  y = t(x)%*%t(A)%*%A%*%x + 0.1*sum(sapply(A%*%x,function(v) v^3)) + 0.01*sum(sapply(A%*%x,function(v) v^4))
  return(y)
}
nablaf <- function(x) {
  nablay = 2*t(A)%*%A%*%x + 0.3*t(A)%*%sapply(A%*%x, function(v) v^2) + 0.04*t(A)%*%sapply(A%*%x, function(v) v^3)
  return(nablay)
}

alfa = 1
plo = 1e-4
phi = 0.99
x = x0
c=0.0001

bls <- function(f,grad_f,x,p) {
  ro=0.5
  while(f(x + alfa * p) > f(x) + c * alfa * t(gradf_x) %*% p){
    alfa = ro * alfa
    if (f(x + alfa * p) > f(x) + c * alfa * t(gradf_x) %*% p) {
      ro <- max(plo, ro * 0.9)  
    } else {
      ro <- min(phi, ro * 1.1)
    }
    
  }
  return(alfa)
}
euclidean <- function(a, b) sqrt(sum((a - b)^2))
grad <- function(x = x0, j = 100, deltax=1e-7) {
  p = - nablaf(x)
  
  for (i in 1:j) {
    alfa <- bls(f(x),nablaf(x),x,p)
    prev <- x; prev
    x <- x + alfa * p; x
    p <- -nablaf(x)
    if (euclidean(x,prev) < deltax) {
      i
      break
    } 
  }
  return(x)
}

grad() -> nula
round(nula)
round(f(nula))
