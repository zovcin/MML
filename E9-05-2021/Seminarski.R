euclidean <- function(a, b) sqrt(sum((a - b)^2))

optimm <- function(x, fn, nablafn, j, deltax) {
  p = -nablafn(x)
  
  for (i in 1:j) {
    alfa <- bls(fn,nablafn,x,p)
    prev <- x; prev
    x <- x + alfa * p; x
    p <- -nablafn(x)
    if (euclidean(x,prev) < deltax) {
      break
    } 
  }
  print(i)
  return(x)
}

bls <- function(fn,nablafn=nablafn,x,p) {
  alfa = 1
  plo = 1e-4
  phi = 0.99
  c=0.0001
  ro=0.5
  while(fn(x + alfa * p) > fn(x) + c * alfa * t(nablafn(x)) %*% p){
    alfa = ro * alfa
    if (fn(x + alfa * p) > fn(x) + c * alfa * t(nablafn(x)) %*% p) {
      ro <- max(plo, ro * 0.9)  
    } else {
      ro <- min(phi, ro * 1.1)
    }
  }
  return(alfa)
}

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
  return(as.numeric(y))
}

nablafn <- function(x) {
  g=as.vector(c(
    (2*x[1]+2*x[2]+2*x[3]+2*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3),
    (2*x[1]+4*x[2]+4*x[3]+4*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(3/10)*(x[2]+x[3]+x[4])^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3+(1/25)*(x[2]+x[3]+x[4])^3),
    (2*x[1]+4*x[2]+6*x[3]+6*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(3/10)*(x[2]+x[3]+x[4])^2+(3/10)*(x[3]+x[4])^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3+(1/25)*(x[2]+x[3]+x[4])^3+(1/25)*(x[3]+x[4])^3),
    (2*x[1]+4*x[2]+6*x[3]+8*x[4]+(3/10)*(x[1]+x[2]+x[3]+x[4])^2+(3/10)*(x[2]+x[3]+x[4])^2+(3/10)*(x[3]+x[4])^2+(3/10)*x[4]^2+(1/25)*(x[1]+x[2]+x[3]+x[4])^3+(1/25)*(x[2]+x[3]+x[4])^3+(1/25)*(x[3]+x[4])^3+(1/25)*x[4]^3)
  ))
  return(matrix(g,ncol=1))
}

x0=c(0.2,0.2,0.2,0.2);
nula <- optimm(x0, fn=fn, nablafn=nablafn, j=1000, deltax=1e-7)
nula
fn(nula)
