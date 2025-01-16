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

ef <- function(x){
  y <- exp( 0.1 * sqrt( x[1]^2 + x[3]*x[2]^2 ) ) + exp( -1*0.1 * sqrt( x[1]^2 + x[3]*x[2]^2 ) ) - 2
  return(y)
}

fn <- function(x) {
  b <- 1
  x <- c(x, b)
  
  y <- y <- 200* ef(x) / ef(c(10,10,b))
  return(as.numeric(y))
}

efMinusB <- function(x){
  y <- exp( 0.1 * sqrt( x[1]^2 + x[3]*x[2]^2 ) ) - exp( -1*0.1 * sqrt( x[1]^2 + x[3]*x[2]^2 ) )
  return(y)
}

nablafn <- function(x) {
  b <- 1
  x <- c(x, b)
  dfdx1 <- 200/ef(c(10,10,b)) * 0.1*x[1]*efMinusB(x)/sqrt(x[1]^2+x[3]*x[2]^2)
  dfdx2 <- 200/ef(c(10,10,b)) * 0.1*b*x[2]*efMinusB(x)/sqrt(x[1]^2+x[3]*x[2]^2)
  g=as.vector(c(dfdx1, dfdx2))
  return(matrix(g,ncol=1))
}

x0=c(0.5,0.5);
nula <- optimm(x0, fn=fn, nablafn=nablafn, j=1000, deltax=1e-7)
nula
fn(nula)
