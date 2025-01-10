ef <- function(x1,x2,b){
  y <- exp( 0.1 * sqrt( x1^2 + b*x2^2 ) ) + exp( -1*0.1 * sqrt( x1^2 + b*x2^2 ) ) - 2
  return(y)
}
efMinusB <- function(x1,x2,b){
  y <- exp( 0.1 * sqrt( x1^2 + b*x2^2 ) ) - exp( -1*0.1 * sqrt( x1^2 + b*x2^2 ) )
  return(y)
}
efPlusB <- function(x1,x2,b){
  y <- exp( 0.1 * sqrt( x1^2 + b*x2^2 ) ) + exp( -1*0.1 * sqrt( x1^2 + b*x2^2 ) )
  return(y)
}
f <- function(x1,x2,b){
  y <- 200* ef(x1,x2,b) / ef(10,10,b)
  return(y)
}
pravac <- function(x1,x2,b){
  dfdx1 <- 200/ef(10,10,b) * 0.1*x1*efMinusB(x1,x2,b)/sqrt(x1*x1+b*x2*x2)
  dfdx2 <- 200/ef(10,10,b) * 0.1*b*x2*efMinusB(x1,x2,b)/sqrt(x1*x1+b*x2*x2)
  
  dfdx1dx1 <- 200/ef(10,10,b)*0.1*( -1/sqrt( (x1*x1+b*x2*x2)^3 )*x1*x1*efMinusB(x1,x2,b) + 1/sqrt(x1*x1+b*x2*x2)*(efMinusB(x1,x2,b)+0.1/sqrt(x1*x1+b*x2*x2)*x1*x1*efPlusB(x1,x2,b)))
  dfdx1dx2 <- 200/ef(10,10,b)*0.1*x1*( -1/sqrt((x1*x1+b*x2*x2)^3)*b*x2*efMinusB(x1,x2,b)+0.1*b/(x1*x1+b*x2*x2)*x2*efPlusB(x1,x2,b))
  dfdx2dx2 <- 200/ef(10,10,b)*0.1*b*(-1*b/sqrt((x1*x1+b*x2*x2)^3)*x2*x2*efMinusB(x1,x2,b)+1/sqrt(x1*x1+b*x2*x2)*(efMinusB(x1,x2,b)+0.1*b/sqrt(x1*x1+b*x2*x2)*x2*x2*efPlusB(x1,x2,b)))
  dfdx2dx1 <- 200/ef(10,10,b)*0.1*b*x2*(-1/sqrt((x1*x1+b*x2*x2)^3)*x1*efMinusB(x1,x2,b)+1/(x1*x1+b*x2*x2)*0.1*efPlusB(x1,x2,b))
  
  
  grad2 <- matrix(c(dfdx1dx1,dfdx1dx2,dfdx2dx1,dfdx2dx2),nrow = 2, ncol = 2, byrow = TRUE)
  grad1 <- matrix(c(dfdx1, dfdx2), nrow = 2, ncol = 1)
  
  grad <- -1*(solve(grad2)%*%grad1)
  
  return(grad)
}

izvod <- function(x1,x2,b){
  dfdx1 <- 200/ef(10,10,b) * 0.1*x1*efMinusB(x1,x2,b)/sqrt(x1*x1+b*x2*x2)
  dfdx2 <- 200/ef(10,10,b) * 0.1*b*x2*efMinusB(x1,x2,b)/sqrt(x1*x1+b*x2*x2)
  grad1 <- matrix(c(dfdx1, dfdx2), nrow = 2, ncol = 1)
  return(grad1)
}

backtrackingLineSearch <- function(x1,x2,b,p){
  alpha <- list(1,1)
  c <- 0.0001
  grad <- izvod(x1,x2,b)
  y0 <- f(x1+alpha[[1]]*p[1,1],x2+alpha[[2]]*p[2,1],b)
  y1 <- f(x1,x2,b)+c*( t(grad)%*%p)
  while(y0>y1){
    alpha = alpha%*%p
    
    y0 <- f(x1+alpha[[1]]*p[1,1],x2+alpha[[2]]*p[2,1],b)
    y1 <- f(x1,x2,b)+c*( t(grad)%*%p)
  }
  return(alpha)
}

pronalazakMinimuma <- function(x1,x2,b){
  NUMBER_OF_ITERATIONS <- 4
  i <- 0
  while(i<NUMBER_OF_ITERATIONS){
    p <- pravac(x1,x2,b)
    alpha = backtrackingLineSearch(x1,x2,b,p)
    x1 <- x1+alpha[[1]]*p[1,1]
    x2 <- x2+alpha[[2]]*p[2,1]
    i <- i+1
  }
  return(list(x1,x2))
}

x1 <- 0.5
x2 <- 0.5
b <- 1
x<-pronalazakMinimuma(x1,x2,b)
print(x)
