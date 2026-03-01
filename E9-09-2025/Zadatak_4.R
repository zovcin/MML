#Zadatak 5.1

library(matlib)                         

n <- 4                                  

A <- matrix(c(1,1,1,1,                  
              0,1,1,1,
              0,0,1,1,
              0,0,0,1),
            ncol = n, byrow = TRUE)     
print(A)



#F-ja
f <- function(x) {                      
  y <- A %*% x   
  print(y)
  term1 <- t(x) %*% t(A) %*% A %*% x    
  term2 <- 0.1 * sum(y^3)               
  term3 <- 0.01 * sum(y^4)              
  return(as.numeric(term1 + term2 + term3)) #as.numeric ce vratiti broj
}



#Gradijent
grad_f <- function(x) {                 
  y <- A %*% x                          
  print(y)
  term1 <- 2 * t(A) %*% A %*% x         
  term2 <- 0.3 * t(A) %*% (y^2)         
  term3 <- 0.04 * t(A) %*% (y^3)        
  return(term1 + term2 + term3)         
}



gamma <- 0.01                           
n_iter <- 100                           

x <- 0.2 * matrix(1, n, 1)              




for (k in 1:n_iter) {                 
  x <- x - gamma * grad_f(x)          
}



#Rezultat
x                                       
f(x)                                    

