library(matlib)


#####
# 4 # 
#####

n = 4

# Gornja trougaona matrica sa jedinicama
A = matrix(c(1,1,1,1,
             0,1,1,1,
             0,0,1,1,
             0,0,0,1),
           ncol = n, byrow = TRUE)

# Funkcija cilja

f = function(x){
  y = A %*% x
  term1 = t(x) %*% t(A) %*% A %*% x
  term2 = 0.1 * sum(y^3)
  term3 = 0.01 * sum(y^4)
  return(term1 + term2 + term3)
}


# Gradijent

grad_f = function(x){
  y = A %*% x
  
  term1 = 2 * t(A) %*% A %*% x
  term2 = 0.3 * t(A) %*% (y^2)
  term3 = 0.04 * t(A) %*% (y^3)
  
  return(term1 + term2 + term3)
}

# Parametri

gamma = 0.01
n_iter = 100

# Početna tačka
x = 0.2 * matrix(1, n, 1)


# Gradient Descent

for(k in 1:n_iter){
  x = x - gamma * grad_f(x)
}

# Rezultat

x
f(x)
