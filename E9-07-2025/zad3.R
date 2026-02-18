library(matlib)

#####
# 3 # 
#####

# Gradient Descent (100 iteracija)
# Matrica Q i vektor b
Q = matrix(c(2, 1,
             1, 20), 
           ncol = 2, byrow = TRUE)

b = matrix(c(5, 3), ncol = 1)

# Funkcija f(x)

f = function(x){
  return( 0.5 * t(x) %*% Q %*% x - t(b) %*% x )
}

# Gradijent

grad_f = function(x){
  return( Q %*% x - b )
}

# Parametri

gamma = 0.085
n_iter = 100

# Pcetna tacka
x = matrix(c(-3, -1), ncol = 1)

# Cuvanje iteracija
X_iteracije = matrix(0, nrow = 2, ncol = n_iter + 1)
X_iteracije[,1] = x

# Gradient Descent petlja

for(k in 1:n_iter){
  x = x - gamma * grad_f(x)
  X_iteracije[,k+1] = x
}

# Rezultat nakon 100 iteracija

x
f(x)
