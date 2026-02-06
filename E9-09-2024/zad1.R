library(matlib)

svd_fun <- function(A) {
  m <- nrow(A)
  n <- ncol(A)
  
  AtA <- t(A) %*% A
  
  eg <- eigen(AtA)
  
  V = GramSchmidt(eg$vectors, normalize = TRUE)
  
  singular_values <- sqrt(pmax(eg$values, 0))

  U_tmp = matrix(0, m, m)
  for (i in 1:n) {
    if (singular_values[i] > 1e-30) {
      U_tmp[, i] <- (A %*% V[, i]) / singular_values[i]
    }
  }
  
  U_tmp = cbind(U_tmp, diag(m))
  
  U = GramSchmidt(U_tmp, normalize = TRUE)
  
  return(list(
    d = singular_values,
    u = U,
    v = V
  ))
}

A=matrix(c(
  1,  3,  1, 2, 1,  6,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 2, 1,  6,
  2,  6,  2, 4, 2, 12),ncol=6,byrow=T)

svd_test <- svd_fun(A)
svd_test

r_svd <- svd(A, nu=nrow(A), nv=ncol(A))
r_svd

# kada se zaokruze vrednosti matrice B, kako bi se eliminisale greske nastale 
# prilikom izracunavanje eigen vektora i vrenosti, i porede sa originalnim   
# vrednostima matrice A mozemo videti da razlika nema, odnosno da je SVD tacno izracunat

Sigma <- matrix(0, nrow(A), ncol(A))
diag(Sigma) <- svd_test$d

# Now the dimensions match: (7x7) %*% (7x6) %*% (6x6)
B = round(svd_test$u %*% Sigma %*% t(svd_test$v))
B == A

# Pri poredjenju singularnih vrednosti mozemo videti da su razlike minimalne, a 
# uzrok samih razlika je kompleksnost funkcije eigen, koju R-ov SVD 
# verovatno nije koristio

svd_test$d - r_svd$d




