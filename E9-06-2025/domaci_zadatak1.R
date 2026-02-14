library(matlib)
my_svd <- function(A) {
  A <- as.matrix(A)
  
  ATA <- t(A) %*% A
  eg <- eigen(ATA, symmetric = TRUE)
  
  ord <- order(eg$values, decreasing = TRUE)
  lambda <- eg$values[ord]
  V <- eg$vectors[, ord, drop = FALSE]
  
  lambda[lambda < 0] <- 0
  d <- sqrt(lambda)
  
  r <- sum(d > 0)
  
  Vr <- V[, 1:r, drop = FALSE]
  dr <- d[1:r]
  
  U <- A %*% Vr
  U <- sweep(U, 2, dr, "/")
  
  list(d = dr, u = U, v = Vr)
}


#-------------------------
#TESTIRANJE

set.seed(1)

#primer m x n (5x3 matrica)
A <- matrix(rnorm(5*3), 5, 3)

cat("Ulazna matrica A:\n")
print(round(A, 4))

#poziv funkcije
res <- my_svd(A)

cat("\nSingularne vrednosti (moja funkcija):\n")
print(round(res$d, 6))

#rekonstrukcija
A_hat <- res$u %*% diag(res$d) %*% t(res$v)

cat("\nRekonstruisana matrica A_hat:\n")
print(round(A_hat, 4))

#greška rekonstrukcije
cat("\nReconstruction error (Frobenius norm):\n")
print(norm(A - A_hat, "F"))

#poređenje sa ugrađenom svd funkcijom
d_base <- svd(A)$d

cat("\nSingularne vrednosti (R svd):\n")
print(round(d_base, 6))

cat("\nMaksimalna razlika singularnih vrednosti:\n")
print(max(abs(res$d - d_base)))

cat("\nU^T U, provera da je U ortonormalna:\n")
print(round(t(res$u) %*% res$u, 6))

cat("\nV^T V, provera da je U ortonormalna:\n")
print(round(t(res$v) %*% res$v, 6))
