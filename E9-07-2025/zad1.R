library(matlib)

#####
# 1 #
#####


# 1) Zadata matrica 

A = matrix(c(
  1,  3,  1, 2, 1,  6,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 2, 1,  6,
  2,  6,  2, 4, 2, 12),
  ncol = 6, byrow = TRUE)

# 2) SVD funkcija

mySVD = function(A, tol = 1e-8){
  
  m = nrow(A)
  n = ncol(A)
  
  AtA = t(A) %*% A
  eig = eigen(AtA)
  
  V = eig$vectors
  lambda = eig$values
  
  sigma = sqrt(pmax(lambda, 0))
  
  #sortirani indeksi
  indeksi = order(sigma, decreasing = TRUE)
  sigma = sigma[indeksi]
  V = V[, indeksi]
  
  # rang
  r = sum(sigma > tol)
  
  sigma_r = sigma[1:r]
  V_r = V[, 1:r]
  
  # primena Gram-Schmidt na V
  V_r = GramSchmidt(V_r)
  
  # računanje U
  U_r = matrix(0, m, r)
  
  for(i in 1:r){
    U_r[, i] = (A %*% V_r[, i]) / sigma_r[i]
  }
  
  Sigma_r = diag(sigma_r)
  
  list(U = U_r, Sigma = Sigma_r, V = V_r, d = sigma_r)
}


res = mySVD(A)
print(res)

# 4) Provera (rezultat treba da bude 0)

A_reconstructed = res$U %*% res$Sigma %*% t(res$V)

round(A - A_reconstructed, 6)

# 5) Poređenje sa ugrađenim SVD


builtin = svd(A)

builtin$d
res$d
