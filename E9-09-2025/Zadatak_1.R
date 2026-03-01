library(matlib)

A = matrix(c(
  1,  3,  1, 2, 1,  6,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 2, 1,  6,
  2,  6,  2, 4, 2, 12),
  ncol = 6, byrow = TRUE)
print(A)



#SVD f-ja

SVD = function(A, tol = 1e-8){
  
  m = nrow(A)
  n = ncol(A)
  
  AtA = t(A) %*% A
  eig = eigen(AtA)
  
  V = eig$vectors
  lambda = eig$values
  
  sigma = sqrt(pmax(lambda, 0))
  
  #Sortiranje indeksa
  indeks = order(sigma, decreasing = TRUE)
  sigma = sigma[indeks]
  V = V[, indeks]
  
  #Rang
  r = sum(sigma > tol)
  
  sigma_r = sigma[1:r]
  V_r = V[, 1:r]
  
  

  
  
  U_r = matrix(0, m, r)
  
  for(i in 1:r){
    U_r[, i] = (A %*% V_r[, i]) / sigma_r[i]
  }
  
  Sigma_r = diag(sigma_r)
  
  list(U = U_r, Sigma = Sigma_r, V = V_r, d = sigma_r)
}


res = SVD(A)
print(res)


#Provera

A_reconstructed = res$U %*% res$Sigma %*% t(res$V)

round(A - A_reconstructed, 6)


#Ugradjen SVD

builtin = svd(A)

builtin$d
res$d
