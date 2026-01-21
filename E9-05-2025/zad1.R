library(matlib)

bole_SVD = function(A) {
  AtA = t(A) %*% A
  eig = eigen(AtA)

  lambda = eig$values
  d = sqrt(lambda)
  d[is.nan(d)] = 0

  V = GramSchmidt(eig$vectors, normalize = TRUE)

  # u_i = A%*%v_i / sigma_i
  U = matrix(0, nrow(A), ncol(A))
  for (i in 1:ncol(A)) {
    if (d[i] > sqrt(.Machine$double.eps)) {
      U[, i] = (A %*% V[, i]) / d[i]
    }
  }

  list(d = d, u = U, v = V)
}

A = matrix(c(
  1,  3,  1, 2, 1,  6,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 2, 1,  6,
  2,  6,  2, 4, 2, 12), ncol = 6, byrow = TRUE)

print(A)

result <- bole_SVD(A)

# poredjenje singularnih vrednosti
result_d <- zapsmall(result$d); result_d
svd_d <- zapsmall(svd(A)$d); svd_d
all(result_d == svd_d)

# poredjenje matrica
A2 = result$u %*% diag(result$d) %*% t(result$v)
A_z <- zapsmall(A); A_z
A2_z <- zapsmall(A2); A2_z
all(A_z == A2_z)
