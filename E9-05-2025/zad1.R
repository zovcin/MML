library(matlib)

bole_SVD = function(A) {
  AtA = t(A) %*% A
  eig = eigen(AtA)
  rank_A <- R(A)

  lambda = eig$values
  lambda[(rank_A + 1):length(lambda)] = 0
  d = sqrt(lambda)

  V = GramSchmidt(eig$vectors, normalize = TRUE)
  V = V[, 1:rank_A]

  # u_i = A%*%v_i / sigma_i
  U = matrix(0, nrow(A), rank_A)
  for (i in 1:rank_A) {
      U[, i] = (A %*% V[, i]) / d[i]
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
u = result$u
d = diag(result$d)[1:R(A),1:R(A)]
vt = t(result$v)
cat("U:", dim(u), "\nD:", dim(d), "\nVt:", dim(vt), "\n")
A2 = u %*% d %*% vt
A_z <- zapsmall(A); A_z
A2_z <- zapsmall(A2); A2_z
all(A_z == A2_z)
