library(matlib)

#####
# 2 #
#####

#########################
# Zadatak 1.a) Kolokvijum
#########################

A = matrix(c(
  1, 0, -1,  2,
  2, 1, -2,  5,
  0, 1,  1,  0,
  3, 1, -3,  7,
  1, 2,  0,  1
), ncol = 4, byrow = TRUE)

svd_A = svd(A)

U = svd_A$u
Sigma = svd_A$d
V = svd_A$v

tol = 1e-10
r = sum(Sigma > tol)

# Slika

dim_Im = r
basis_Im = U[, 1:r]

# Kernel

n = ncol(A)
dim_Ker = n - r

if(dim_Ker > 0){
  basis_Ker = V[, (r+1):n]
} else {
  basis_Ker = 0
}

# Rezultati

dim_Im
basis_Im

dim_Ker
basis_Ker


#########################
# Zadatak 1.b) Kolokvijum
#########################

# dati vektor y (iz kolokvijuma)
y = matrix(c(1, 2, 0, 3, 1), ncol = 1)

#  SVD iz prethodnog dela
U_r = U[, 1:r]

# Projekcija na sliku

P = U_r %*% t(U_r)

y_proj = P %*% y

# Provera

round(y - y_proj, 10)

# y jeste u slici

#########################
# Zadatak 1.c) Kolokvijum
#########################

# SVD već imamo:
# U, Sigma (svd_A$d), V


# Konstrukcija Sigma^+

Sigma_inv = diag(1 / Sigma[1:r])

# Pseudo-inverz

A_pinv = V[,1:r] %*% Sigma_inv %*% t(U[,1:r])

# Posebno rešenje

x_resenje = A_pinv %*% y

x_resenje


# Provera

round(A %*% x_resenje - y, 10)

#########################
# Zadatak 1.d) Kolokvijum
#########################

if(dim_Ker == 0){
  x_resenje
} else {
  list(
    posebno_resenje = x_resenje,
    baza_kernela = basis_Ker
  )
}


########################
# Zadatak 2 - KOlokvijum
########################

# Definicija vektora
u1 = c(0, -1, 2, 0, 2)
u2 = c(1, -3, 1, -1, 2)
u3 = c(-3, 4, 1, 2, 1)
u4 = c(-1, -3, 5, 0, 7)

x = c(-9, -1, -1, 4, 1)

# Formiramo matricu A

A = cbind(u1, u2, u3, u4)

# SVD

svd_A = svd(A)

U = svd_A$u
Sigma = svd_A$d

tol = 1e-10
r = sum(Sigma > tol)

# Ortonormalna baza potprostora

U_r = U[,1:r]

# (a) Projekcija

x_proj = U_r %*% t(U_r) %*% x
x_proj

# (b) Rastojanje

dist = sqrt( sum( (x - x_proj)^2 ) )
dist


########################
# Zadatak 3 - Kolokvijum
########################

A = matrix(c(4, 2, 3,
             -3, -1, -3,
             -2, -2, -1),
           ncol = 3, byrow = TRUE)

svd_A = svd(A)

U = svd_A$u
Sigma_val = svd_A$d
V = svd_A$v

Sigma = diag(Sigma_val)

U
Sigma
V

