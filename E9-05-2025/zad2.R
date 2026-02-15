library(matlib)

# 1
A = matrix(c(
  2, -1,  0,  1,  3,
  3,  1,  2, -1,  0,
  0,  2,  1, -4, -2,
  1,  0,  1,  2, -1
), nrow = 5, ncol = 4)

y = matrix(c(-8, 3, -5, -17, 1), ncol = 1)

s = svd(A, nu = nrow(A), nv = ncol(A))
U = s$u; d = s$d; V = s$v
rank_a = sum(zapsmall(d) > 0)

d; rank_a

# 1a
# dim Im F = rank_a = 3 

# dim Ker F = 4 - 3 = 1, baza iz poslednje kolone V
kernel_A = V[, 4]

print(all(round(A %*% kernel_A, digits = 13) == 0))

# iz kernel_A: -0.5*a1 + 0.5*a2 - 0.5*a3 - 0.5*a4 = 0 -> a4 = -a1 + a2 - a3
# baza Im F = prve tri kolone:
A[, 1:3]

# 1b
# A %*% x = y => t(A) %*% A %*% x = t(A) %*% y => x = inv(t(A) %*% A) %*% t(A) %*% y
# proj = A %*% x = A %*% inv(t(A) %*% A) %*% t(A) %*% y
# A = U %*% S %*% t(V)
# t(V) %*% V = I, S %*% inv(S) = I, sve se skrati:
# proj = U_R %*% t(U_R) %*% y
proj_y = U[, 1:rank_a] %*% t(U[, 1:rank_a]) %*% y
round(sqrt(t(y - proj_y) %*% (y - proj_y)), 13)  # = 0 => proj_y == y => y pripada Im F

# y = x1*a1 + x2*a2 + x3*a3, koordinate y u bazi Im F
img_F = A[, 1:3]
svd_img_F = svd(img_F, nu = nrow(img_F), nv = ncol(img_F))
u_img_F = svd_img_F$u[, 1:rank_a]
u_img_F_inv = t(u_img_F)
v_img_F = svd_img_F$v[, 1:rank_a]
sigma_inv_img_F = diag(1/svd_img_F$d[1:rank_a])
# img_A %*% x = y
# U %*% S %*% t(V) %*% x = y
# x = inv(t(V)) %*% inv(S) %*% inv(U) %*% y   (inv(t(V)) = V, inv(S) = diag(1/S), inv(U) = t(U))
# x = V * diag(1/S) * t(U) * y
x = v_img_F %*% sigma_inv_img_F %*% u_img_F_inv %*% y; x
all(zapsmall(img_F %*% x) == y)  # Kordinate y u bazi Im F su x = [5, -6, 7] -> y = 5*a1 - 6*a2 + 7*a3

# 1c
# x = [5, -6, 7] za kolone 1,2,3; kolona 4 zavisna => 0
x1 = c(x[1], x[2], x[3], 0)
x1; zapsmall(A %*% x1)  # x1 = [5, -6, 7, 0]

# x2 = x1 + L * kernel_A (dim Ker = 1) za svako L e R
x2 = x1 + kernel_A
x2; zapsmall(A %*% x2)

# 1d
# razlika bilo koja dva resenja = lambda * kernel_A (1D)
# ne postoji x3 linearno nezavisno od x2 - x1


# 2
A = matrix(c(0,-1,2,0,2,1,-3,1,-1,2,-3,4,1,2,1,-1,-3,5,0,7), ncol = 4); A
x = matrix(c(-9,-1,-1,4,1), ncol = 1)

s2 = svd(A, nu = nrow(A), nv = ncol(A))
rank_a = R(A)
s2$d; rank_a

# 2a
# A %*% x_proj = x => t(A) %*% A %*% x_proj = t(A) %*% x => x_proj = inv(t(A) %*% A) %*% t(A) %*% x
# A = U %*% S %*% t(V) (SVD)
# x_proj = inv(V %*% t(S) %*% t(U) %*% U %*% S %*% t(V)) %*% V %*% t(S) %*% t(U) %*% x
#        = inv(V %*% t(S) %*% S %*% t(V)) %*% V %*% t(S) %*% t(U) %*% x       (t(U) %*% U = I)
#        = V %*% inv(S) %*% inv(t(S)) %*% t(V) %*% V %*% t(S) %*% t(U) %*% x
#        = V %*% inv(S) %*% t(U) %*% x
# proj = A %*% x_proj = U %*% S %*% t(V) %*% V %*% inv(S) %*% t(U) %*% x = U_R %*% t(U_R) %*% x
u = s2$u[, 1:rank_a]
v = s2$v[, 1:rank_a]
sigma_inv = diag(1/s2$d[1:rank_a])

x_proj = v %*% sigma_inv %*% t(u) %*% x; x_proj # koordinate projekcije x u bazi A = [-7.18, 3.14, 2.28, 1.38]

proj = A %*% x_proj; proj # ortogonalna projekcija x na A = [-5.10, 2.75, -2.02, 1.43, 3.89]

all(round(matrix(Proj(x, A)),13) == round(proj, 13)) # dobije se isto sto i sa funckijom Proj

# 2b
d = x - proj
sqrt(t(d) %*% d)  # d(x, A) = 6.728287


# 3
A = matrix(c(4,2,3, -3,-1,-3, -2,-2,-1), nrow = 3)
v1 = c(1,0,1); v2 = c(1,1,0); v3 = c(1,1,1)

# 3a
A %*% v1  # = 2*v1 => lambda1 = 2
A %*% v2  # = 1*v2 => lambda2 = 1
A %*% v3  # = -1*v3 => lambda3 = -1

# 3b
P = cbind(v1, v2, v3)
svd_P = svd(P)
sum(zapsmall(svd_P$d) > 0)  # = 3 => P punog ranga => A3 dijagonalizabilna
D = diag(c(2, 1, -1))
P_inv = inv(P); P_inv

all((P %*% D %*% P_inv) == A) # dobija se matrica A
