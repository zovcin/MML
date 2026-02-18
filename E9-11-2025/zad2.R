library(matlib)

svd_manual = function(A) {
  
  # Pretvaranje u matricu (za sigurnost)
  A = as.matrix(A)
  
  # Dimenzije matrice A
  m = nrow(A)   # broj vrsta
  n = ncol(A)   # broj kolona
  
  tol = 1e-8    # prag numeričke nule
  
  # 1. Formiramo matricu A^T A
  # Ovo je simetrična matrica dimenzije n x n
  # Njena eigen dekompozicija daje DESNE singularne vektore (V)
  # Ista funkcija je mogla biti kreirana i prvo trazeci LEVE singularne vektore (U)
  AtA = t(A) %*% A
  
  # Eigen dekompozicija
  eig = eigen(AtA)
  
  # Sopstvene vrednosti 
  lambda = eig$values
  
  # Korekcija numeričkih grešaka za svaki slucaj
  lambda[lambda < 0] = 0
  
  # Singularne vrednosti
  sigma = sqrt(lambda)
  
  
  # Sortiranje singularnih vrednosti (za svaki slucaj, verovatno nije neophodno)
  ord = order(sigma, decreasing = TRUE)
  
  sigma = sigma[ord]                # sortirane singularne vrednosti
  V = eig$vectors[, ord]            # odgovarajući eigen vektori
  
  # Ortonormalizacija baze 
  V = GramSchmidt(V, normalize = TRUE)
  
  
  # Rang matrice
  r = sum(sigma > tol)
  
  # 2. Računanje leve matrice U
  # Koristimo relaciju:
  # u_i = (A v_i) / sigma_i
  U = matrix(0, m, m)   
  
  for (i in 1:r) {
    
    # Direktna formula iz teorije SVD-a
    U[, i] = (A %*% V[, i]) / sigma[i]
  }
  
  # Dopuna ortonormirane baze ako je rang < m
  if (r < m) {
    
    # Dodajemo standardnu bazu kao potencijalne clanove
    c = cbind(U[, 1:r, drop = FALSE], diag(m))
    
    # Gram–Schmidt pravi kompletnu ortonormiranu bazu
    Q = GramSchmidt(c, normalize = TRUE)
    
    # Uzimamo prvih m ortonormiranih vektora
    U = Q[, 1:m]
  }
  
  
  # 3. Dijagonalna matrica singularnih vrednosti
  D = matrix(0, m, n)
  
  for (i in 1:min(m, n)) {
    
    # Singularne vrednosti na dijagonalu
    D[i, i] = sigma[i]
  }
  
  
  # Povratne vrednosti
  list(
    d = sigma[1:min(m, n)],  # singularne vrednosti
    u = U,                   # levi singularni vektori
    v = V,                   # desni singularni vektori
    D = D
  )
}


# 1. ZADATAK — Linearna transformacija 


A1 <- matrix(c(
  2, 3, 0, 1,
  -1, 1, 2, 0,
  0, 2, 1, 1,
  1,-1,-4, 2,
  3, 0,-2,-1
), nrow = 5, byrow = TRUE)

y1 <- matrix(c(-8, 3, -5, -17, 1), ncol = 1)

svd1 <- svd_manual(A1)

U1 <- svd1$u
V1 <- svd1$v
d1 <- svd1$d

tol <- 1e-8
r1 <- sum(d1 > tol)

cat("Rang matrice A1:", r1, "\n\n")

# ----- Im -----
ImPhi <- U1[, 1:r1, drop = FALSE]
cat("Baza Im :\n")
print(ImPhi)
cat("\n")

# ----- Ker -----
KerPhi <- V1[, (r1+1):ncol(V1), drop = FALSE]
cat("Baza Ker:\n")
print(KerPhi)
cat("\n")

# ----- Provera -----
y_proj <- ImPhi %*% t(ImPhi) %*% y1
cat("||y - proj|| =", norm(y1 - y_proj, type = "2"), "\n\n")

# ----- Rešavanje Ax = y preko pseudoinverza -----
Sigma_inv <- diag(1 / d1[1:r1])
A1_pinv <- V1[, 1:r1] %*% Sigma_inv %*% t(U1[, 1:r1])

x1 <- A1_pinv %*% y1
cat("Jedno rešenje x1:\n")
print(x1)

B <- A1[, 1:3]    

# SVD baza
svd_res <- svd_manual(A1)
Ur <- svd_res$u[, 1:3]

# Matrica transformacije
C <- solve(t(B) %*% B) %*% t(B) %*% Ur

#sa svim nulama, ovo je dokaz da imamo valjanu matricu baze slike
#iako nije standardna baza sa "lepim brojevima"
round(B %*% C - Ur, 8)


# 2. ZADATAK — Ortogonalna projekcija korišćenjem SVD


v1 <- c(0,-1,2,0,2)
v2 <- c(1,-3,1,-1,2)
v3 <- c(-3,4,1,2,1)
v4 <- c(-1,-3,5,0,7)

x <- matrix(c(-9,-1,-1,4,1), ncol = 1)

# Matrica generisana vektorima potprostora
Umat <- cbind(v1, v2, v3, v4)

# SVD dekompozicija matrice Umat
svd2 <- svd_manual(Umat)

U <- svd2$u
d <- svd2$d

tol <- 1e-8
r <- sum(d > tol)

cat("Dimenzija potprostora U:", r, "\n\n")

# Ortonormirana baza potprostora = prve r kolona U
U_basis <- U[, 1:r, drop = FALSE]

cat("Ortonormirana baza potprostora:\n")
print(U_basis)

# ---------------------------------------------------------
# (a) Ortogonalna projekcija
# ---------------------------------------------------------
x_proj <- U_basis %*% t(U_basis) %*% x

cat("\nProjekcija π_U(x):\n")
print(x_proj)

# ---------------------------------------------------------
# (b) Udaljenost
# ---------------------------------------------------------
dist <- norm(x - x_proj, type = "2")

cat("\nUdaljenost d(x,U):", dist, "\n")


# 3. ZADATAK — Dijagonalizacija uz SVD inverz


A3 <- matrix(c(
  4,-3,-2,
  2,-1,-2,
  3,-3,-1
), nrow = 3, byrow = TRUE)

x1 <- matrix(c(1,0,1), ncol = 1)
x2 <- matrix(c(1,1,0), ncol = 1)
x3 <- matrix(c(1,1,1), ncol = 1)

lambda1 <- 2
lambda2 <- 1
lambda3 <- -1

# ----- Provera karakterističnih vektora -----
cat("Provera karakterističnih vektora:\n")
print(A3 %*% x1 - lambda1 * x1)
print(A3 %*% x2 - lambda2 * x2)
print(A3 %*% x3 - lambda3 * x3)

# ----- Rekonstrukcija A preko SVD -----
svdA3 <- svd_manual(A3)
cat("\nRekonstrukcija A3 iz SVD:\n")
print(round(svdA3$u %*% svdA3$D %*% t(svdA3$v), 6))

# ----- Matrice za dijagonalizaciju -----
P <- cbind(x1, x2, x3)
D <- diag(c(lambda1, lambda2, lambda3))

# ----- Inverz matrice P preko SVD -----
svdP <- svd_manual(P)

SigmaP_inv <- matrix(0, 3, 3)
for(i in 1:3) {
  SigmaP_inv[i,i] <- 1 / svdP$d[i]
}

P_inv <- svdP$v %*% SigmaP_inv %*% t(svdP$u)

cat("\nPseudoinverz (inverz) matrice P:\n")
print(P_inv)

# ----- Provera dijagonalizacije -----
A3_check <- P %*% D %*% P_inv

cat("\nProvera A = P D P^{-1}:\n")
print(round(A3_check, 6))

cat("\nRazlika (A3_check - A3):\n")
print(round(A3_check - A3, 8))
