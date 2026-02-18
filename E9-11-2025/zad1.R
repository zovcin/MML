library(matlib)   # zbog GramSchmidt()

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

#provera nad matricom A
A = matrix(c(
  1,  3,  1, 2, 1,  6,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 5, 1,  7,
  1,  3,  1, 2, 1,  6,
  2,  6,  2, 4, 2, 12
), ncol = 6, byrow = TRUE)

#objekat koji cuva povratne vrednsoti manualno pisanog svd-a
rezultat_manual = svd_manual(A)

rezultat_classic = svd(A)

zapsmall(rezultat_manual$d)
zapsmall(rezultat_classic$d)

#rekonstrukcija matrice A koriscenjem povratnih vrednosti manualno pisanog svd-a
A_new = rezultat_manual$u %*% rezultat_manual$D %*% t(rezultat_manual$v)
print(A_new)
print(A)
