library(matlib)

# =========================================================
# Spall’s Polynomial — Gradient Descent (n = 4) zadatak 5.1
# =========================================================

n <- 4

# Gornja trougaona matrica jedinica
A <- matrix(0, n, n)
for (i in 1:n) {
  for (j in i:n) {
    A[i,j] <- 1
  }
}

# Početna tačka x0 = 0.2 (1,1,1,1)
x <- matrix(rep(0.2, n), ncol = 1)

alpha <- 0.001     # mali korak (bitno zbog polinoma)
n_iter <- 100

cat("Matrica A:\n")
print(A)

cat("\nPocetna tacka:\n")
print(x); cat("\n")

gradient_f <- function(x) {
  
  y <- A %*% x        # y = Ax
  
  # Kvadratni deo
  grad_quad <- 2 * t(A) %*% A %*% x
  
  # Kubni deo
  grad_cubic <- 0.1 * 3 * t(A) %*% (y^2)
  
  # Četvrti stepen
  grad_quartic <- 0.01 * 4 * t(A) %*% (y^3)
  
  return(grad_quad + grad_cubic + grad_quartic)
}

for (k in 1:n_iter) {
  
  grad <- gradient_f(x)
  x <- x - alpha * grad
  
  cat("Iteracija", k, ":\n")
  print(x); cat("\n")
}

# Konacni komentar resenja koje smo dobili je da x konvergira ka minimumu koji je x=0
# Verovatno je nedovoljan broj iteracija da se dodje do malo preciznijeg resenja

cat("=====================================\n")
cat("Konacna aproksimacija:\n")
print(x)
