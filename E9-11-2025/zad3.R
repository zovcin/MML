library(matlib)

# =========================================================
# Example 7.1 — MML-Book (gradijentni postupak)
# =========================================================

# Definicija matrice iz funkcije
A <- matrix(c(
  2, 1,
  1, 20
), nrow = 2, byrow = TRUE)

# Definicija vektora b
b <- matrix(c(5, 3), ncol = 1)

# Početna tačka x0
x <- matrix(c(-3, -1), ncol = 1)

# Korak učenja (mora biti dovoljno mali)
alpha <- 0.05

# Broj iteracija
n_iter <- 100

cat("Pocetna tacka x0:\n")
print(x)
cat("\n")

# ---------------------------------------------------------
# Gradijentni postupak
# ---------------------------------------------------------
for (k in 1:n_iter) {
  
  # Gradijent ∇f(x) = A x − b
  grad <- A %*% x - b
  
  # Iterativni korak
  x <- x - alpha * grad
  
  cat("Iteracija", k, ":\n")
  print(x)
  cat("\n")
}

cat("=====================================\n")
cat("Konacna aproksimacija minimuma:\n")
print(x)

# ---------------------------------------------------------
# Tačno rešenje (analitički minimum)
# ---------------------------------------------------------
x_star <- solve(A) %*% b

cat("\nTacno resenje (A^{-1} b):\n")
print(x_star)

cat("\nRazlika (numericko - tacno):\n")
print(x - x_star)
