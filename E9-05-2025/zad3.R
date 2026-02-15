A = matrix(c(2, 1, 1, 20), nrow = 2)
b = matrix(c(5, 3), ncol = 1)
x = matrix(c(-3, -1), ncol = 1)
gamma = 0.085

f = function(x) {
  0.5 * t(x) %*% A %*% x - t(b) %*% x
}

hist = data.frame(i = integer(), x1 = numeric(), x2 = numeric(), fx = numeric())
for (i in 1:100) {
  grad = t(t(x) %*% A - t(b))
  x = x - gamma * grad
  if (i %% 10 == 0) {
    hist = rbind(hist, data.frame(i = i, x1 = x[1], x2 = x[2], fx = as.numeric(f(x))))
  }
}
print(hist, digits = 13)

hist[nrow(hist), "x1"]
hist[nrow(hist), "x2"]
# x = [2.49, 0.03]
hist[nrow(hist), "fx"] 
# f(x) = -6.26