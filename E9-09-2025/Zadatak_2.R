library(matlib)


#Zadatak 1
#a

A = matrix(c(
  2,  3,  0,  1,
  -1,  1,  2,  0,
  0,  2,  1,  1,
  1, -1, -4,  2,
  3,  0, -2, -1
), ncol = 4, byrow = TRUE)
print(A)

#SVD
svd_A =svd(A)

U = svd_A$u
print(U)

Sigma = svd_A$d

V = svd_A$v
print(V)

tol = 1e-10
r = sum(Sigma > tol)



dim_Im = r
basis_Im = U[, 1:r]



#ker
#ker
n = ncol(A)                               
dim_Ker = n - r                           

if(dim_Ker > 0){
  basis_Ker = V[, (r+1):n, drop = FALSE]  
} else {
  basis_Ker = NULL                        
}





#b
#da li je y u Im i projekcija

y = matrix(c(-8, 3, -5, -17, 1), ncol = 1)
print(y)

U_r = U[, 1:r]


#Projekcija na Im(A)
P = U_r %*% t(U_r)       
y_proj = P %*% y

round(y - y_proj, 10) 


#y jeste u slici




#c

Sigma_inv = diag(1 / Sigma[1:r])

A_pinv = V[,1:r] %*% Sigma_inv %*% t(U[,1:r])

x1 = A_pinv %*% y         
x1



#Provera
round(A %*% x1 - y, 10)







#d

if(dim_Ker == 0){
  x1                                       
} else {
  k = basis_Ker[, 1, drop = FALSE]         
  x2 = x1 + k                              
  x2
  
  
  #Provera
  round(A %*% x2 - y, 10)                  
  
  list(posebno_resenje = x1, drugo_resenje = x2, baza_kernela = basis_Ker)
}







#zadatak 2
#Projekcija i udaljenost

u1 = c(0, -1, 2, 0, 2)
print(u1)
u2 = c(1, -3, 1, -1, 2)
print(u2)
u3 = c(-3, 4, 1, 2, 1)
print(u3)
u4 = c(-1, -3, 5, 0, 7)
print(u4)

x = matrix(c(-9, -1, -1, 4, 1), ncol = 1)
print(x)

A = cbind(u1, u2, u3, u4)
print(A)

svd_A = svd(A)

U = svd_A$u
Sigma = svd_A$d

tol = 1e-10
r = sum(Sigma > tol)

U_r = U[,1:r]

x_proj = U_r %*% t(U_r) %*% x
x_proj

dist = sqrt( sum( (x - x_proj)^2 ) )
dist










#zadatak 3

A3 = matrix(c(
  4, -3, -2,
  2, -1, -2,
  3, -3, -1
), ncol = 3, byrow = TRUE)
print(A3)


v1 = matrix(c(1,0,1), ncol = 1)
print(v1)
v2 = matrix(c(1,1,0), ncol = 1)
print(v2)
v3 = matrix(c(1,1,1), ncol = 1)
print(v3)



A3 %*% v1     
A3 %*% v2     
A3 %*% v3     



#lambda
lambda1 = 2
lambda2 = 1
lambda3 = -1




P = cbind(v1, v2, v3)
D = diag(c(lambda1, lambda2, lambda3))


P_inv = solve(P)




#Provera dijagonalizacije
round(A3 - P %*% D %*% P_inv, 10)   

