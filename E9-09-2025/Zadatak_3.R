library(matlib)                                   





Q = matrix(c(2, 1,                                 
             1, 20),                               
           ncol = 2, byrow = TRUE)  

print(Q)
b = matrix(c(5, 3), ncol = 1) 
print(b)




#F-ja f                                   

f = function(x){                                   
  return( 0.5 * t(x) %*% Q %*% x - t(b) %*% x )    
}                                                  




#Graijent                                       

grad_f = function(x){                               
  return( Q %*% x - b )                             
}                                                 
print(grad_f)



# Parametri                                       

gamma = 0.085  


n_iter = 100                                       

                                   
x = matrix(c(-3, -1), ncol = 1)  
print(x)




# Cuvanje iteracija                              
X_iteracije = matrix(0, nrow = 2, ncol = n_iter + 1) 
X_iteracije[,1] = x                                



for(k in 1:n_iter){                                
  x = x - gamma * grad_f(x)                        
  X_iteracije[,k+1] = x                            
}                                                  




#Rezultati

x                                                  
f(x)                                               
