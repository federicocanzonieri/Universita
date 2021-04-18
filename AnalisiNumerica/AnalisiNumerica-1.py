# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np


##METODO DELLE SOSTITUZIONI IN AVANTI (FORWARD)
def lower_system_resolve(lower,b,verbose=False):
    
    solution=np.zeros((lower.shape[0]))
    solution[0]=(b[0]/lower[0][0])
    
    for i in range(1,lower.shape[0]):
    
        tmp=np.dot(solution[0:i],lower[i][0:i])    
        solution[i]=(b[i]-tmp)/lower[i][i]
    
    if verbose:
        print("Lower:\n",lower)
        print("Termini noti:\n",b)
        print("Solution:\n",solution)
        print("Verifying:\n",np.dot(lower,solution))
    
    return solution


##METODO SOSTITUZIONI IN INDIETRO(BACKWARD)
def upper_system_resolve(upper,b,verbose=False):
    
    solution=np.zeros((upper.shape[0]))
    solution[-1]=(b[-1]/upper[-1][-1])
    
    for i in range(upper.shape[0]-2,-1,-1):
        tmp=np.dot(solution[i:],upper[i][i:])    
        solution[i]=(b[i]-tmp)/upper[i][i]
    
    if verbose:
        print("Upper:\n",upper)
        print("Termini noti:\n",b)
        print(solution)
        print(np.dot(upper,solution))
    
    return solution







#####TESTING LOWER RESOLUTION

lower=np.array([[4,0,0],[5,2,0],[8,2,3]])
b=np.array([[1],[2],[3]])
solution=np.array([[0.25],[3/8],[1/12]])
solution_=lower_system_resolve(lower, b,verbose=True)

assert np.allclose(solution,solution)
print("Test lower passed\n")



#####TESTING UPPER RESOLUTION

upper=np.array([[8,2,3],[0,5,2],[0,0,4]])
b=np.array([[1],[2],[3]])
solution=np.array([[0.25],[3/8],[1/12]])
solution_=upper_system_resolve(upper, b,verbose=True)

assert np.allclose(solution,solution)
print("Test upper passed\n")




#####CHOLESKY
#####RETURN THE L(matrix), A=L*L^T
def cholesky(matrix,verbose=False):
    
    lower=np.zeros((matrix.shape))
    lower[0][0]=np.sqrt(matrix[0][0])
    
    for i in range(1,lower.shape[0]):
        for j in range(0,i+1):
            tmp=0
            if i!=j:
                for k in range(0,j):
                    
                        tmp+=lower[i][k]*lower[j][k]
                        
                lower[i][j]=(matrix[i][j]-tmp)/lower[j][j]
            else:
                for k in range(0,i):
                    tmp+=lower[i][k]**2
                lower[i][j]=np.sqrt(matrix[i][j]-tmp)
                  
    if verbose:
        print("Lower:\n",lower)
        print("A:\n",matrix)
        print("Verifying:\n",np.dot(lower,lower.T))
        
    return lower
    
#####TESTING CHOLESKY DECOMPOSITION
A=np.array([[6,15,55],[15,55,225],[55,225,979]])
lower=np.array([[2.44948974, 0.,0],[6.12372436 ,4.18330013,0],[22.45365598,20.91650066,6.11010093]])
lower_=cholesky(A,True)

assert np.allclose(lower,lower_)
print("Test Cholesky passed\n")




#####L=
#### RIFORMULAZIONE GAUSS A=L~U
#### L~=L1^-1*L2^1*.....LN_1^-1
def gauss_decomposition(A,b,verbose=False,solve=False):
    
    ##LIST MATRIX L_i
    L=[]
    L_inverse=[]
    A_=A.copy()
    
    for i in range(A.shape[0]-1):
        
        #CRITERIO DI CROUT (diagonali pari a 1)
        L_=np.diag([1 for _ in range(A.shape[0])])
        L_inv=np.diag([1 for _ in range(A.shape[0])])
        
        for j in range(i+1,A.shape[0]):
            L_[j][i]=-A_[j,i]/A_[i][i]
            L_inv[j][i]=A_[j,i]/A_[i][i]
            
            
        L.append(L_)
        L_inverse.append(L_inv)
        A_=np.dot(L_,A_)    
        
        if verbose:
            print("L{}".format(i+1))
            print(L_)
            print("A{}".format(i+2))
            print(A_)
        
   
    
    L_tilde=L_inverse[0]
    for L_i in L_inverse[1:]:
        L_tilde=np.dot(L_tilde,L_i)
    
    if verbose:
        print("L_tilde")
        print(L_tilde)
        print("U matrix")
        print(A_)
    
    if solve:
        #RISOLVO IL NUOVO SISTEMA
        y=lower_system_resolve(L_tilde, b)
        x=upper_system_resolve(A_,y)
        
        if verbose:
            print("X:\n",x)
        #restituisco anche la sol
        return L_tilde,A_,x
    #solo la decomposition
    return L_tilde,A_

#####TESTING GAUSS DECOMPOSITION/RESOLUTION
A=np.array([[1,2,3],[4,5,9],[7,8,9]])
b=np.array([[1],[2],[3]])
solution=np.array([[-0.33333333,0.66666667,-0]])
L,U,x=gauss_decomposition(A, b,solve=True,verbose=True)

assert np.allclose(np.dot(L,U),A)
assert np.allclose(solution,x)

print("Test Decomposition Gauss passed")







