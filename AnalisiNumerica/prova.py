# -*- coding: utf-8 -*-
"""
Created on Thu May 27 11:00:31 2021

@author: feder
"""


import numpy as np 




###### METODO NEWTON A PIU VARIABILI

print(" METODO NEWTON A PIU VARIABILI")


def jacobian(x):

    x1=float(x[0])
    x2=float(x[1])
    x3=float(x[2])
    
    row_0=5*x1**4,3*x2**2*x3**4,4*x2**3*x3**3
    row_1=2*x1*x2*x3,x1**2*x3,x1**2*x2
    row_2=0,0,4*x3**3
    
    
    J=np.vstack([row_0,row_1,row_2])
    
  
    
    return J

def F(x0,f1,f2,f3):
    
    x1=x0[0]
    x2=x0[1]
    x3=x0[2]
    
    return np.array([f1(x1,x2,x3),f2(x1,x2,x3),f3(x1,x2,x3)])


### FUNZIONI
f1=lambda x1,x2,x3:x1**5+x2**3*x3**4+1
f2=lambda x1,x2,x3:x1**2*x2*x3
f3=lambda x1,x2,x3:x3**4-1

## SOLUZIONE
sol_ZERO=np.array([-1,0,-1]).reshape(-1,1) 
##GUESS INIZIALE
x_0=np.array([-100.0,0,-100]) 
# print(x_0)
# print(s_2)


z_k=np.zeros_like(x_0)-x_0
print("Z_k:\n",z_k)



def Jacobi(A,b,guess_iniziale=np.array([[100.0],[0.0],[0.0]]),iterazioni_max=50,verbose=False):

    #print ("JACOBI")
    x=guess_iniziale ##guess iniziale
    iterazioni=iterazioni_max
    
    #print("Determinante di A",np.linalg.det(A),"\n")
    
    for i in range(iterazioni):
        x_=x.copy() # x_k 
        #print("ITER",i)
        A=jacobian(x_)
        b= -F(x_,f1,f2,f3)
        for j in range(A.shape[0]):
            tmp=0
            for k in range(A.shape[0]):
                if k!=j:
                    tmp+=A[j][k]*x_[k]
            x[j]=(b[j]-tmp)/A[j][j]+x_[j]
            #print("B:\n",b)
            #print("tmp;\n",tmp)
        # print("X_k+1",x)
        # print("x_k",x_)
        #     #assert False
            #x[j]=-x[j]
        #print(x)
        ##Criterio arresto
        if np.allclose(np.dot(A,x)-b,np.zeros_like(b)):
            break
            
    
    if verbose:
        print("Convergenza in {} iterazioni".format(i))
        print("A:\n",A,"\nb:\n",b)          
        print("x:\n",x)
        print("Errore:\n",np.dot(A,x)-b)
        
    return x.reshape(-1)

# def Jacobi(J,Fm,guess_iniziale=np.array([[100.0],[0.0],[0.0]]),iterazioni_max=150,verbose=False):

#     print ("JACOBI")
#     x=guess_iniziale ##guess iniziale
#     iterazioni=iterazioni_max
    
#     print("Determinante di J",np.linalg.det(J))
#     x=np.array([-100,0,-100]).reshape(-1)
    
#     for i in range(iterazioni):
#         x_=x.copy() # x_k 
#         J=jacobian(x_)
#         Fm= -F(x_,f1,f2,f3)
#         for j in range(J.shape[0]):
#             tmp=0
#             for k in range(J.shape[0]):
#                 if k!=j:
#                     tmp+=J[j][k]*x_[k]
#             x[j]=(Fm[j]-tmp)/J[j][j]+x_[j]
        
        
#         #x=x+x_
        
#         ##Criterio arresto
#         if np.allclose(np.dot(J,x)-Fm,np.zeros_like(Fm)):
#             break
            
    
#     if verbose:
#         print("Convergenza in {} iterazioni".format(i))
#         print("A:\n",J,"\nb:\n",F)          
#         print("x:\n",x)
#         print("Errore:\n",np.dot(J,x)-F)
        
#     return x.reshape(-1)
#     #print(np.linalg.det(A))


def prova_method(x_k):
    
    x_old=x_k.copy()
    
    lst=[]
    print("X_k:\n",x_k)
    
    # print("quii")
    # for i in range(25):
    #     J=jacobian(x_k)
    #     F_=F(x_k,f1,f2,f3)
    #     x_new=x_k-np.dot(np.linalg.inv(J),F_)
        
    #     x_k=x_new.copy()
    #     print(x_k)
    #     lst.append(x_k)
    
    
    print("SASDSADADAAD")
    
    #JACOBI SOLITO
    for i in range(100):
        x_=x_old.copy() # x_k 
        #J=jacobian(x_old)
        Fm= -F(x_old,f1,f2,f3)
        for j in range(J.shape[0]):
            tmp=0
            for k in range(J.shape[0]):
                if k!=j:
                    tmp+=J[j][k]*x_[k]
            x_old[j]=(Fm[j]-tmp)/J[j][j]-x_[j]
        
            print(x_old)
        
        #x=x+x_
        
        ##Criterio arresto
        if np.allclose(np.dot(J,x_old)-Fm,np.zeros_like(Fm)):
            break
    
    
    
    print(x_old)
    return x_old,lst


def jacobii(A,b,N=50,x=None):
    """Solves the equation Ax=b via the Jacobi iterative method."""
    # Create an initial guess if needed                                                                                                                                                            
    if x is None:
        x = np.zeros(len(A[0]))
        
    lst=[]
    x=np.array([1,1000,0.1]).reshape(-1)
    # Create a vector of the diagonal elements of A                                                                                                                                                
    # and subtract them from A                                                                                                                                                                     
    D = np.diag(A)
    R = A - np.diagflat(D)

    # Iterate for N times                                                                                                                                                                          
    for i in range(50):
        x_=x.copy()
        
        A=jacobian(x)
        b=-F(x,f1,f2,f3)
        D = np.diag(A)
        R = A - np.diagflat(D)
        # print(R)
        # print(A)
        # print(D)
        # print(np.diagflat(D))
        # assert False
        #print(b)
        #print(np.dot(R,x_))
        # assert False
        x =(b - np.dot(R,x_)) / D +x_
        #x=-x
        #print(x)
        lst.append(x.copy())
        # if np.allclose(np.dot(A,x)-b,np.zeros_like(b)):
        #     print(np.dot(A,x)-b)
        #     break
    
        
    
    return x,lst



x_0=np.array([-1000.0,-1000.0,-1000]) 
x,lise=jacobii(jacobian(x_0),-F(x_0,f1,f2,f3),x=x_0)
print("ererere",x)

# J=jacobian(x_0)

# print("F_x:\n",F_x)
# print("J:\n",J)
# print(x_0)
q=100

x_0=np.array([-100.0,0,-100]) 
F_x=-F(x_0,f1,f2,f3)
x=Jacobi(jacobian(x_0),F_x,guess_iniziale=x_0,iterazioni_max=q,verbose=False)
print("X:\n",x)

x_0=np.array([-1000.0,-1000.0,-1000]) 
F_x=-F(x_0,f1,f2,f3)
x=Jacobi(jacobian(x_0),F_x,guess_iniziale=x_0,iterazioni_max=q,verbose=False)
print("X:\n",x)

x_0=np.array([0.1,0.1,0.1]) 
F_x=-F(x_0,f1,f2,f3)
x=Jacobi(jacobian(x_0),F_x,guess_iniziale=x_0,iterazioni_max=q,verbose=False)
print("X:\n",x)

x_0=np.array([-100.0,0,100]) 
F_x=-F(x_0,f1,f2,f3)
x=Jacobi(jacobian(x_0),F_x,guess_iniziale=x_0,iterazioni_max=q,verbose=False)
print("X:\n",x)





# print(f1(x[0],x[1],x[2]))
# print(f2(x[0],x[1],x[2]))
# print(f3(x[0],x[1],x[2]))


# print("LISE")
# print(lise[:50])
# print("LST")
# print(lst[:50])



###CONTINUAREEEEEEEEEEEEEEEEEEE


