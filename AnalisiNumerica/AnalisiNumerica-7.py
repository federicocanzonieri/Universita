# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 16:28:21 2021

@author: feder
"""

####AUTOVALORI 


import numpy as np



def L2(x):
    return np.sqrt(np.sum(x**2))


def rayleigh(x_k,x_k_1):
    return np.dot(x_k.T,x_k_1)/L2(x_k)**2


##METODO POTENZE
def potenze_autovalori(A,MAX_ITER=10,norm=True):
    
    
    x_=np.array([3,10,9])
    x_=x_/L2(x_)
    for i in range(MAX_ITER):
        x=np.dot(A,x_)
        #print(x)
        eigvals=(rayleigh(x_,x))
        
        ##NORM SENNO NON CONVERGE
        if norm:
            x=x/L2(x)
        x_=x.copy()
    
    #AUTOVETTORI ASSOCIATI ALL'AUTOVALORE MASSIMO
    return -x,eigvals
        
print("TESTING METODO POTENZE")
A=np.array([[1,7,8],[2,9,4],[3,1,9]])

print("Matrice:\n",A,"\nShape:",A.shape)
print("Autovettore:",np.linalg.eig(A))
l=np.linalg.eig(A)[0][1]


eigvec,eigval=potenze_autovalori(A)
print("Autovettore Stimato:\n",eigvec)
print("Autovettore Stimato:\n",eigval)


assert np.allclose(np.abs(np.linalg.eig(A)[1][:,1]),np.abs(eigvec))
assert np.allclose(np.abs(l),np.abs(eigval))
print("\nTEST METODO POTENZE PASSED\n")




def potenze_autovalori_inverse(A,MAX_ITER=10,norm=True):

    return potenze_autovalori(np.linalg.inv(A),MAX_ITER,norm)


##TEST METODO POTENZE INV

print("TESTING METODO POTENZE INVERSE")
#A=np.array([[-40,100,18],[13,55,4],[13,51,9]])
A_inv=np.linalg.inv(A)
print("Matrice:\n",A_inv,"\nShape:",A_inv.shape)
# print("Autovettore:",np.linalg.eig(A))
l=np.linalg.eig(A_inv)[0][0]
print(np.linalg.eig(A_inv))


eigvec,eigval=potenze_autovalori_inverse(A)
print("Autovettore Stimato:\n",eigvec)
print("Autovettore Stimato:\n",eigval)
print("Autovettore Stimato(A orig):\n",1/eigval)


assert np.allclose(np.abs(np.linalg.eig(A_inv)[1][:,0]),np.abs(eigvec))
assert np.allclose(np.abs(l),np.abs(eigval))
print("\nTEST METODO POTENZE INVERSE PASSED\n")


def J_matrix(i,k,n,x):
    
    
    r=(np.sum(x[i:,i]**2)**0.5)
    # print(x[i:,i])
    # print(r)
    s=x[k,i]/r
    c=x[i,i]/r
    
    I=np.eye(n,dtype=float)
    I[i][i]=c
    I[i][k]=s
    I[k][i]=-s
    I[k][k]=c
    
    # print(I)
    
    return I

def givens(i,k,n,A):
    
    J=J_matrix(i,k,n,A)
    
    A_=np.dot(J,A)
    return A_,J


print("TESTING METODO GIVENS")
#x=np.array([1,4,2,9]).reshape(-1,1)
# A=np.array([[6,5,0],[5,1,4],[0,4,3]])
# A2,G1=givens(0,1,3,A)
# A3,G2=givens(1,2,3,A2)
# Q=np.dot(G1.T,G2.T)


# for _ in range(100):
#     A3=np.dot(A3,Q)

#     A2,G1=givens(0,1,3,A3)
#     A3,G2=givens(1,2,3,A2)
    
#     Q=np.dot(G1.T,G2.T)
    

# print(np.diag(np.dot(A3,Q)))
# print(np.linalg.eigvals(A))

def norm(x):
    return np.sqrt(np.sum(x**2))

def prod_scalare(x,y):
    return np.abs(np.sum(np.multiply(x,y)))

def Householder(A):
    
    print("A:\n",A)
    A_=A.copy()
    Q_matrix=[]
    for idx in range(A.shape[0]-1):
       
        x=A_[:,idx].copy().reshape(-1,1)
        x[:idx]=0
        print("x:\n",x)    
        e=np.zeros(A.shape[0]).reshape(-1,1)
        e[idx]=1
        print("e:\n",e)
        u=x-norm(x)*e
        v=u/norm(u)
        print("u:\n",u)
        print("v:\n",v)
        Q_=np.eye(A.shape[0])-2*v*v.T
        R_=np.dot(Q_,A_)
        A_=R_
        print("Q:\n",Q_)
        print("R:\n",R_)
        #tt=np.dot(Q_,R_)
        print("QR:\n",np.dot(Q_,R_))
        #print(A_)
        Q_matrix.append(Q_)
    
    
    #print("ere",np.dot(Q_matrix[0],tt))
    
    return R_,Q_matrix

def sign(x):
    if x>0:
        return 1
    else:
        return -1

def Householder_mod(A):
    print("MOd")
    print("A:\n",A)
    A_=A.copy()
    Q_matrix=[]
    for idx in range(A.shape[0]-2):
        print(idx)
        print(A_)
        x=A_[idx+1:,idx].copy().reshape(-1,1)
        # x[:idx]=0
        print("x:\n",x)    
        e=np.zeros(A.shape[0]-idx-1).reshape(-1,1)
        e[0]=1
        print("e:\n",e)
        u=x+norm(x)*e*sign(x[0])
        v=u/norm(u)
        print("u:\n",u)
        print("v:\n",v)
        Q_=np.eye(A.shape[0])
        Q_[0][0]=1
        print(Q_)
        
        Q_[idx+1:,idx+1:]=np.eye(A.shape[0]-idx-1)-2*v*v.T
        print(Q_)
        R_=np.dot(Q_,A_)
        #A_=np.dot(R_,Q_)
        print("Q:\n",Q_)
        print("R:\n",R_)
        #tt=np.dot(Q_,R_)
        print("QR:\n",np.dot(Q_,R_))
        print("RQ:\n",np.dot(R_,Q_))
        A_=np.dot(R_,Q_.T) ##simile
        #assert False
        #print(A_)
        Q_matrix.append(Q_)
    
    
    #print("ere",np.dot(Q_matrix[0],tt))
    
    return R_,Q_matrix
       


with np.printoptions(2):
    print("QUIIIIIIIIIIIIIIIII")
    A=np.array([[1,2,3,4],[4,4,4,4],[0,1,-1,1],[0,0,2,3]])
    A=np.array([[1,2,2,1],[2,3,7,1],[2,7,5,2],[1,1,2,6]])
    A=np.array([[1,2,1,2],[3,1,3,2],[6,1,2,2],[6,1,3,4]])
    #A=np.array([[5,-2,2],[4,-3,4],[3,-6,7]])
    
    R,Q=Householder_mod(A)
    print("R:\n",R)
    print("Q[0]:\n",Q[0])
    print("Q[1]:\n",Q[1])
    #print("Q[2]:\n",Q[2])
    hh=np.dot(np.dot(Q[0],A),Q[0].T)
    print("HESSIN\n",np.dot(np.dot(Q[1],hh),Q[1].T))
    Q__=(np.dot(Q[1],Q[0]))
    print("RQ:\n",np.dot(R,Q__))
    pre=np.dot(R,Q__)
    print(pre)
    print("QHQ:\n",np.dot(np.dot(Q__,R),Q__.T))
    assert False
    

print("\n\n\n\\n\n")
A=np.array([[1,4,3],[4,1,2],[3,2,1]])
#A=np.array([[1,2,1],[1,7,3],[1,1,10]])
print("EIG NUMPY",np.linalg.eigvals(A))

A1=A.copy()
for _ in range(1):
    
    R,Q=Householder(A1)
    #print("\n\n\n\n\n")
    #print(A1)

    q0=Q[0]
    q1=Q[1]
    #q2=Q[2]
    pr=(np.dot(q0.T,q1.T))
    print(np.linalg.eigvals(R))
    
    #print(np.linalg.det(pr.T))
    
    print(np.linalg.eigvals(A))
    
    print(np.diag(R))
    A1=(np.dot(R,pr))
    #print(A1)
    
print(A1)
print(np.diag(A1))
print("EIG NUMPY",np.linalg.eigvals(A))
print("EIG NUMPY",np.linalg.det(A))

assert False


##CONTINUARE



print("GRAM SCHMIDT")
#A=np.array([[-1,-1,1],[1,3,3],[-1,-1,5],[1,3,7]])
def gram_schmidt(A):

    print("A:\n",A)

    R=np.zeros((3,3))
    q=np.zeros(A.shape)
    
    for idx in range(3):
        v=A[:,idx].copy()
        b=A[:,idx].copy()
        for j in range(idx):     
            R[j][idx]=prod_scalare(q[:,j],b )/prod_scalare(q[:,j], q[:,j])
            
            v=v-q[:,j]*R[j][idx]
        
        R[idx][idx]=norm(v)
        q_=v/norm(v)
        #q.append(q_)
        q[:,idx]=q_
        #print(q)
        #assert False
    print(q)
    print(R)
    
    return q,R
    

A=np.array([[1,0,0],[0,1,0],[0,0,1]])
qq,rr=gram_schmidt(A)


v1=qq[:,0]
v2=qq[:,1]
v3=qq[:,2]

print(v1,v2,v3)
print(prod_scalare(v1, v1))
print(prod_scalare(v1, v3))
print(prod_scalare(v3, v2))


# r,q=Householder(A)

# print("we",np.dot(q,r))
# print("ereeeee",np.dot(qq,rr))

# q0=q[0]
# q1=q[1]
# print(np.dot(np.dot(q0,q1),r))

# print(np.dot(q0,q1))
# print(r)


