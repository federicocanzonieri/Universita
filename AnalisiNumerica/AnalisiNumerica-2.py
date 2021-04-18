import numpy as np



def sum_row_col(i,j,L,U):
    tmp=0
    for k in range(min(i,j)):
        tmp+=L[i][k]*U[k][j]
    return tmp


##DOOLITTLE DECOMPOSITION
def Doolittle(A,verbose=False):
    
    L=np.zeros_like(A)
    U=np.zeros_like(A)
    U[0]=A[0]
    print("Doolittle")
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if i==j:
                L[i][i]=1   ##dolittle 
                U[i][i]=A[i][i]-sum_row_col(i,j,L,U)
            elif i>j:
                L[i][j]=(A[i][j]-sum_row_col(i,j,L,U))/U[j][j]

            elif j>i:
                U[i][j]=(A[i][j]-sum_row_col(i,j,L,U)) ##/L[i][i]
             
    if verbose:
        print("U:\n",U)
        print("L:\n",L)
        print("A:\n",A)
        print("Decomposizione:\n",np.dot(L,U))   
        # print("Inv L",np.linalg.inv(L))
    
    return L,U


##TESTING DOOLITTLE
A=np.array([[1,2,3],[4,5,6],[7,8,9]])
L,U=Doolittle(A,verbose=True)

assert np.allclose(np.dot(L,U),A)
print("Doolittle test passed\n")



## CROUT DECOMPOSITION

def Croout(A,verbose=False):

    L=np.zeros_like(A)
    U=np.zeros_like(A)
    print("Croout")
    L[:,0]=A[:,0]

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if i==j:
                U[i][i]=1 ##crout 
                L[i][i]=A[i][i]-sum_row_col(i,j,L,U)
            elif i>j:
                L[i][j]=(A[i][j]-sum_row_col(i,j,L,U))

            elif j>i:
                U[i][j]=(A[i][j]-sum_row_col(i,j,L,U))/L[i,i]
    
    if verbose:
        print("U:\n",U)
        print("L:\n",L)
        print("A:\n",A)
        print("Decomposizione:\n",np.dot(L,U))   
        
    return L,U


##TESTING CROOUT
A=np.array([[1,2,3],[4,5,6],[7,8,9]])
L,U=Croout(A,verbose=True)

assert np.allclose(np.dot(L,U),A)
print("Croout test passed\n")


###THOMAS DECOMPOSITION
def Thomas(A,diagonal_only=False,verbose=False):
    
    print("THOMAS")
    ##PRENDO GLI ELEMENTI DELLE DIAGONALI
    a=[]
    b=[]
    c=[]
    for i in range(A.shape[0]):
        a.append(A[i][i])
        try:
            b.append(A[i+1][i])
           
        except IndexError:
            pass
        try:
            c.append(A[i][i+1])
        except IndexError:
            pass
            
    a=np.array(a)
    b=np.array(b)
    c=np.array(c)
    alpha=np.zeros_like(A[0])
    gamma=np.zeros(A[0].shape[0]-1)
    
    alpha[0]=a[0]
    gamma[0]=c[0]/alpha[0]
    
    for i in range(1,A.shape[0]):
        
        alpha[i]=a[i]-b[i-1]*gamma[i-1]
        if i<A.shape[0]-1:
            gamma[i]=c[i]/alpha[i]
            
    ##Solo le diagonali
    if diagonal_only:
        if verbose:
            print(alpha,b,gamma)
        return alpha,b,gamma
    
    
    #Costruisco le matrici
    l1=np.zeros_like(A)
    l2=np.zeros_like(A)
    for i in range(A.shape[0]):
        l1[i][i]=alpha[i]
        l2[i][i]=1
        try:
            l2[i][i+1]=gamma[i]
           
        except IndexError:
            pass
            
        try:
            l1[i+1][i]=b[i]
        except IndexError:
            pass
    
    if verbose:
        print("L1:")
        print(l1)
        print("L2:")
        print(l2)
        
    return l1,l2
    
    # print(np.dot(l1,l2))
    # print(A)


   
###TESTING THOMAS 
#Trig matrix
A[2,0]=0
A[0,2]=0
alpha,b,gamma=Thomas(A,diagonal_only=True)
l1,l2=Thomas(A,verbose=True)

assert np.allclose(np.dot(l1,l2),A)
print("Thomas test passed\n")



### ITERATIVI
### JACOBI
def Jacobi(A,b,guess_iniziale=np.array([[100.0],[0.0],[0.0]]),iterazioni_max=50,verbose=False):

    print ("JACOBI")
    x=guess_iniziale ##guess iniziale
    iterazioni=iterazioni_max
    
    print("Determinante di A",np.linalg.det(A))
    
    for i in range(iterazioni):
        x_=x.copy() # x_k 
        for j in range(A.shape[0]):
            tmp=0
            for k in range(A.shape[0]):
                if k!=j:
                    tmp+=A[j][k]*x_[k]
            x[j]=(b[j]-tmp)/A[j][j]
            
        ##Criterio arresto
        if np.allclose(np.dot(A,x)-b,np.zeros_like(b)):
            break
            
    
    if verbose:
        print("Convergenza in {} iterazioni".format(i))
        print("A:\n",A,"\nb:\n",b)          
        print("x:\n",x)
        print("Errore:\n",np.dot(A,x)-b)
        
    return x.reshape(-1)
    #print(np.linalg.det(A))


###TESTING JACOBI
A=np.array([[5,-2,3],[-3,9,1],[3,-1,-7]],dtype=np.float64)
b=np.array([[-1],[4],[3]],np.float64)
solution=np.array([93/346 , 100/173 , -137/346])

x_=Jacobi(A,b,verbose=True)


assert np.allclose(x_,solution)
print("Jacobi test passed\n")


def GaussSeidol(A,b,guess_iniziale=np.array([[100.0],[100.0],[100.0]]),iterazioni_max=50,verbose=False):

    print("GAUSS SEIDOL")
    
    x=guess_iniziale ##guess iniziale
    iterazioni=iterazioni_max
    
    print("Determinante di A",np.linalg.det(A))
        
    for i in range(iterazioni):
        
        for j in range(A.shape[0]):
            tmp=0
            for k in range(A.shape[0]):
                if k!=j:
                    tmp+=A[j][k]*x[k]
            x[j]=(b[j]-tmp)/A[j][j]
            
        ##Criterio arresto
        if np.allclose(np.dot(A,x)-b,np.zeros_like(b)):
            break ##fuori dal ciclo

    if verbose:
        print("Convergenza in {} iterazioni".format(i))
        print("A:\n",A,"\nb:\n",b)          
        print("x:\n",x)
        print("Errore:\n",np.dot(A,x)-b)
                 
    return x.reshape(-1)
    #print(np.linalg.det(A))


###TESTING JACOBI
A=np.array([[5,-2,3],[-3,9,1],[3,-1,-7]],dtype=np.float64)
b=np.array([[-1],[4],[3]],np.float64)
solution=np.array([93/346 , 100/173 , -137/346])

x_=GaussSeidol(A,b,verbose=True)


assert np.allclose(x_,solution)
print("Gauss Seidol test passed\n")


###CONTINAURE CON LE MATRICI DI JACOBI E GAUSS SEIDOL ...COND CONVERGENZA







