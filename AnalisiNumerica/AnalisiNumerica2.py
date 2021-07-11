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
# A=np.array([[1,2,3],[4,5,6],[7,8,9]])
# L,U=Doolittle(A,verbose=True)

# assert np.allclose(np.dot(L,U),A)
# print("Doolittle test passed\n")



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
# A=np.array([[1,2,3],[4,5,6],[7,8,9]])
# L,U=Croout(A,verbose=True)

# assert np.allclose(np.dot(L,U),A)
# print("Croout test passed\n")


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


   
# ###TESTING THOMAS 
# #Trig matrix
# A[2,0]=0
# A[0,2]=0
# alpha,b,gamma=Thomas(A,diagonal_only=True)
# l1,l2=Thomas(A,verbose=True)

# assert np.allclose(np.dot(l1,l2),A)
# print("Thomas test passed\n")



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


# ###TESTING JACOBI
# A=np.array([[5,-2,3],[-3,9,1],[3,-1,-7]],dtype=np.float64)
# b=np.array([[-1],[4],[3]],np.float64)
# solution=np.array([93/346 , 100/173 , -137/346])

# x_=Jacobi(A,b,verbose=True)


# assert np.allclose(x_,solution)
# print("Jacobi test passed\n")


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


# ###TESTING GAUSS SEIDOL
# A=np.array([[5,-2,3],[-3,9,1],[3,-1,-7]],dtype=np.float64)
# b=np.array([[-1],[4],[3]],np.float64)
# solution=np.array([93/346 , 100/173 , -137/346])

# x_=GaussSeidol(A,b,verbose=True)


assert np.allclose(x_,solution)
print("Gauss Seidol test passed\n")


### MATRICE DI JACOBI

def test_diag_dominante(A,strett=True):
    
    for i in range(A.shape[0]):
        a_ii=A[i][i]
        # print(np.abs(a_ii))
        # print(np.abs(A[i]))
        # print(np.sum(np.abs(A[i]))-np.abs(a_ii))
        if strett:
            if np.abs(a_ii)<=np.sum(np.abs(A[i]))-np.abs(a_ii):
               return False
        
        
    return True


def Jacobi_matrix(A,b,iterazioni_max=50,guess_iniziale=np.array([[0],[0],[0]]),verbose=False):
    
    
    F=-np.triu(A,k=1)
    E=-np.tril(A,k=-1)
    D=A+E+F
    D_inv=np.linalg.inv(D)
    M_J=np.dot(D_inv,(E+F))
    
   
    
    if not test_diag_dominante(A):
        print("Non converge (CS strett diag dominante)")
    
    if (np.max(np.abs(np.linalg.eigvals(M_J))))>1:
        print("Non converge (Raggio spettrale >1)")
    
    x=guess_iniziale
    x_=x.copy()
    
    for it in range(iterazioni_max):
        x=np.dot(M_J,x_)+np.dot(D_inv,b)
        x_=x.copy()

    if verbose:
        print("Matrice di Jacobi:\n",M_J)
        #print(D_inv)
        
    return x
        

# ###TESTING MATRICE DI JACOBI
# print("MATRICE JACOBI")
# A=np.array([[3,-1,1],[2,6,2],[-1,-1,6]],dtype=np.float64)
# b=np.array([[4],[12],[0]],np.float64)
# x=Jacobi_matrix(A,b,verbose=False)
# print("x:\n",x,"\n\n")


def GaussSeidol_matrix(A,b,iterazioni_max=50,guess_iniziale=np.array([[0],[0],[0]]),verbose=False):
    
    F=-np.triu(A,k=1)
    E=-np.tril(A,k=-1)
    D=A+E+F
    D_E_inv=np.linalg.inv(D-E)
    M_GS=np.dot(D_E_inv,(F))
    
    # print(D-E)
    # print(D_E_inv)
    # print(F)
    # print(M_GS)
    
    x=guess_iniziale
    x_=x.copy()
    
    if not test_diag_dominante(A):
        print("Non converge (CS A strett diag dominante)")
    
    if (np.max(np.abs(np.linalg.eigvals(M_GS))))>1:
        print("Non converge (Raggio spettrale >1)")
    else:
        print("Converge (Raggio spettrale <1)")
    ##SIMM DEF POSITIVA DA IMPLEMENTARE
    
    
    for it in range(iterazioni_max):
        x=np.dot(M_GS,x_)+np.dot(D_E_inv,b)
        x_=x.copy()

    if verbose:
        print("Matrice GaussSeidol:\n",M_GS)
        #print(D_E_inv)
        
    return x
    
    
# ###TESTING MATRICE DI GAUSS SEIDOL
# print("MATRICE GAUSS SEIDOL")
# A=np.array([[2,-1,1],[2,2,2],[-1,-1,2]],dtype=np.float64)
# b=np.array([[4],[12],[0]],np.float64)

# x=GaussSeidol_matrix(A,b,verbose=False)
# print("x:\n",x,"\n\n")






###METODO SOR

def SOR_method(A,b,iterazioni_max=50,omega=0.1,guess_iniziale=np.array([[0],[0],[0]]),verbose=False):
    
    print("SOR method:\n Omega:{}".format(omega))
    F=-np.triu(A,k=1)
    E=-np.tril(A,k=-1)
    D=A+E+F
    D_E_inv=np.linalg.inv(D-omega*E)
    D_inv=np.linalg.inv(D)
    L=np.dot(D_inv,E)
    R=np.dot(D_inv,F)
    
    H=np.dot(np.linalg.inv(np.eye(A.shape[0])-omega*L),((1-omega)*np.eye(A.shape[0])+omega*R))
    
    
    x=guess_iniziale
    x_=x.copy()
    
    #convergenza_SOR=(np.max(np.abs(np.linalg.eigvals(H))))
    if omega>0 and omega<2:
        print("CN Convergenza")
    
    
    if test_diag_dominante(A):
        if omega>0 and omega<=1:
            
            print("Converge (CNS A strett diag dominante)")
    

    
    for it in range(iterazioni_max):
        x=np.dot(H,x_)+omega*np.dot(D_E_inv,b)
        x_=x.copy()

    if verbose:
        print("Matrice GaussSeidol:\n",H)
        #print(D_E_inv)
        
    return x
    

# print("METODO SOR")
# A=np.array([[2,-1,1],[2,2,2],[-1,-1,2]],dtype=np.float64)
# b=np.array([[4],[12],[0]],np.float64)

# x=SOR_method(A,b,verbose=False)
# print("x:\n",x)









