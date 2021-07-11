# -*- coding: utf-8 -*-
"""
Created on Tue Apr 27 22:59:53 2021

@author: feder
"""
import numpy as np
import matplotlib.pyplot as plt

##INTERPOLAZIONE POLINOMIALE
##COEFFICENTI INDETERMINATI

def coefficenti_indeterminati(x,y,verbose=False,punti_valutare=None):
    
    if verbose:
        plt.figure()
        plt.scatter(x,y,label="Punti reali",s=30,marker='x')
        plt.grid()
        plt.title("Distribuzione punti")

    ##VANDERMONDE MATRIX
    V=[]
    for i in range(0,x.shape[0]):
        V.append(x**i)
        
    V=np.array(V).T
    ##Risolvo il sistema 
    coef=np.dot(np.linalg.inv(V),y)
    if verbose:
        print("Vandermonde:\n",V)
        print("Coefficenti:\n",coef)
    
    res=[]
    poly_string=""
    if not np.allclose(coef[0],0):
        poly_string=str(float(round(coef[0],4)))+"+"
    
    for index,cd in enumerate(coef[1:]):
        if cd==1:
            poly_string+="x^"+str(index+1)+"+"
        elif np.allclose(cd,0):
            pass
        else:    
            poly_string+=str(round(float(cd),4))+"x^"+str(index+1)+"+"
    
    print("Polinomio: {}".format(poly_string[:-1]))

    if punti_valutare is not None:
        punti=punti_valutare
        for x_ in punti:
            poly=[]
            for grado in range(0,x.shape[0]):
                 poly.append(x_**grado)
            
            res.append(np.dot(coef.reshape(-1),np.array(poly).reshape(x.shape[0],1)))
    
        if verbose:
            plt.scatter(punti,res,label="Punti interpolati")
            plt.legend()
            plt.show()
        
        return np.dot(V,coef),coef,res   
    #print("Interpolata:\n",np.array(res))
    #ritorniamo le y interpolate e i coef e le y nei punti dati in input
    return np.dot(V,coef),coef
    
    
##TESTING COEF INDETERMINATI

x=np.array([0,1,2,3,4,5,6])
y=3*x**2+1

y_stimate,_=coefficenti_indeterminati(x, y)
assert np.allclose(y_stimate,y)
print("Test coefficenti indeterminati passed\n\n")


###INTERPOLAZIONE LAGRANGIANA

def L_i(x_,i,x):
    
    res=1
    x_i=x[i]
    for idx,x_j in enumerate(x):
        if idx!=i:
            res*=(x_-x_j)/(x_i-x_j)
    return res
        
def polinomi_lagrange(x,f,a,verbose=False):
    
    if verbose:
        plt.figure()
        plt.scatter(a,f,marker='x',label='Reali')
        plt.grid()
        
    y_interpolated=[]
    for x_ in x:
        sum_=0
        for i in range(len(a)):
            sum_+=f[i]*L_i(x_,i,a)
            #print(f[i]*L_i(x_,i,a))
            
        y_interpolated.append(sum_)
    if verbose:
        plt.scatter(x,y_interpolated,label='Stimat')
        plt.legend()
        
    return y_interpolated
    

a=np.arange(1,30)
f=a*np.sin(a)+a**2+np.log(a+1)

x_inter=polinomi_lagrange(a, f, a,verbose=False)
x_=polinomi_lagrange(a+0.5, f, a,verbose=False)

assert np.allclose(x_inter,f)
print("Test Lagrangian interpolazione\n\n")


##METODO DIFFERENZE DIVISE RIVEDER

def calculate_nodi_ug_spaziati(l,h,step_size,i):
    
    steps=[]
    passo=(h-l)/i
    for i_ in range(i+1):
        steps.append(l+i_*passo)

    return steps,passo

def differenze_divise(x,l,h,punti):
    
    steps,passo=np.array(calculate_nodi_ug_spaziati(l,h,0.1,punti))
    f=np.cos(steps)
    
    print((steps))
    print(f)
    tmp=f.copy()
    poly=[]
    poly.append(tmp[0])
    for idx in range(1,punti+1):
        lst=np.zeros(tmp.shape[0]-1)
        for i in range(0,lst.shape[0]):
            
            lst[i]=-(tmp[i]-tmp[i+1]/(steps[i]-steps[i+1]))
        
        #print(lst)
        poly.append(lst[0])
        tmp=lst.copy()
        
    print("Polinomio:\n")
    print(poly)
    print("r:",(x-l)/passo)
    r=(x-l)/passo
    print("Valuto polinomio:\n")
    pol=poly[0]
    for idx,coef in enumerate(poly[1:]):
       
        print(coef*r_j(r,idx+1))
        pol+=coef*r_j(r,idx+1)
    
    print(pol)

def r_j(r,j):
    
    prod=1
    for j_ in range(0,j):
        prod*=r-(j_)
    
    return prod



##TESTING DIFFERENZE DIVISE

x=0.44

l=0.3
h=0.6

differenze_divise(x,l,h,3)


#####RIVEDERE DIFF DIVISE




###RUNGE, NODI CHEBYCHEV E INSTABILITA POLINOMIALE

def runge(x):
    return 1/(1+x**2)

def testing_runge():
    ##TESTING RUNGE, CHEBYCHEV
    punti=6
    x=np.linspace(-1,1,punti)
    x_=np.linspace(-1,1,80)
    
    
    plt.figure()
    plt.scatter(x_,runge(x_),label='runge',s=10)
    
    
    ##PROVIAMO CON UN POL INT di grado 5
    y_,c,y__=coefficenti_indeterminati(x,runge(x),False,x_)
    plt.scatter(x,y_,label='punti inter',marker='x',s=25)
    plt.scatter(x_,y__,label='p5',s=10)
    
    
    ##PROVIAMO CON UN POL INT di grado 10
    punti=11
    x=np.linspace(-1,1,punti)
    y_,c,y__=coefficenti_indeterminati(x,runge(x),False,x_)
    #plt.scatter(x,y_,label='p10_',s=10,marker='+')
    plt.scatter(x_,y__,label='p10',s=10)
    
    
    plt.grid()
    plt.title("Runge/p5/p10")
    #plt.legend()
    
    ##PROVA LAGRANGIANA
    
    # x_=np.linspace(-5,5,20)
    # plt.figure()
    # y_lagran=polinomi_lagrange(x_, runge(x), x,verbose=False)
    # plt.scatter(x,runge(x),s=20,marker='x')
    # plt.scatter(x_,y_lagran)
    # plt.title("Runge lagrangiana")


testing_runge()




###POLINOMIO CHEBYCHEV

def T_k(x,k,specie='1'):
    
    if k<0:
        raise ("K cannot be lower than 0")
    
    if k==0:
        return np.ones(x.shape[0])
    elif k==1:
        if specie=='1':
            return x
        else:
            return 2*x
    else:
        return T_k(x,k-1,specie)*2*x-T_k(x,k-2,specie)


def zeri_T_k(k):

    return np.array([ np.cos((2*i+1)*np.pi/(2*k)) for i in range(k)])
    

##TESTING CHEBYCHEV

def testing_chebychev():
    
    x=np.linspace(-5,5,20)
    y=T_k(x,0)
    
    plt.figure()
    plt.scatter(x,T_k(x,0),label='T_0')
    plt.scatter(x,T_k(x,1),label='T_1')
    plt.scatter(x,T_k(x,2),label='T_2')
    plt.scatter(x,T_k(x,3),label='T_3')
    plt.scatter(x,T_k(x,4),label='T_4')
    plt.grid()
    plt.title("Polinomi chebychev")
    plt.legend()
    
    
    y_2=2*x**2-1
    y_3=4*x**3-3*x
    y_4=8*x**4-8*x**2+1
    
    assert np.allclose(T_k(x,2),y_2)
    assert np.allclose(T_k(x,3),y_3)
    assert np.allclose(T_k(x,4),y_4)
    
    
    print("TESTING CHEBYCHEV PASSED")
    
def testing_chebychev_2():
    
    x=np.linspace(-5,5,20)
    y=T_k(x,0)
    
    plt.figure()
    plt.scatter(x,T_k(x,0,'2'),label='T_0')
    plt.scatter(x,T_k(x,1,'2'),label='T_1')
    plt.scatter(x,T_k(x,2,'2'),label='T_2')
    plt.scatter(x,T_k(x,3,'2'),label='T_3')
    plt.scatter(x,T_k(x,4,'2'),label='T_4')
    plt.grid()
    plt.title("Polinomi chebychev")
    plt.legend()
    
    
    y_2=4*x**2-1
    
    assert np.allclose(T_k(x,2,'2'),y_2)
   
    
    print("TESTING CHEBYCHEV PASSED")


##NODI CHEBYCHEV
#RUNGE E NODI CHEBYCHEV

x=np.linspace(-1,1,10)
x_=np.linspace(-1,1,40)

#plt.figure()
plt.scatter(x_,runge(x_))
x_=zeri_T_k(5)
plt.scatter(x_,runge(x_),marker='+',s=50,label='cheby')

y_,coef,y__=coefficenti_indeterminati(x_, runge(x_),False,x)
plt.scatter(x,y__,"cheby inter")
plt.legend()


plt.grid()
plt.title("Chebychev")






