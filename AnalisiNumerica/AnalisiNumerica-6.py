# -*- coding: utf-8 -*-
"""
Created on Wed May 26 23:02:32 2021

@author: feder
"""
####EQUAZIONE NON LINEARI 

####LIBRERIE
import numpy as np
import matplotlib.pyplot as plt


####METODO BISEZIONE
def bisezione(x,y,f,a_,b_,epsilon=0.00001,MAX_ITER=20):
    
    #MOSTRA FUNZIONE
    plt.figure()
    plt.plot(x,y,label='Funz')
    plt.scatter(a_,f(a_),marker='o',c='red')
    plt.text(a_,f(a_)+0.1,"0")
    plt.scatter(b_,f(b_),marker='o',c='red')
    plt.text(b_,f(b_)+0.1,"0")
    plt.plot(np.linspace(0,5,20),np.zeros_like(np.linspace(0,3.6,20)),c='black')
    
    a=a_
    b=b_
    print("Intervallo iniziale:[{},{}] MAX_ITER:{} TOLLERANZA:{}".format(a,b,MAX_ITER,epsilon))

    for i in range(1,MAX_ITER):
        
        #CONDIZIONE BISEZIONE
        assert f(a)*f(b)<0 
        
        #CALCOLO PUNTO MEDIO
        mid_point=(a+b)*0.5    
        mid_point_f=f(mid_point)
        
        plt.scatter(mid_point,mid_point_f,facecolors="none",edgecolors='r',s=40)
        
        ##CONVERGENZA CONDITION
        # if(np.allclose(mid_point_f,0)):
        #     break
        
        if mid_point_f*f(a)<0:
            err_relativo=(mid_point-a)/mid_point
            b=mid_point
            plt.text(b,f(b),str(i))
    
        else:
            err_relativo=(mid_point-b)/mid_point
            a=mid_point
            plt.text(a,f(a),str(i))
    
        #CONVERGENZA < SOGLIA EPSILON
    
        print("Iterazione:{} Intervallo:[{:.3f},{:.3f}], ErroreMAX {:.3f}, Err relativo {:.6f}".format(i,a,b,(b_-a_)/2**(i),abs(err_relativo)))
        if abs(err_relativo)<epsilon:
            break
    
        
    plt.title("Bisezione")
    plt.legend()
    plt.grid(True)
        
    return mid_point
    
    

##TESTING BISEZIONE
print("######### BISEZIONE:\n")

a=1
b=3
punti=20
x=np.linspace(a,b,punti)
f = lambda x:x**3-x**2-8
#f = lambda x:np.e**-x-x
y=f(x)
sol_ZERO=2.39

sol_STIMATA=bisezione(x, y, f, a, b)
print("\nZERO STIMATO:",sol_STIMATA)
print("ZERO EFFETTIVO:",sol_ZERO)

assert np.allclose(sol_STIMATA,sol_ZERO,atol=1.e-2)
print("TEST BISEZIONE PASSED\n")



####METODO CORDE
##UTILITY
##(x,y) e m disegna la retta e da lo zero
def retta_2_punti(m,y_1,x_1,start=0.5,end=3.5):
    #print(m,y_1,x_1)
    x=np.linspace(start,end,50)
    ### X,Y,ZERO
    return x, m*(x-x_1)+y_1, x_1-y_1/m


def metodo_corde(x,y,f,a_,b_,MAX_ITER=20,epsilon=0.001):
    
    
    print("Intervallo iniziale:[{},{}] MAX_ITER:{} TOLLERANZA:{}".format(a_,b_,MAX_ITER,epsilon))

    coef=((f(b_)-f(a_))/(b_-a_))
    
    plt.figure()
    plt.plot(x,y,label='Funz')
    plt.scatter(a_,f(a_),marker='o',c='red')
    plt.text(a_,f(a_),"0")
    plt.scatter(b_,f(b_),marker='o',c='red')
    plt.text(b_,f(b_),"1")
    
    #ASSE X
    plt.plot(np.linspace(0,5,20),np.zeros_like(np.linspace(0,3.6,20)),c='black')
    temp=b_
    
    for i in range(MAX_ITER):              
        x_,y_,zero=retta_2_punti(coef, f(temp), temp,end=b)
        plt.plot(x_,y_)
        plt.scatter(temp,0,marker='o',s=30)
        
        #print(temp-f(temp)/coef)
        prec=temp
        temp=temp-f(temp)/coef
       
        plt.scatter(temp,f(temp),marker='o',c='red')
        plt.text(temp,0,str(i+2))
        plt.text(temp,f(temp),str(i+2))
        
        err_relativo=np.abs((temp-prec)/temp)
        print("Iterazione:{} X_n:{:.3f}, Err relativo:{:.2f}".format(i+1,temp,err_relativo))
        if err_relativo<epsilon:
            break
        
        
    plt.title("Metodo Corde")
    plt.legend()
    plt.grid(True)

    return temp

####TESTING METODO CORDE
print("\n\n######### TESTING METODO CORDE:\n")

a=1
b=4
punti=15
x=np.linspace(a,b,punti)
f = lambda x:x**3-x**2-8
y=f(x)
sol_ZERO=2.39
sol_STIMATA=metodo_corde(x, y, f, a, b)
print("\nZERO STIMATO:",sol_STIMATA)
print("ZERO EFFETTIVO:",sol_ZERO)



assert np.allclose(sol_ZERO,sol_STIMATA,atol=1.e-2)
print("\nTESTING METODO CORDE PASSED")



##REGULA FALSI

def metodo_falsa_posizione(x,y,f,a_,b_,MAX_ITER=15,epsilon=0.001):
    
    
    plt.figure()
    plt.plot(x,y,label='Funz')
    plt.scatter(a_,f(a_),marker='o',c='red')
    plt.text(a_,f(a_),"0")
    plt.scatter(b_,f(b_),marker='o',c='red')
    plt.text(b_,f(b_),"0")
    plt.plot(np.linspace(a,b,20),np.zeros_like(np.linspace(0,3.6,20)),c='black')
    plt.grid()
    
    temp=a_
    temp_2=b_
    coef=((f(b_)-f(a_))/(b_-a_)) ## primo coef
    
    print("Intervallo iniziale:[{},{}] MAX_ITER:{} TOLLERANZA:{}".format(a_,b_,MAX_ITER,epsilon))

    
    for i in range(MAX_ITER):              
        x_,y_,zero=retta_2_punti(coef, f(temp), temp,start=a,end=b)
        plt.plot(x_,y_)
        plt.scatter(temp,0,marker='o',s=30)
        
        prec=temp
        temp=temp-f(temp)/coef ##x_n+1
        ##temp_2==b
        if  (f(temp)*f(b_)<0):
            temp_2=b_
        else:
            temp_2=a_
            
        err_relativo=np.abs((temp-prec)/temp)
        coef=(f(temp)-f(temp_2))/(temp-temp_2) ##coef
        
        plt.scatter(temp,f(temp),marker='o',c='red')
        plt.text(temp,0,str(i+1))
        plt.text(temp,f(temp),str(i+1))
        
        print("Iterazione:{} Intervallo:[{:.4f},{:.4f}] X_n:{:.3f}, Err relativo:{:.4f}".format(i+1,min(temp,temp_2),max(temp_2,temp),temp,err_relativo))
        #print(f(temp),f(temp_2))
        if err_relativo<epsilon:
            break
        
    plt.title("Metodo Falsa Posizione")
    plt.legend()

    return temp
    

print("\n\n######### TESTING METODO FALSA POSIZIONE:\n")

a=0.5
b=3
punti=15
x=np.linspace(a,b,punti)
#f = lambda x:x**3-x**2-8
f = lambda x:np.log(x)*x+x**3-2*x**2
y=f(x)
#sol_ZERO=2.39
sol_ZERO=1.689
sol_STIMATA=metodo_falsa_posizione(x, y, f, a,b,MAX_ITER=20)
print("\nZERO STIMATO:",sol_STIMATA)
print("ZERO EFFETTIVO:",sol_ZERO)

assert np.allclose(sol_ZERO,sol_STIMATA,atol=1.e-2)
print("\nTESTING METODO FALSA POSIZIONE PASSED")



##METODO NEWTON
def newton(x,y,f,f_derivato,a_,b_,MAX_ITER=2,x_0=None,epsilon=0.0001):
    
    plt.figure()
    plt.plot(x,y,label='Funzione')
    plt.plot(np.linspace(-3,11,20),np.zeros_like(np.linspace(0,3.6,20)),c='black')
    
    plt.grid()
    
    if x_0 is None:
        x_0=2 #DEFAULT
    
    lst=[]
    print("Intervallo:[{},{}]Punto iniziale:[{}] MAX_ITER:{} TOLLERANZA:{}".format(a_,b_,x_0,MAX_ITER,epsilon))

    for i in range(MAX_ITER):

        y_plot=r(f_derivato(x_0),x_0,f(x_0),x)
        plt.plot(x,y_plot)
        plt.scatter(x_0,f(x_0),marker='+',s=50)
        plt.text(x_0,f(x_0),s=str(i))
        
        x_=x_0-f(x_0)/f_derivato(x_0)
        plt.text(x_,0,str(i))
        err=abs(x_0-x_)
        x_0=x_
        lst.append(x_)
        
        print("Iterazione:{}  X_n:{:.3f}, Err:{:.4f}".format(i+1,x_,err))
        if err<epsilon:
            break
        
        
    ##print(lst)
    plt.title("Metodo Newton")
    plt.legend()
    #plt.ylim([-50,150])

    return x_    

print("\n\n######### TESTING METODO NEWTON\n")


r=lambda m,x_0,y_0,x:m*(x-x_0)+y_0
f=lambda x: x**3-18*x+3+2*x**2
f_derivato=lambda x:3*x**2-18+4*x
a=0.5
b=11
punti=35
x=np.linspace(a,b,punti)
y=f(x)
sol_ZERO=3.2518
sol_STIMATA=newton(x,y,f,f_derivato,a,b,MAX_ITER=20,x_0=2)

print("\nZERO STIMATO:",sol_STIMATA)
print("ZERO EFFETTIVO:",sol_ZERO)


assert np.allclose(sol_STIMATA,sol_ZERO,atol=1.e-3)
print("\nTESTING METODO NEWTON PASSED")



def secanti(x,y,f,a_,b_,MAX_ITER=10,x_1=None,x_0=None,epsilon=0.0001):
    
    plt.figure()
    plt.plot(x,y,label='Funz')
    plt.plot(np.linspace(0.5,5,20),np.zeros_like(np.linspace(0,3.6,20)),c='black')
    plt.grid()
    
    if x_1 is None:
        x_1=2
    if x_0 is None:
        x_0=3
    
    lst=[]
    print("Punti iniziale:[{},{}] MAX_ITER:{} TOLLERANZA:{}".format(x_0,x_1,MAX_ITER,epsilon))

    for i in range(MAX_ITER):
        
        q_k=(f(x_1)-f(x_0))/(x_1-x_0)
        x=np.linspace(a,b,punti)

        y_plot=r(q_k,x_0,f(x_0),x)
        plt.plot(x,y_plot)
        plt.scatter(x_0,f(x_0),marker='o',s=50)
        plt.text(x_0,f(x_0),s=str(i))
        plt.text(x_1,f(x_1),s=str(i+1))
        
        x_=x_1-f(x_1)/q_k
        x_0=x_1
        err=abs(x_-x_1)
        x_1=x_
        lst.append(x_)
        print("Iterazione:{}  X_n:{:.3f}, Err:{:.4f}".format(i+1,x_,err))
        if err<epsilon:
            break
        
        
    #print(lst)
    plt.title("Metodo Secanti")
    plt.legend()
    plt.ylim([-20,20])
    
    return x_
    

###METODO SECANTI
print("\n\n######### TESTING METODO SECANTI\n")

a=0.5
b=5
punti=15
x=np.linspace(a,b,punti)
#f=lambda x:np.e**x-1
f=lambda x: x**3-18*x+3+2*x**2

y=f(x)
#sol_ZERO=0.0
sol_ZERO=3.2518
sol_STIMATA=secanti(x, y, f, a, b)
print("\nZERO STIMATO:",sol_STIMATA)
print("ZERO EFFETTIVO:",sol_ZERO)


assert np.allclose(sol_ZERO,sol_STIMATA,atol=1.e-3)
print("\nTESTING METODO SECANTI PASSED")




