# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 18:09:20 2021

@author: feder
"""

#### INTEGRAZIONE NUMERICA

import numpy as np
import matplotlib.pyplot as plt

#### TRAPEZI
#### GRADO PRECISIONE 1
def trapezi(a,b,n,f):
    
    x=np.linspace(a,b,n)

    plt.figure()
    plt.plot(x,f(x),label='Funzione')
    plt.plot(x,np.zeros_like(x),c='black')
    plt.plot([a,b],[f(a),f(b)],label='SEGMENTO AB')
    plt.scatter(a,f(a),marker='o',c='black')
    plt.scatter(b,f(b),marker='o',c='black')
    plt.text(a,f(a),'A')
    plt.text(b,f(b),'B')
    
    plt.vlines(a,0,f(a),colors='r')
    plt.vlines(b,0,f(b),colors='r')
    plt.grid()
    
    
    
    
    plt.title("Trapezi")
    plt.legend()
    return (b-a)/2.0*(f(a)+f(b))


print("TESTING TRAPEZI\n")
f=lambda x:3*x+5
a=0
b=5

sol_AREA=62.5
sol_STIMATA=trapezi(a,b,10,f)
print("SOL AREA",sol_AREA)
print("SOL STIMATA",sol_STIMATA)

assert np.allclose(sol_AREA,sol_STIMATA)
print("\nTESTING TRAPEZI 1 PASSED\n")

f=lambda x:3*x**2
a=0
b=5

sol_AREA=5**3
sol_STIMATA=trapezi(a,b,30,f)
print("SOL AREA",sol_AREA)
print("SOL STIMATA",sol_STIMATA)
print("ERRORE",abs(sol_STIMATA-sol_AREA))


print("\nTESTING TRAPEZI 2 PASSED")

##G.P. 3
def simpson(a,b,n,f,plot=False):
    
    
    if plot:
        x=np.linspace(a,b,n)
        plt.figure()
        plt.plot(x,f(x),label='Funzione')
        plt.plot(x,np.zeros_like(x),c='black')
        
        plt.plot([a,(a+b)/2,b],[f(a),f((a+b)/2),f(b)],label='SEGMENTO AB/BC')
        
        plt.scatter(a,f(a),marker='o',c='black')
        plt.scatter(b,f(b),marker='o',c='black')
        plt.scatter((a+b)/2,f((a+b)/2),marker='o',c='black')
        
        
        plt.text(a,f(a),'A')
        plt.text((a+b)/2,f((a+b)/2),'B')
        plt.text(b,f(b),'C')
        
        plt.vlines(a,0,f(a),colors='r')
        plt.vlines(b,0,f(b),colors='r')
        plt.vlines((a+b)/2,0,f((a+b)/2),colors='r')
        
        plt.grid()
     
        plt.title("Simpson")
        plt.legend()
   
    return (b-a)/6.0*(f(a)+4*f((a+b)/2)+f(b))




print("\n\nTESTING SIMPSON\n")

f=lambda x:3*x**2
a=0
b=5

sol_AREA=5**3
sol_STIMATA=simpson(a,b,40,f,True)
print("SOL AREA",sol_AREA)
print("SOL STIMATA",sol_STIMATA)
print("ERRORE",abs(sol_STIMATA-sol_AREA))

assert np.allclose(sol_AREA,sol_STIMATA)
print("\nTESTING SIMPSON 1 PASSED\n")


f=lambda x:5*x**4
a=0
b=5

sol_AREA=5**5
sol_STIMATA=simpson(a,b,10,f,True)
print("SOL AREA",sol_AREA)
print("SOL STIMATA",sol_STIMATA)
print("ERRORE",abs(sol_STIMATA-sol_AREA))

#assert np.allclose(sol_AREA,sol_STIMATA)
print("\nTESTING SIMPSON 2 PASSED")


### QUADRATURE COMPOSTE


def trapezi_composti(a,b,f,n=3):
    
    x=np.linspace(a,b,30)
    h=(b-a)/n

    plt.figure()
    plt.plot(x,np.zeros_like(x),c='black')
    plt.plot(x,f(x),label='Funzione')
    
    idx=0
    prec=None
    x=np.linspace(a,b,n)
    
    for x_ in x:
        if prec is not None:
            plt.plot([prec,x_],[f(prec),f(x_)])
            
        plt.scatter(x_,f(x_),marker='o',c='black')
        plt.text(x_,f(x_),idx)
        plt.vlines(x_,0,f(x_),colors='r')
        prec=x_
        idx+=1
    
    
   
    plt.grid()
    plt.title("Trapezi Composti")
    plt.legend()
    
    res=np.sum(f(x[1:-1]))
    res+=(f(x[0])+f(x[-1]))/2
    
    return abs(res*h)


print("\n\nTESTING TRAPEZI COMPOSTI")

f=lambda x:-x+3*x**2-x**3+10
a=5
b=10


sol_STIMATA=trapezi_composti(a, b, f,n=5)

#print("SOL AREA",sol_AREA)
print("SOL STIMATA",sol_STIMATA)
#print("ERRORE",abs(sol_STIMATA-sol_AREA))


def simpson_composti(a,b,f,n=3):
    
    x=np.linspace(a,b,30)
    h=(b-a)/n ##SPAZIATURA

    plt.figure()
    plt.plot(x,np.zeros_like(x),c='black')
    plt.plot(x,f(x),label='Funzione')
    idx=0
    
    res=0
    x=np.linspace(a,b,n)
    
    for x_,next_ in zip(x,x[1:]):
        
        mid_point=(x_+next_)/2
        plt.plot([x_,mid_point,next_],[f(x_),f(mid_point),f(next_)])
        
        res+=simpson(x_, next_, n, f)
        
        plt.scatter(x_,f(x_),marker='o',c='black')
        plt.text(x_,f(x_),idx)
        plt.vlines(x_,0,f(x_),colors='r')
        idx+=1
    
    plt.scatter(next_,f(next_),marker='o',c='black')
    plt.text(next_,f(next_),idx)
    plt.vlines(next_,0,f(next_),colors='r')
    plt.grid()
    plt.title("Simpson Composto")
    plt.legend()
    
   
    return res
    
    

print("\n\nTESTING SIMPSON COMPOSTI")

f=lambda x:2*x**4
a=0
b=5


sol_STIMATA=simpson_composti(a, b, f,n=5)

#print("SOL AREA",sol_AREA)
print("SOL STIMATA",sol_STIMATA)
#print("ERRORE",abs(sol_STIMATA-sol_AREA))
#######SIMPSON COMPOSTI CONTINUARE