# -*- coding: utf-8 -*-
"""
Created on Mon May 17 09:40:25 2021

@author: feder
"""

###SPLINES
import numpy as np
import matplotlib.pyplot as plt


def psi_i(x_,i,x):
    
    
    if i==0:
        if x_>=x[0] and x_<=x[1]:
            return (x[1]-x_)/(x[1]-x[0])
        else:
            return 0
        
    if i==len(x):
        if x_>=x[-2] and x_<=x[-1]:
            return (x_-x[-2])/(x[-1]-x[-2])
        else:
            return 0
    
    
    if x_>=x[i-1] and x_<=x[i]:
        return (x_-x[i-1])/(x[i]-x[i-1])
    
    
    if x_>=x[i] and x_<=x[i+1]:
        return (x[i+1]-x_)/(x[i+1]-x[i])

    else:
        return 0

def splines_lin(x_,x,y):
    
    lst=[]
    for x__ in x_:
        s=0    
        for i in range(len(x)):
           s+=y[i]*psi_i(x__, i, x)       
        lst.append(s)
    return np.array(lst)



x=np.linspace(-2,2,25)
y=x**3-2*x**2+4*x


x_=np.linspace(-2,2,35)
y_=(splines_lin(x_, x, y))


plt.figure()
plt.scatter(x,y,label='punti reali',marker='x',s=35)
plt.scatter(x_,y_,label='punti stimati',facecolor='none',edgecolors='r',s=70)
plt.grid()
plt.title("Splines punti")
plt.legend()

plt.figure()
plt.scatter(x,y,label='punti reali',marker='x',s=35)
print(x_.shape,y_.shape)
plt.plot(x_,y_,label='punti stimati')
plt.grid()
plt.title("Splines retta")
plt.legend()


