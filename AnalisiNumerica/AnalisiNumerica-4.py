# -*- coding: utf-8 -*-
"""
Created on Tue May 11 23:38:41 2021

@author: feder
"""

####MINIMI QUADRATI

import numpy as np
import matplotlib.pyplot as plt

def Minimi_quadrati(A,b):
    
    print("Matrice A:\n",A)
    print("Matrice A^T:\n",A.T)
    print("Matrice A^T*A:\n",np.dot(A.T,A))
    print("Matrice A^T*b:\n",np.dot(A.T,b))
    x=np.dot(np.linalg.inv(np.dot(A.T,A)),np.dot(A.T,b))
    print("x:(Coefficenti)\n",x)
    
    plt.figure()
    plt.scatter(A[:,-2],b)
    plt.grid()
    plt.title("Punti + retta regressione")
    # plt.scatter(A[:,0],A[:,0]*x[1]+x[0])
    points=np.linspace(0,5,25)
    plt.scatter(points,points**3*x[0]+points**2*x[1]+points*x[2]+x[3])
    


A=np.array([[1,-1],[1,1],[2,1]])
b=np.array([[2],[4],[8]])
#Minimi_quadrati(A, b)

##Esempio retta/parabola regressione

x=np.array([[0,0,0,1],[1,1,1,1],[8,4,2,1],[27,9,3,1]])
y=np.array([[1],[3.6],[2.4],[1]])
#y=np.array([[1],[2.1],[2.9],[3.2]])

Minimi_quadrati(x, y)


#x=np.array([[0,1],[1,1],[2,1],[3,1],[4,1],[5,1]])
#y=np.array([[10],[25],[51],[66],[97],[118]])
#Minimi_quadrati(x, y)


####CONTINUARE DA QUIII








