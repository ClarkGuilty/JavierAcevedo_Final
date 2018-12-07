#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 09:34:17 2018

@author: Javier Alejandro Acevedo Barroso
"""


import numpy as np
import matplotlib.pyplot as plt

def func(x):
    return (1.0/np.sqrt(2*np.pi))*np.exp(-x*x/2)

def B(proms):
    rta = (proms - proms.mean())**2
    return 1000/(len(proms)-1) *rta.sum()

def W(desvs):
    rta = desvs**2
    return rta.sum() / len(desvs)

procesos = 4

promedios = np.zeros((procesos,1000))
desv = np.zeros((procesos,1000))
V = np.zeros(1000)

for i in range(4):
    datos = np.loadtxt("gauss{:d}.txt".format(i))
    x = np.linspace(datos.min()*1.1, datos.max(),100)
    for j in range(1,1000):
        promedios[i,j] = np.mean(datos[0:j])
        desv[i,j] = np.std(datos[0:j])
    plt.figure()
    plt.hist(datos,bins=50, density=True)
    plt.plot(x,func(x))
    plt.savefig("p{:d}.png".format(i))
        
        
        
for j in range(1,1000):
    V[j] = (len(promedios[:,j])+1)/(len(promedios[:,j])*1000)*B(promedios[:,j]) + (1000-1)/1000 * W(desv[:,j])
    
    
plt.figure()
N = np.arange(1000)
plt.plot(N,V)
plt.savefig("Gelman.png")
        
        
        
        
        