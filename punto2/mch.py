#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 11:08:10 2018

@author: Javier Alejandro Acevedo Barroso
"""


import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt("datos_observacionales.dat")
treal = datos[:,0]
xreal = datos[:,1]
yreal = datos[:,2]
zreal = datos[:,3]

plt.plot(treal,xreal,label = 'x')
plt.plot(treal,yreal,label = 'y')
plt.plot(treal,zreal,label = 'z')
plt.legend()
plt.xlabel("tiempo")
plt.ylabel("x y z")

t = np.linspace(0, 3,300)
dt = t[1]-t[0]


ii = treal*10
ii = ii.astype(int)*10
ii[-1] -=1
def solveDiff(sigma, rho, beta):
    x = np.linspace(treal.min(), treal.max(), 300)
    y = np.linspace(treal.min(), treal.max(), 300)
    z = np.linspace(treal.min(), treal.max(), 300)
    x[0] = xreal[0]
    y[0] = yreal[0]
    z[0] = zreal[0]
    for i in range(1,300):
        x[i] = x[i-1]+dt*(sigma*(y[i-1]-x[i-1])) 
        y[i] = y[i-1]+dt*( x[i-1]*(rho-z[i-1]) - y[i-1])
        z[i] = z[i-1]+dt*(x[i-1]*y[i-1] - beta*z[i-1] )
    return x,y,z
    
nx,ny,nz = solveDiff(10,2,3)


######         MCH


def loglikelihood(x_obs, y_obs, z_obs, params):
    x,y,z = solveDiff(params[0],params[1],params[2])
    dx = x[ii]-xreal
    dy = y[ii]-yreal
    dz = z[ii]-zreal
    d = -0.5 * np.sum(dx**2+dy**2+dz**2 )
    return d



def logprior(sigma,rho,beta):
    if sigma <30 and rho <30 and beta <30:
        p = 0
    else:
        p = -np.inf
    return p


def hamiltonian(x_obs, y_obs, z_obs, param, param_momentum):
    m = 100.0
    K = 0.5 * np.sum(param_momentum**2)/m
    V = -loglikelihood(x_obs, y_obs, z_obs, param)
    return K + V

def diff_potential(x_obs, y_obs, z_obs, param):
    n_param = len(param)
    div = np.ones(n_param)
    delta = 1E-5
    for i in range(n_param):
        delta_parameter = np.zeros(n_param)
        delta_parameter[i] = delta
        div[i] = loglikelihood(x_obs, y_obs, z_obs, param + delta_parameter) + logprior(param[0],param[1],param[2])
        div[i] = div[i] - loglikelihood(x_obs, y_obs, z_obs, param - delta_parameter)
        div[i] = div[i]/(2.0 * delta)
    return div

#Se actualiza medio dt en p, se actualiza 
def leapfrog_proposal(x_obs, y_obs, z_obs, param, param_momentum):
    N_steps = 5
    delta_t = 1E-3
    m = 10.0
    new_param = param.copy()
    new_param_momentum = param_momentum.copy()
    for i in range(N_steps):
        new_param_momentum = new_param_momentum + diff_potential(x_obs, y_obs, z_obs, param) * 0.5 * delta_t
        new_param = new_param + (new_param_momentum/m) * delta_t
        new_param_momentum = new_param_momentum + diff_potential(x_obs, y_obs, z_obs, param) * 0.5 * delta_t
    new_param_momentum = -new_param_momentum
    return new_param, new_param_momentum


def monte_carlo(x_obs, y_obs, z_obs, N=2000):
    param = [np.random.random(3)]
    param_momentum = [np.random.normal(size=3)]
    for i in range(1,N):
        propuesta_param, propuesta_param_momentum = leapfrog_proposal(x_obs, y_obs, z_obs, param[i-1], param_momentum[i-1])
        energy_new = hamiltonian(x_obs, y_obs, z_obs, propuesta_param, propuesta_param_momentum)
        energy_old = hamiltonian(x_obs, y_obs, z_obs, param[i-1], param_momentum[i-1])
   
        r = min(1,np.exp(-(energy_new - energy_old)))
        alpha = np.random.random()
       # print(r,alpha)
        if(alpha<r):
            param.append(propuesta_param)
        else:
            param.append(param[i-1])
        param_momentum.append(np.random.normal(size=3))    

    param = np.array(param)
    return param

fit = monte_carlo(xreal,yreal,zreal)
best=[]
for i in range(len(fit[0])):
    best.append(np.mean(fit[:,i]))




plt.plot(treal,zreal,label = 'z')
plt.legend()
plt.xlabel("tiempo")
plt.ylabel("x y z")

nx,ny,nz = solveDiff(fit[0],fit[1],fit[2])

plt.figure()
plt.plot(t,nx,label = 'x')
plt.plot(treal,xreal,label = 'x datos')
plt.legend()
plt.figure()
plt.plot(t,ny,label = 'y')
plt.plot(treal,yreal,label = 'y datos')
plt.legend()
plt.figure()
plt.plot(t,nz,label = 'z')
plt.plot(treal,zreal,label = 'z datos')
plt.legend()
    
    
    
    
    
    