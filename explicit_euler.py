# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 21:59:12 2023

@author: YHU
"""

# solve diffusion equation
import numpy  as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

#%matplotlib notebook

def explicit_euler(a,b,D,alpha,beta,t):
    C = 0.3 # C < 0.5 to make sure the explicit euler method is stable

    N_space = 20
    dx = (b-a)/N_space  
    # create grid
    #fx = np.zeros(N-1)
    dt = (C*(dx**2))/D
    N_time = (t-0)/dt 

    #x = np.linspace(a,b,N+1)[1:-1]
    #print(x)
    #fx = np.sin((math.pi*(x-a)/(b-a)))
    
    fx = np.zeros(len(N_space))
    
    u = np.zeros((int(N_time)+1,N_space+1))
    
    for i in range(0,N_space):
        u[0,i] = fx[i]

    print(x)
    print(N_time)
    print(u)
    #The outer for loop increments over time (increases n )
    for n in range(int(n_tsteps)):
    
        for i  in range(N-1):
        
            if i == 0:
                u[n+1,0] = u[n,0] + C* (u[n,1]-2*u[n,0]+alpha)
            elif i > 0 and i < N-2:
                u[n+1,i] = u[n,i] + C* (u[n,i+1]-2*u[n,i]+u[n,i-1])
            else:
                u[n+1,N-2] = u[n,N-2] + C* (beta -2*u[n,N-2]+u[n,N-3])
    
    return u
    
    
u = explicit_euler(0,1,0.5,1,0,1)    
print(u)   


#%%
fig,ax = plt.subplots()
ax.set_ylim(0,1)
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')

line, = ax.plot(x,u[0,:])

def animate(i):
    line.set_data((x,u[i,:]))
    return line

ani = animation.FuncAnimation(fig,animate,frames=int(100),blit=True,interval = 100)
plt.show()

#%%
plt.figure()
plt.plot(x,u[4,:])

t =dt* np.arange(n_tsteps) 
print(t)

f = np.exp(-D*((math.pi)**2)*(dt*4)/(b-a)**2)*np.sin(math.pi*(x-a)/(b-a))

plt.plot(x,f)


