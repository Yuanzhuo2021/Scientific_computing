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
    
    x = np.linspace(a,b,N+1)
    fx = 4*x*(3-x)
    
    
    u = np.zeros((int(N_time)+1,N_space+1))
    
    for i in range(0,N_space):
        u[0,i] = fx[i]
        
    for k in range(0,int(N_time)):
        u[k,0] = alpha
        u[k,N_space] = beta 

    print(x)
    print(N_time)
    print(u)
    #The outer for loop increments over time (increases n )
    for n in range(int(N_time)):
    
        for i  in range(N_space-1):
        
            if i == 0:
                u[n+1,1] = u[n,1] + C* (u[n,2]-2*u[n,1]+alpha)
            elif i > 0 and i < N_space:
                u[n+1,i+1] = u[n,i+1] + C* (u[n,i+2]-2*u[n,i+1]+u[n,i])
            else:
                u[n+1,N_space-1] = u[n,N_space-1] + C* (beta -2*u[n,N-1]+u[n,N_space-2])
    
    
    #fig,ax = plt.subplots()

    #ax.set_xlabel('x')
    #ax.set_ylabel('u(x,t)')

    #line, = ax.plot(x,u[0,:])

    #def animate(i):
        #line.set_data((x,u[i,:]))
      #  return line

    #ani = animation.FuncAnimation(fig,animate,frames=int(N_time),blit=False,interval = 100)
   # plt.show()
    
    return [x,u,N_time]

z = explicit_euler(0,3,1,0,0,3)    
print(z[1])   
print(len(z[1]))

plt.plot(z[0],z[1][0,:])
plt.plot(z[0],z[1][100,:])
plt.plot(z[0],z[1][200,:])

#%%


#%%
plt.figure()
plt.plot(x,u[4,:])

t =dt* np.arange(n_tsteps) 
print(t)

f = np.exp(-D*((math.pi)**2)*(dt*4)/(b-a)**2)*np.sin(math.pi*(x-a)/(b-a))

plt.plot(x,f)


