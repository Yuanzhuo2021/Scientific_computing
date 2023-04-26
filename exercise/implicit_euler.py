# -*- coding: utf-8 -*-
"""
Implicit Euler

@author: YHU
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.sparse as ss 
import matplotlib.animation as animation


#%matplotlib notebook
def explicit_euler(a,b,alpha,beta,D,t,method):
    """
    This function uses explicit euler to solve PDE

    Parameters
    ----------
    a : int/float
        Left boundary 
    b : int/float
        Right boundary 
    alpha : int/float
        Left boundary condition e.g. u(a,0) = alpha 
    beta : int/float
        Right boundary condition e.g. u(b,0) = beta
    D : float/int
        Parameter of the diffusion term 
    t : int/float
        end time solve to
    method:string
        user choose explicit method or RK4 to solve pde
        
    Returns
    -------
    list
        Return a list, the first element is the range of x, the second element is 
        the matrix containing solution u(x,t), the third element is number of steps of 
        time

    """
    
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
    
    x = np.linspace(a,b,N_space+1)
    fx = np.sin(math.pi*(x-a)/(b-a))
    
    u = np.zeros((int(N_time)+1,N_space+1))
    
    # set initial conditions
    for i in range(0,N_space+1):
        u[0,i] = fx[i]
     
    #set boundary conditions
    for k in range(0,int(N_time)+1):
        u[k,0] = alpha
        u[k,N_space] = beta 
        
    # use explicit euler
    if method == 'euler':
        #The outer for loop increments over time (increases n )
        for n in range(int(N_time)):
            # the inner for loop increments over space
            for i in range(N_space-1):
                u[n+1,i+1] = u[n,i+1] + C* (u[n,i+2]-2*u[n,i+1]+u[n,i])
             
    # Use RK4               
    else:
        #The outer for loop increments over time (increases n )
        for n in range(int(N_time)):
            # the inner for loop increments over space
            for i in range(N_space-1):
                 k1 = C*(u[n,i+2]-2*u[n,i+1]+u[n,i])
                 k2 = C*(u[n,i+2]-2*u[n,i+1]+u[n,i]) + dt*D*k1/2
                 k3 = C*(u[n,i+2]-2*u[n,i+1]+u[n,i]) + dt*D*k2/2
                 k4 = C*(u[n,i+2]-2*u[n,i+1]+u[n,i]) + dt*D*k3
                 u[n+1,i+1] = u[n,i+1] + (k1 + 2*k2 + 2*k3 + k4)/6
    
    return [x,u,N_time]




def implicit_euler(a,b,alpha,beta,D,t):
    """
    This function uses implicit euler method to solve pde

    Parameters
    ----------
    a : int
        Left boundary 
    b : int
        Right boundary 
    alpha : float
        Left boundary condition e.g. u(a,0) = alpha 
    beta : float
        Right boundary condition e.g. u(b,0) = beta
    D : float/int
        Parameter in pde
    t : int
        Final time solve to

    Returns
    -------
    list
    The first element in list is a list of x values, the second element is a numpy array of solution u(x,t) of pde

    """
    # set number of steps of space and time
    N_space = 200
    N_time = 1000
    
    #step size
    dx = (b-a)/N_space
    dt = (t-0)/N_time
    
    C = D*dt/(dx**2)

    # identity matrix
    I = ss.identity(N_space-1).toarray()
    
    #stepwise x values
    x = np.linspace(a,b,N_space+1)

    # use numpy to solve the system of equations. There are N-1 equations as there are N-1 unknowns.

    # build A-DD
    A = ss.identity(N_space-1).toarray()
    # the first row and the last row are quite different.We deal with them seperately
    A[0,0]=-2
    A[0,1]=1
    A[N_space-2,N_space-2]=-2
    A[N_space-2,N_space-3]=1
    # for the rest of the row, coefficient for u[i+1] is 1, u[i] is -2,u[i-1] is 1
    for i in range(N_space-3):
        A[i+1,i]=1
        A[i+1,i+1]=-2
        A[i+1,i+2]=1
    
    #build b_DD
    b_DD = np.zeros(N_space-1)
    #put initail and end conditions
    b_DD[0]=alpha
    b_DD[N_space-2]=beta


    # implement the equation
    u0 = np.zeros(N_space-1)
    U = [u0]
    
    for n in range(N_time):
        u1 = ss.linalg.spsolve(I-C*A,u0+C*b_DD)
        U.append(u1)
        u0 = u1 
        
    # add boundary conditions
    U = np.hstack((U,beta*np.ones((N_time+1,1)))) 
    U = np.hstack((alpha*np.ones((N_time+1,1)),U))
    
    return [x,U]


def crank_nicolson(a,b,alpha,beta,D,t):
    """
    This function uses Crank Nicolson method to solve PDE

    Parameters
    ----------
    a : int
        Left boundary 
    b : int
        Right boundary 
    alpha : float
        Left boundary condition e.g. u(a,0) = alpha 
    beta : float
        Right boundary condition e.g. u(b,0) = beta
    D : float/int
        Parameter in pde
    t : int
        Final time solve to

    Returns
    -------
    list
    The first element in list is a list of x values, the second element is a numpy array of solution u(x,t) of pde

    """
    
    # set number of steps of space and time
    N_space = 200
    N_time = 1000
    
    #step size
    dx = (b-a)/N_space
    dt = (t-0)/N_time
    
    C = D*dt/(dx**2)

    # identity matrix
    I = np.identity(N_space-1)
    
    #stepwise x values
    x = np.linspace(a,b,N_space+1)

    # use numpy to solve the system of equations. There are N-1 equations as there are N-1 unknowns.

    # build A-DD
    A = np.zeros((N_space-1,N_space-1))
    # the first row and the last row are quite different.We deal with them seperately
    A[0,0]=-2
    A[0,1]=1
    A[N_space-2,N_space-2]=-2
    A[N_space-2,N_space-3]=1
    # for the rest of the row, coefficient for u[i+1] is 1, u[i] is -2,u[i-1] is 1
    for i in range(N_space-3):
        A[i+1,i]=1
        A[i+1,i+1]=-2
        A[i+1,i+2]=1
    
    #build b_DD
    b_DD = np.zeros(N_space-1)
    #put initail and end conditions
    b_DD[0]=alpha
    b_DD[N_space-2]=beta


    # implement the equation
    u0 = (np.zeros(N_space-1))
    U = [u0]
    
    for n in range(N_time):
        u1 = np.linalg.solve(I-(C/2)*A,np.dot((I+(C/2)*A),u0)+C*b_DD)
        U.append(u1)
        u0 = u1 
        
    # add boundary conditions
    U = np.hstack((U,beta*np.ones((N_time+1,1)))) 
    U = np.hstack((alpha*np.ones((N_time+1,1)),U))
    
    return [x,U]






#%%
if __name__ == '__main__':
    
    a = 0
    b = 1
    alpha = 0
    beta = 0
    D = 0.5
    t = 1
    
    #solve linear diffusion equation
    z = explicit_euler(0,1,0,0,0.5,1,'euler')    
    print(z[1])

   
    plt.plot(z[0],z[1][0])
    plt.plot(z[0],z[1][int(z[2]/10),:])
    plt.plot(z[0],z[1][int(z[2]/5),:])
    plt.plot(z[0],z[1][int(z[2]/2),:])
    plt.plot(z[0],z[1][-1])
    
    plt.legend(['t=0','t=0.1','t=0.2','t=0.5','t=1'])
    plt.xlabel('x')
    plt.ylabel('u(x,t)')


    # ref : blackboard week20demo
    #fig,ax = plt.subplots()
    #ax.set_xlabel('x')
    #ax.set_ylabel('u(x,t)')

    #line, = ax.plot(z[0],z[1][0,:])

    #def animate(i):
        #line.set_data((z[0],z[1][i,:]))
        #return line

    #ani = animation.FuncAnimation(fig,animate,frames = int(z[2]),blit =True,interval = 10)
    #plt.show()

    # test the result
    x = z[0]
    # real solution
    f = np.exp(-D*((math.pi)**2)*(1)/((b-a)**2))*np.sin(math.pi*(x-a)/(b-a))    
    
    test = np.allclose(f,z[1][-1],1e-2)
    print("The solution passes the test: ", test)
    
    
    
    
#%%    
    # solve a PDE using implicit euler method
    z = implicit_euler(0,1,1,0,0.5,1)

    plt.plot(z[0],z[1][0])
    plt.plot(z[0],z[1][200])
    plt.plot(z[0],z[1][400])
    plt.plot(z[0],z[1][999])
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.legend(['t=0','t=0.2','t=0.4','t=1'])


#%%
    # solve a PDE using crank nicolson method
    z = crank_nicolson(0,1,1,0,0.5,1)

    plt.figure()
    plt.plot(z[0],z[1][0])
    plt.plot(z[0],z[1][200])
    plt.plot(z[0],z[1][400])
    plt.plot(z[0],z[1][999])
    plt.xlabel('x')
    plt.ylabel('u(x,t)')
    plt.legend(['t=0','t=0.2','t=0.4','t=1'])


