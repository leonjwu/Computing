import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import curve_fit
from m1 import flow as fl #assumes p41.f90 has been compiled with: f2py --f90flags='-fopenmp' -c p41.f90 -m m1 -lgomp


def jacobi(n,kmax=10000,tol=1.0e-8,s0=0.1,display=False):
    """ Solve liquid flow model equations with
        jacobi iteration.
        Input:
            n: number of grid points in r and theta
            kmax: max number of iterations
            tol: convergence test parameter
            s0: amplitude of cylinder deformation
            display: if True, plots showing the velocity field and boundary deformation
            are generated
        Output:
            w,deltaw: Final velocity field and |max change in w| each iteration
    """

    #-------------------------------------------
    #Set Numerical parameters and generate grid
    Del_t = 0.5*np.pi/(n+1)
    Del_r = 1.0/(n+1)
    Del_r2 = Del_r**2
    Del_t2 = Del_t**2
    r = np.linspace(0,1,n+2)
    t = np.linspace(0,np.pi/2,n+2) #theta
    tg,rg = np.meshgrid(t,r) # r-theta grid

    #Factors used in update equation (after dividing by gamma)
    rg2 = rg*rg
    fac = 0.5/(rg2*Del_t2 + Del_r2)
    facp = rg2*Del_t2*fac*(1+0.5*Del_r/rg) #alpha_p/gamma
    facm = rg2*Del_t2*fac*(1-0.5*Del_r/rg) #alpha_m/gamma
    fac2 = Del_r2*fac #beta/gamma
    RHS = fac*(rg2*Del_r2*Del_t2) #1/gamma

    #set initial condition/boundary deformation
    w0 = (1-rg**2)/4 #Exact solution when s0=0
    s_bc = s0*np.exp(-10.*((t-np.pi/2)**2))/Del_r
    fac_bc = s_bc/(1+s_bc)

    deltaw = []
    w = w0.copy()
    wnew = w0.copy()

    #Jacobi iteration
    for k in range(kmax):
        #Compute wnew
        wnew[1:-1,1:-1] = RHS[1:-1,1:-1] + w[2:,1:-1]*facp[1:-1,1:-1] + w[:-2,1:-1]*facm[1:-1,1:-1] + (w[1:-1,:-2] + w[1:-1,2:])*fac2[1:-1,1:-1] #Jacobi update

        #Apply boundary conditions
        wnew[:,0] = wnew[:,1] #theta=0
        wnew[:,-1] = wnew[:,-2] #theta=pi/2
        wnew[0,:] = wnew[1,:] #r=0
        wnew[-1,:] = wnew[-2,:]*fac_bc #r=1s

        #Compute delta_p
        deltaw += [np.max(np.abs(w-wnew))]
        w = wnew.copy()
        # if k%1000==0: print("k,dwmax:",k,deltaw[k])
        #check for convergence
        if deltaw[k]<tol:
            # print("Python converged,k=%d,dw_max=%28.16f " %(k,deltaw[k]))
            break

    print(f"Python: k={k}, dw_max={round(deltaw[k])}")
    deltaw = deltaw[:k+1]

    return w,deltaw

if __name__=='__main__':
    #performance()

