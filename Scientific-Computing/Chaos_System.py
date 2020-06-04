import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.linalg import solve_banded
from scipy.sparse import diags
from scipy.spatial.distance import pdist
from scipy.signal import welch
import time


def microbes(phi,kappa,mu,L=1024,Nx=1024,Nt=1201,T=600,display=False):
    """
    Simulate microbe competition model

    Input:
    phi,kappa,mu: model parameters
    Nx: Number of grid points in x
    Nt: Number of time steps
    T: Timespan for simulation is [0,T]
    Display: Function creates contour plot of f when true

    Output:
    f,g: Nt x Nx arrays containing solution
    """

    #generate grid
    x = np.linspace(0,L,Nx)
    dx = x[1]-x[0]
    dx2inv = 1/dx**2

    def RHS(y,t,k,r,phi,dx2inv):
        #RHS of model equations used by odeint

        n = y.size//2

        f = y[:n]
        g = y[n:]

        #Compute 2nd derivatives
        d2f = (f[2:]-2*f[1:-1]+f[:-2])*dx2inv
        d2g = (g[2:]-2*g[1:-1]+g[:-2])*dx2inv

        #Construct RHS
        R = f/(f+phi)
        dfdt = d2f + f[1:-1]*(1-f[1:-1])- R[1:-1]*g[1:-1]
        dgdt = d2g - r*k*g[1:-1] + k*R[1:-1]*g[1:-1]
        dy = np.zeros(2*n)
        dy[1:n-1] = dfdt
        dy[n+1:-1] = dgdt

        #Enforce boundary conditions
        a1,a2 = -4/3,-1/3
        dy[0] = a1*dy[1]+a2*dy[2]
        dy[n-1] = a1*dy[n-2]+a2*dy[n-3]
        dy[n] = a1*dy[n+1]+a2*dy[n+2]
        dy[-1] = a1*dy[-2]+a2*dy[-3]

        return dy


    #Steady states
    rho = mu/kappa
    F = rho*phi/(1-rho)
    G = (1-F)*(F+phi)
    y0 = np.zeros(2*Nx) #initialize signal
    y0[:Nx] = F
    y0[Nx:] = G + 0.01*np.cos(10*np.pi/L*x) + 0.01*np.cos(20*np.pi/L*x)

    t = np.linspace(0,T,Nt)

    #compute solution
    print("running simulation...")
    y = odeint(RHS,y0,t,args=(kappa,rho,phi,dx2inv),rtol=1e-6,atol=1e-6)
    f = y[:,:Nx]
    g = y[:,Nx:]
    print("finished simulation")
    if display:
        plt.figure()
        plt.contour(x,t,f)
        plt.xlabel('x')
        plt.ylabel('t')
        plt.title('Contours of f')

    return f,g


def newdiff(f,h):
    """
    Input:
        f: array whose 2nd derivative will be computed
        h: grid spacing
    Output:
        d2f: second derivative of f computed with compact fd scheme
    """

    N = len(f)

    # Coefficients for compact fd scheme
    alpha = 9/38
    a = (696-1191*alpha)/428
    b = (2454*alpha-294)/535
    c = (1179*alpha-344)/2140

    d = np.array([145/12, -76/3, 29/2, -4/3, 1/12])
    e = d[::-1]
    g = np.array([c/9, b/4, a, -2*c/9-2*b/4-2*a, a, b/4, c/9])

    # Construct banded matrix ab
    ab = np.ones((3, N))
    ab[0, 0] = 0
    ab[0, 1] = 10
    ab[0, 2:] = alpha
    ab[2, :-2] = alpha
    ab[2, -2] = 10
    ab[2, -1] = 0

    # Construct RHS b
    b = np.zeros(N)
    b[0] = np.sum(d*f[:5])
    b[1] = np.sum(g[2:]*f[:5]) + np.sum(g[:2]*f[-3:-1])
    b[2] = np.sum(g[1:]*f[:6]) + g[0]*f[-2]
    b[3:-3] = g[0]*f[:-6] + g[1]*f[1:-5] + g[2]*f[2:-4] \
        + g[3]*f[3:-3] + g[4]*f[4:-2] + g[5]*f[5:-1] + g[6]*f[6:]
    b[-3] = np.sum(g[:-1]*f[-6:]) + g[-1]*f[1]
    b[-2] = np.sum(g[:-2]*f[-5:]) + np.sum(g[-2:]*f[1:3])
    b[-1] = np.sum(e*f[-5:])

    b /= h**2

    # Enable options to enhance performance
    d2f = solve_banded((1, 1), ab, b, overwrite_ab=True, overwrite_b=True, check_finite=False)
    return d2f


def diff(f,h):
    d2f = np.zeros_like(f)
    d2f[0] = f[-2] - 2*f[0] + f[1]
    d2f[1:-1] = f[:-2] -2*f[1:-1] + f[2:]
    d2f[-1] = d2f[0]
    d2f /= h**2
    return d2f


if __name__=='__main__':
    #dymanics()
