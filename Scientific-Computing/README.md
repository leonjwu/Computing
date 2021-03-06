# [Project 2: Scientific Computing: Multidimensional SVD, Image Repair Algorithm, Chaotic Systems and Path-Finding Algorithms (Python)](https://github.com/leonwu4951/Computing/blob/master/Scientific-Computing/)

---
---
- Designed multidimensional SVD techniques to reduced dataset storage by 95%

Original Image |  Reduced Image, using 5% of original memory
:-------------------------:|:-------------------------:
![](https://github.com/leonjwu/Computing/blob/master/Scientific-Computing/Full_Image.png)  |  ![](https://github.com/leonjwu/Computing/blob/master/Scientific-Computing/Reduced_Image.png)

<details><summary><big><big><big><b>Click to view code</b></big></big></big></summary>
<p>
    
```python
import numpy as np
import matplotlib.pyplot as plt


def reduce(H, inputs=(12, 69, 12), display=False):
    """
    Construct one or more arrays from H
    that can be used by reconstruct
    Input:
        H: 3-D data array
        inputs: can be used to provide other input as needed
    Output:
        arrays: a tuple containing the arrays produced from H
    """
    M, N, T = H.shape

    if not display:
        p, q, r = inputs

    # Compute p if display
    if display:
        xm, xn, xt = plot_SVD3(H)
        # print(f'xm, xn, xt = {xm}, {xn}, {xt}')
        p = xm #12
    print(f'p={p}')

    # Reduce H, slicing in first dimension
    K = min(N, T)
    UM = np.zeros((M, N, K))
    SM = np.zeros((M, K))
    VhM = np.zeros((M, K, T))
    # Compute SVD for each slice
    for i in range(M):
        UM[i,:,:], SM[i,:], VhM[i,:,:] = np.linalg.svd(H[i,:,:], full_matrices=False)
    # Keep only relevant portions of arrays
    UM = UM[:,:,:p]
    SM = SM[:,:p]
    VhM = VhM[:,:p,:]

    # Reduce UM, slicing in p dimension
    if display:
        x1, x2, x3 = plot_SVD3(UM)
        # print(f'x1, x2, x3 = {x1}, {x2}, {x3}')
        q = x3 #69
    print(f'q={q}')
    K = min(M, N)
    UMU = np.zeros((p, M, K))
    UMS = np.zeros((p, K))
    UMVh = np.zeros((p, K, N))
    # Compute SVD for each slice
    print(UM.shape)
    for i in range(p):
        UMU[i,:,:], UMS[i,:], UMVh[i,:,:] = np.linalg.svd(UM[:,:,i], full_matrices=False)
    # Keep only relevant portions of arrays
    UMU = UMU[:,:,:q]
    UMS = UMS[:,:q]
    UMVh = UMVh[:,:q,:]

    # Reduce VhM, slicing in p dimension
    if display:
        x1, x2, x3 = plot_SVD3(VhM)
        # print(f'x1, x2, x3 = {x1}, {x2}, {x3}')
        r = x2 #12
    print(f'r={r}')
    K = min(M, T)
    VhMU = np.zeros((p, M, K))
    VhMS = np.zeros((p, K))
    VhMVh = np.zeros((p, K, T))
    # Compute SVD for each slice
    for i in range(p):
        VhMU[i,:,:], VhMS[i,:], VhMVh[i,:,:] = np.linalg.svd(VhM[:,i,:], full_matrices=False)
    # Keep only relevant portions of arrays
    VhMU = VhMU[:,:,:r]
    VhMS = VhMS[:,:r]
    VhMVh = VhMVh[:,:r,:]

    if display:
        print("Original Shape:")
        print(H.shape)
        a = H.size
        print(a)
        print("New Shapes:")
        print(UMU.shape, UMS.shape, UMVh.shape, VhMU.shape, VhMS.shape, VhMVh.shape, SM.shape)
        b = np.sum((UMU.size, UMS.size, UMVh.size, VhMU.size, VhMS.size, VhMVh.size, SM.size))
        print(b)
        print("Size of output compared to input:")
        print(f'{round(100*b/a, 3)}%')

    arrays = (UMU, UMS, UMVh, VhMU, VhMS, VhMVh, SM)
    return arrays


def reconstruct(arrays,inputs=()):
    """
    Generate matrix with same shape as H (see reduce above)
    that has some meaningful correspondence to H
    Input:
        arrays: tuple generated by reduce
        inputs: can be used to provide other input as needed
    Output:
        Hnew: a numpy array with the same shape as H
    """
    # Get arrays from reduce function
    UMU, UMS, UMVh, VhMU, VhMS, VhMVh, SM = arrays
    p, M, q = UMU.shape
    N = UMVh.shape[2]
    T = VhMVh.shape[2]

    # Construct UM
    UM = np.zeros((M, N, p))
    for i in range(p):
        UM[:,:,i] = np.matmul(UMU[i,:,:], np.matmul(np.diag(UMS[i,:]), UMVh[i,:,:]))

    # Construct Vh
    VhM = np.zeros((M, p, T))
    for i in range(p):
        VhM[:,i,:] = np.matmul(VhMU[i,:,:], np.matmul(np.diag(VhMS[i,:]), VhMVh[i,:,:]))

    # Construct H
    Hnew = np.zeros((M,N,T))
    for i in range(M):
        Hnew[i,:,:] = np.matmul(UM[i,:,:], np.matmul(np.diag(SM[i,:]), VhM[i,:,:]))

    return Hnew


if __name__=='__main__':
    #reduce(X)
```
</p>
</details>

---
---

- Optimised code for repairing images through vectorisation
- Reformulated matrix computations for most efficient use of Python vectorisation

Broken Image  |  Repaired Image
:-------------------------:|:-------------------------:
![](https://github.com/leonwu4951/Computing/blob/master/Scientific-Computing/Broken_Image.png)  |  ![](https://github.com/leonwu4951/Computing/blob/master/Scientific-Computing/Repaired_Image.png)

<details><summary><big><big><big><b>Click to view code</b></big></big></big></summary>
<p>
    
```python
    
import numpy as np
import matplotlib.pyplot as plt


def repair(R,p,l=1.0,niter=10,inputs=()):
    """Repair corrupted data stored in input
    array, R. Efficient and complete version of repair1.
    Input:
        R: 2-D data array (should be loaded from data1.npy)
        p: dimension parameter
        l: l2-regularization parameter
        niter: maximum number of iterations during optimization
        inputs: can be used to provide other input as needed
    Output:
        A,B: a x p and p x b numpy arrays set during optimization
    """
    #problem setup
    a,b = R.shape
    iK,jK = np.where(R != -1000) #indices for valid data

    #Set initial A,B
    A = np.ones((a,p))
    B = np.ones((p,b))
    # R_mean = np.mean(R)
    # A[:,:] = R_mean
    # B[:,:] = R_mean

    #Create lists of indices used during optimization
    mlist = [[] for i in range(a)]
    nlist = [[] for j in range(b)]

    for i,j in zip(iK,jK):
        mlist[i].append(j)
        nlist[j].append(i)

    dA = np.zeros(niter)
    dB = np.zeros(niter)

    # np.random.seed(1)
    for z in range(niter):
        Aold = A.copy()
        Bold = B.copy()

        #Loop through elements of A and B in different
        #order each optimization step
        for m in np.random.permutation(a):
            for n in np.random.permutation(b):
                if n < p: # Update A[m,n]
                    Bfac = np.sum(B[n, mlist[m]]**2)
                    Rsum = np.dot(A[m, :], B[:, mlist[m]]) - A[m, n]*B[n, mlist[m]]
                    Asum = np.dot(R[m, mlist[m]] - Rsum, B[n, mlist[m]])
                    A[m,n] = Asum/(Bfac+l) #New A[m,n]
                if m < p: # Update B[m,n]
                    Afac = np.sum(A[nlist[n], m]**2)
                    Rsum = np.dot(A[nlist[n], :], B[:, n]) - A[nlist[n], m]*B[m, n]
                    Bsum = np.dot(R[nlist[n], n] - Rsum, A[nlist[n], m])
                    B[m,n] = Bsum/(Afac+l) #New B[m,n]
        dA[z] = np.sum(np.abs(A-Aold))
        dB[z] = np.sum(np.abs(B-Bold))
        if z%10==0: print("z,dA,dB=",z,dA[z],dB[z])

    return A,B


if __name__ == "__main__":
    #repair(X)
```
</p>
</details>

---
---

- Simulated Chaotic Bacteria Interactions through PDEs and Finite Differences
- Implemented banded matrix solvers using SciPy to massively speed up simulations 

Chaos Contours  |  Finite-Difference Scheme Errors
:-------------------------:|:-------------------------:
![](https://github.com/leonwu4951/Computing/blob/master/Scientific-Computing/Chaos.png)  |  ![](https://github.com/leonwu4951/Computing/blob/master/Scientific-Computing/FDerror.png)

<details><summary><big><big><big><b>Click to view code</b></big></big></big></summary>
<p>
    
```python
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
```
</p>
</details>

---
---

- Designed fast path-finding algorithms utilising dequeues, hash tables and binary heaps
- Modified BFS, DFS and Dijkstra algorithms for specific path-finding problems

<details><summary><big><big><big><b>Click to view code</b></big></big></big></summary>
<p>
    
```python
from collections import deque


def flightLegs(Alist,start,dest):
    """
    Find the minimum number of flights required to travel between start and dest,
    and  determine the number of distinct routes between start and dest which
    require this minimum number of flights.
    Input:
        Alist: Adjacency list for airport network
        start, dest: starting and destination airports for journey.
        Airports are numbered from 0 to N-1 (inclusive) where N = len(Alist)
    Output:
        Flights: 2-element list containing:
        [the min number of flights between start and dest, the number of distinct
        jouneys with this min number]
        Return an empty list if no journey exist which connect start and dest
    """

    # Initialization
    L2 = [0 for l in Alist]  # Labels
    L3 = L2.copy()  # Distances when first discovered (also the minimum distance with bfs)
    L4 = L2.copy()  # Number of shortest paths
    Q = deque([start])
    L2[start] = 1
    L4[start] = 1

    # Loop thorugh all connected nodes to start node
    while Q:
        x = Q.popleft() # Remove node from front of queue

        for v in Alist[x]:
            if L2[v] == 0:
                # Add unexplored neighbors to back of queue
                Q.append(v)
                L2[v] = 1
                L3[v] = 1+L3[x]

            # Check whether this is a shortest path to v
            if L3[x]+1 == L3[v]:
                # Add number of shortest paths to x to number of shortest paths to v
                L4[v] += L4[x]

        # Check whether dest was found during iteration of x
        if L2[dest]:
            depth = L3[dest]
            # Iterate through all neighbours of dest and terminate
            for v in Alist[dest]:
                if L3[v]+1 == depth:
                    L4[dest] += L4[v]

            # Subtract duplicate count from node x, which is a neighbour of dest!!!
            L4[dest] -= L4[x]
            return [L3[dest], L4[dest]]

    # dest node not connected to start node
    return []


def safeJourney(Alist,start,dest):
    """
    Find safest journey from station start to dest
    Input:
        Alist: List whose ith element contains list of 2-element tuples. The first element
        of the tuple is a station with a direct connection to station i, and the second element is
        the density for the connection.
    start, dest: starting and destination stations for journey.

    Output:
        Slist: Two element list containing safest journey and safety factor for safest journey
    """

    # Initialize dictionaries
    dinit = float('inf')
    Edict = {} # Explored nodes
    Udict = {} # Unexplored nodes
    N = len(Alist)

    Udict = {i: [None, dinit] for i in range(N)}
    Udict[start][1] = 0

    # Main search
    while Udict:
        # Find node with min d in Udict and move to Edict
        dmin = dinit
        nmin = None
        for n, v in Udict.items():
            d = v[1]
            if d < dmin:
                dmin = d
                nmin = n

        if nmin is None:
            # Destination node is disconnected
            return []

        # Add node to explored
        Edict[nmin] = Udict.pop(nmin)

        if nmin == dest:
            # Solution is complete
            # Route from start to dest
            path = deque([dest])
            node = dest
            while node != start:
                node = Edict[node][0]
                path.appendleft(node)
            return [list(path), Edict[dest][1]]

        # Update provisional safeties for unexplored neighbors of nmin
        for n, d in Alist[nmin]:
            if n in Udict:
                dcomp = max(d, dmin)
                if dcomp < Udict[n][1]:
                    # Replace distance and last node with better one
                    Udict[n] = [nmin, dcomp]

    return []


def shortJourney(Alist,start,dest):
    """
    Find shortest journey from station start to dest. If multiple shortest journeys
    exist, select journey which goes through the smallest number of stations.
    Input:
        Alist: List whose ith element contains list of 2-element tuples. The first element
        of the tuple is a station with a direct connection to station i, and the second element is
        the time for the connection (rounded to the nearest minute).
    start, dest: starting and destination stations for journey.

    Output:
        Slist: Two element list containing shortest journey and duration of shortest journey
    """

    # Initialize dictionaries
    dinit = float('inf')
    Edict = {} # Explored nodes
    Udict = {} # Unexplored nodes
    N = len(Alist)

    Udict = {i: [None, dinit, 1] for i in range(N)}
    Udict[start][1] = 0

    # Main search
    while Udict:
        # Find node with min d in Udict and move to Edict
        dmin = dinit
        nmin = None
        for n, v in Udict.items():
            d = v[1]
            if d < dmin:
                dmin = d
                nmin = n

        if nmin is None:
            # Destination node is disconnected
            return []

        # Add node to explored
        Edict[nmin] = Udict.pop(nmin)

        if nmin == dest:
            # Solution is complete
            # Route from start to dest
            length = Edict[dest][2]
            path = [dest] * length
            node = dest
            for i in range(length-2, -1, -1):
                node = Edict[node][0]
                path[i] = node
            return [path, Edict[dest][1]]

        # Update times distances for unexplored neighbors of nmin
        for n, d in Alist[nmin]:
            if n in Udict:
                dcomp = dmin + d
                if dcomp < Udict[n][1]:
                    # Replace distance and last node with better one
                    Udict[n] = [nmin, dcomp, Edict[nmin][2]+1]
                elif dcomp == Udict[n][1]:
                    # Choose journey with smallest number of steps
                    new_length = Edict[nmin][2]+1
                    old_length = Udict[n][2]
                    if new_length < old_length:
                        Udict[n] = [nmin, dcomp, new_length]

    return []


def cheapCycling(SList,CList):
    """
    Find first and last stations for cheapest cycling trip
    Input:
        Slist: list whose ith element contains cheapest fare for arrival at and
        return from the ith station (stored in a 2-element list or tuple)
        Clist: list whose ith element contains a list of stations which can be
        cycled to directly from station i
    Stations are numbered from 0 to N-1 with N = len(Slist) = len(Clist)
    Output:
        stations: two-element list containing first and last stations of journey
    """

    N = len(CList)
    V = {i for i in range(N)}
    Q = set()
    total_cost = float('inf')
    nodes = [None, None]
    W = dict() # Connections
    counter = -1

    while V:
        # Take any unexplored node
        y = V.pop()
        if CList[y] == []:
            continue

        Q.add(y)
        counter += 1
        W[counter] = set([y])

        # Run until all nodes have been searched
        while Q:
            x = Q.pop() # Remove any node from set

            for v in CList[x]:
                if v not in W[counter]:
                    W[counter].add(v)
                    V.remove(v)
                    Q.add(v) # Add unexplored neighbors to set

    for k, conn_graph in W.items():
        admin = float('inf')
        admin2 = float('inf')
        anmin = None
        anmin2 = None
        rdmin = float('inf')
        rdmin2 = float('inf')
        rnmin = None
        rnmin2 = None
        # Calculate minimum arrival costs
        for node in conn_graph:
            d = SList[node][0]
            if d < admin2:
                if d < admin:
                    admin2 = admin
                    anmin2 = anmin
                    admin = d
                    anmin = node
                else:
                    admin2 = d
                    anmin2 = node

        # Calculate minimum return costs
        for node in conn_graph:
            d = SList[node][1]
            if d < rdmin2:
                if d < rdmin:
                    rdmin2 = rdmin
                    rnmin2 = rnmin
                    rdmin = d
                    rnmin = node
                else:
                    rdmin2 = d
                    rnmin2 = node

        # Compare cost to other connected graphs
        if anmin == rnmin:
            # Arrival and return nodes not distinct
            new_cost1 = admin2 + rdmin
            new_cost2 = admin + rdmin2

            # Code readability favoured over minute efficiency gains with extra if loops
            if new_cost1 < total_cost:
                total_cost = new_cost1
                nodes = [anmin2, rnmin]

            if new_cost2 < total_cost:
                total_cost = new_cost2
                nodes = [anmin, rnmin2]

        else:
            new_cost = admin + rdmin
            if new_cost < total_cost:
                total_cost = new_cost
                nodes = [anmin, rnmin]

    return nodes


if __name__=='__main__':
    pass
    
```
</p>
</details>

---
