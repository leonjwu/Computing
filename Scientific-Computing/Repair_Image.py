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
