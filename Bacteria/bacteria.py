"""MATH96012 2019 Project 1
Leon Wu
01190736
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial
# --------------------------------


def simulate1(N=64, L=8, s0=0.2, r0=1, A=0.2, Nt=100):
    """Part1: Simulate bacterial colony dynamics

    Looping through the Nt time steps, distances between bacteria are calculated in a matrix.
    Then the valid bacteria to include when calculating the summation term in theta
    are stored in a boolean matrix  - the number of neighbours is also computed from this.

    The terms exp(i*theta) can then be multiplied with the boolean matrix to obtain the
    sum of all thetas for bacteria within radius r0, for each bacterium j.

    The angles, velocities, positions and alignment are then calculated using only operations
    on matrices.

    Also, to reduce the number of calculations in for loops, objects such as the
    random array and 1/N are precomputed outside of the loop.
    """

    # Set initial condition
    phi_init = np.random.rand(N) * (2*np.pi)
    r_init = np.sqrt(np.random.rand(N))
    Xinit, Yinit = r_init*np.cos(phi_init), r_init*np.sin(phi_init)
    Xinit += L/2
    Yinit += L/2

    theta_init = np.random.rand(N)*(2*np.pi)  # initial directions of motion
    # ---------------------

    # Initialize position and alignment matrices for output
    X = np.zeros((Nt+1, N))
    Y = np.zeros((Nt+1, N))
    alpha = np.zeros(Nt+1)

    # Precompute scaled Nt x N random matrix
    AexpR = A * np.exp(1j * np.random.rand(Nt, N) * (2 * np.pi))

    # Precompute 1/N
    invN = 1/N

    # Add initial positions and alignment
    X[0, :] = Xinit % L
    Y[0, :] = Yinit % L
    alpha[0] = invN * np.abs(np.sum(np.exp(1j * theta_init), axis=0))

    # Loop through Nt time steps
    for t in range(Nt):

        # Create Nx2 matrix of positions
        positions = np.column_stack((X[t, :], Y[t, :]))

        # Calculate matrix of distances between all bacteria
        distances = scipy.spatial.distance.cdist(positions, positions)

        # Convert to boolean matrix for bacteria within radius r0 (including itself)
        distances_bool = (distances <= r0).astype(int)

        # Count neighbours for each bacterium (including itself)
        neighbours = np.sum(distances_bool, axis=0)

        # Calculate sum of all thetas for bacteria within radius r0 by multiplying
        # boolean matrix with thetas
        sum_exp = np.matmul(np.exp(1j * theta_init), distances_bool)

        # Calculate new theta
        theta_init = np.angle(sum_exp + neighbours * AexpR[t, :])

        # Calculate u and v
        u = s0 * np.cos(theta_init)
        v = s0 * np.sin(theta_init)

        # Calculate new x and y using previous positions and current velocities
        # Loop bacteria outside square to opposite side
        X[t+1, :] = (X[t, :] + u) % L
        Y[t+1, :] = (Y[t, :] + v) % L

        # Calculate measure of alignment, alpha
        alpha[t+1] = invN * np.abs(np.sum(np.exp(1j * theta_init), axis=0))

    return X, Y, alpha


def analyze():
    """Part 2: Add input variables and modify return as needed
    """

    # Set parameters for all plots
    L = 4
    Nt = 400
    A_values = np.arange(0.2, 0.85, 0.05)
    S = 500  # Number of simulations to run


    # Plot 1
    N = 32
    A_index = 6
    A = A_values[A_index]
    _, _, alpha = simulate2(10, A_values, L=L, Nt=Nt, N=N)

    fig = plt.figure(figsize=(10, 6))
    plt.xlabel('t')
    plt.ylabel(f'alpha')
    plt.title(f'10 paths of alpha for N = {N}, A = {round(A, 3)}')
    for i in range(10):
        plt.plot(range(Nt + 1), alpha[:, i, A_index])
    plt.plot(range(Nt + 1), np.mean(alpha[:, :, A_index], axis=1), '--', color='black', markersize=8, label=f'Average path')

    plt.legend(loc='center')
    plt.show()

    N = 16
    # Plot 2
    _, _, alpha = simulate2(S, A_values, L=L, Nt=Nt, N=N)
    fig = plt.figure(figsize=(10, 6))
    plt.xlabel('t')
    plt.ylabel(f'alpha (averaged)')
    plt.title(f'Paths of alpha for N = {N}, averaged over {S} simulations')
    for i in range(len(A_values)):
        plt.plot(range(Nt + 1), np.mean(alpha[:, :, i], axis=1), label=f'A = {round(A_values[i], 2)}')

    plt.legend(loc='center')
    plt.show()


    # Plot 3
    A_values_full = np.arange(0.0, 0.81, 0.01)
    start_index = 40  # Fixed to be after the start of "stationary" behaviour
    N_values = [16, 32]
    var_alpha = dict()
    indices = np.where(A_values_full >= 0.2)
    A_values = A_values_full[indices]
    index = len(A_values_full) - len(A_values)

    fig = plt.figure(figsize=(10, 6))
    plt.xlabel('A')
    plt.ylabel(f'alpha (averaged)')
    plt.title(f'Average values of alpha across {S} simulations and across {start_index} <= t <= {Nt} against A')

    for N in N_values:
        _, _, alpha_full = simulate2(S, A_values_full, L=L, Nt=Nt, N=N)
        alpha = alpha_full[:, :, index:]

        # Average over all stationary values (after start_index) and also simulations
        alpha_average = np.mean(np.mean(alpha[start_index:, :, :], axis=0), axis=0)
        alpha_average_full = np.mean(np.mean(alpha_full[start_index:, :, :], axis=0), axis=0)

        # Calculate variance of alpha across stationary time period, and average across simulations
        var_alpha[N] = np.mean(np.var(alpha[start_index:, :, :], axis=0), axis=0)

        # Plot point of maximum change in alpha wrt A
        diffs = np.zeros(len(A_values) - 1)
        for i in range(len(A_values) - 1):
            diffs[i] = alpha_average[i+1] - alpha_average[i]

        # Find first index of A_values corresponding to biggest change in alpha for a move in A
        # Note A_values are equally spaced
        diff_index = np.argmax(np.abs(diffs))

        # Calculate A* as the value of A between the two A values at which alpha changes the most
        A_star = round((A_values[diff_index+1] + A_values[diff_index])/2, 3)

        plt.plot(A_values, alpha_average, label=f'N = {N}')
        plt.plot(A_star, (alpha_average[diff_index] + alpha_average[diff_index+1])/2,
                 marker='x', markersize=5, linestyle='None', label=f'A* = {A_star}')

    plt.legend(loc="best")
    plt.show()


    # Plot 4
    fig = plt.figure(figsize=(10, 6))
    plt.xlabel('A')
    plt.ylabel(f'var(alpha) (averaged)')
    plt.title(f'Average values of the variance of alpha across {start_index} <= t <= {Nt}, averaged across {S} simulations against A')

    for N in N_values:
        plt.plot(A_values, var_alpha[N], label=f'N = {N}')

    plt.legend(loc="best")
    plt.show()


    # Plot 5, using N = 32 from previous plot
    fig = plt.figure(figsize=(10, 6))
    plt.xlabel('1 - A/A*')
    plt.ylabel(f'alpha (averaged)')
    plt.title(f'Average values of alpha, N = {N} across {S} simulations and across {start_index} <= t <= {Nt} against 1 - A/A*')

    # Restrict A values depending on A_star
    indices2 = np.where(A_values_full < A_star-0.025)
    plt.plot((1 - A_values_full/A_star)[indices2], alpha_average_full[indices2])
    plt.show()

    return None


def simulate2(S, A_values, N=64, L=8, s0=0.2, r0=1, Nt=100):
    """Part 2: Simulation code for Part 2, add input variables and simulation
    code needed for simulations required by analyze

    This function is similar to simulate1 in the steps taken, with added vectorization
    across the number of simulations S, and the values of A in A_values. This results in
    a much faster runtime than looping through simulations and values of A using simulate1.

    To precompute the random term named AexpR, the array of A values is multiplied across
    the whole ndarray using the corresponding axis for the A values.

    This function returns X and Y with dimensions (Nt+1, N, S, len(A_values)) and alpha
    with dimension (Nt+1, S, len(A_values)).
    """

    # Precompute values
    invN = 1/N
    len_A = len(A_values)

    # Set initial condition
    phi_init = np.random.rand(N, S, len_A) * (2*np.pi)
    r_init = np.sqrt(np.random.rand(N, S, len_A))
    Xinit, Yinit = r_init*np.cos(phi_init), r_init*np.sin(phi_init)
    Xinit += L/2
    Yinit += L/2

    theta_init = np.random.rand(N, S, len_A)*(2*np.pi)  # initial directions of motion
    # ---------------------

    # Initialize position and alignment matrices for output
    X = np.zeros((Nt+1, N, S, len_A))
    Y = np.zeros((Nt+1, N, S, len_A))
    alpha = np.zeros((Nt+1, S, len_A))

    # Precompute random numbers, scaled by A values using np.einsum across axis3 (l)
    expR = np.exp(1j * np.random.rand(Nt, N, S, len_A) * (2 * np.pi))
    AexpR = np.einsum('ijkl,l->ijkl', expR, A_values)

    # Add initial positions and alignment
    X[0, :, :, :] = Xinit % L
    Y[0, :, :, :] = Yinit % L
    alpha[0, :, :] = invN * np.abs(np.sum(np.exp(1j * theta_init), axis=0))

    # Loop through Nt time steps
    for t in range(Nt):

        distances = np.zeros((N, N, S, len_A))
        # Calculate matrix of distances between all bacteria within each simulation and A value
        for s in range(S):
            for a in range(len_A):
                positions = np.stack((X[t, :, s, a], Y[t, :, s, a]), axis=1)
                distances[:, :, s, a] = scipy.spatial.distance.cdist(positions, positions)

        # Convert to boolean matrix for bacteria within radius r0 (including itself)
        distances_bool = (distances <= r0).astype(int)

        # Count neighbours for each bacterium (including itself) within each simulation and A value
        neighbours = np.sum(distances_bool, axis=0)

        # Calculate sum of all thetas for bacteria within radius r0 by multiplying
        # boolean matrix with thetas within each simulation and A value

        sum_exp = np.zeros((N, S, len(A_values)), dtype=complex)
        for s in range(S):
            for a in range(len_A):
                sum_exp[:, s, a] = np.matmul(np.exp(1j * theta_init[:, s, a]), distances_bool[:, :, s, a])

        # Calculate new theta
        theta_init = np.angle(sum_exp + neighbours * AexpR[t, :, :, :])

        # Calculate u and v
        u = s0 * np.cos(theta_init)
        v = s0 * np.sin(theta_init)

        # Calculate new x and y using previous positions and current velocities
        # Loop bacteria outside square to opposite side
        X[t+1, :, :, :] = (X[t, :, :, :] + u) % L
        Y[t+1, :, :, :] = (Y[t, :, :, :] + v) % L

        # Calculate measure of alignment
        alpha[t+1, :, :] = invN * np.abs(np.sum(np.exp(1j * theta_init), axis=0))

    return X, Y, alpha


if __name__ == '__main__':
    analyze()
