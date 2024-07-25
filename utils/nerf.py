import numpy as np

# Natural Extension of Reference Frame
# Parsons, J., Holmes, J. B., Rojas, J. M., Tsai, J. & Strauss, C. E. M.
# J Comput Chem 26, 1063–1068 (2005).

def NeRF(A, B, C, bond, angl, dihd):

    t = np.pi - angl

    D2 = np.array([bond * np.cos(t),
                   bond * np.cos(dihd) * np.sin(t),
                   bond * np.sin(dihd) * np.sin(t)])
  
    AB = B - A
    BC = C - B
    n_bc = BC / np.linalg.norm(BC)

    n = np.cross(AB, n_bc)
    n = n / np.linalg.norm(n)

    Mx = n_bc
    My = np.cross(n, n_bc)
    #My = My / np.linalg.norm(My)   # Not necessary as the norm should be 1
    Mz = n
    M = np.array([Mx, My, Mz]).T
    # Apply transpose becasue vectors (Mx, My, Mz) have to be the columns of M.

    D = np.matmul(M, D2) + C

    return D

