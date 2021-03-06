'''
Code from https://arxiv.org/abs/1902.04751 converted from Matlab to Python
'''

from qutip import *
import numpy as np
from scipy.linalg import fractional_matrix_power

def random_ginibre(n, m):
    return (np.random.normal(size=[n, m]) + 1j*np.random.normal(size=[n, m]))/np.sqrt(2)

def random_haar_POVM(d, k, n):
    ret = np.zeros([k, d, d], dtype=complex)
    S = np.zeros([d, d], dtype=complex)
    for i in range(k):
        xi = random_ginibre(d, n)
        wi = xi @ xi.conj().T
        ret[i,:,:] = wi
        S = S + wi

    S = fractional_matrix_power(S, -1/2)

    for i in range(k):
        wi = np.squeeze(ret[i,:,:])
        ret[i,:,:] = S @ wi @ S

    return ret