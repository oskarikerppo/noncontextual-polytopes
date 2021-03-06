'''
Implement the method of https://arxiv.org/abs/1412.0924 in Python
'''



import random_haar_povm as rand_haar
from qutip import *
import numpy as np
import itertools as it
import picos as pic
import scipy.linalg as la

#print(rand_haar.random_haar_POVM(2,2,1))

def rand_state(dim, rank=1):
	return np.array(rand_dm_ginibre(dim, rank=rank))

#test for exact equality
def arreq_in_list(myarr, list_arrays):
    return next((True for elem in list_arrays if np.array_equal(elem, myarr)), False)

#test for approximate equality (for floating point types)
def arreqclose_in_list(myarr, list_arrays):
    return next((True for elem in list_arrays if elem.size == myarr.size and np.allclose(elem, myarr)), False)

def arreqclose(myarr, somearr):
	return True if np.allclose(myarr, somearr) else False

def ip(A, B):
	return np.trace(A.conj().T @ B)

def proj(v, u):
	return ip(u, v) * u / ip(u, u)

d = 2 #Dimension of quantum systems and measurements
k = 2 #Number of outcome of measurements
m = 6 #number of measurements
n = 3 #Level of relaxation hierarchy
X = 4 #Number of states

state_ranks = 1
mes_ranks = 1

Id = np.eye(d)


#Initialize zero matrix and gamma
T = []
Gamma = np.array([1])
zero_matrix = np.zeros(Gamma.real.shape)
M_basis=[]
elems_in_basis = 0
equals = []
while not arreqclose(Gamma, zero_matrix):
	elems_in_basis += 1
	print(elems_in_basis)
	s1 = rand_state(d, state_ranks)

	M = []
	for i in range(m):
		M.append(rand_haar.random_haar_POVM(d, k, mes_ranks))

	effects = []
	for i in range(m):
		for o in range(k):
			effects.append(M[i][o])


	#create list of operators of corresponding relaxation level n
	
	ops = []
	ops.append(Id)
	for i in range(1, n+1):
		n_level_idx = [x for x in it.product(range(len(effects)), repeat=i)]
		for idx in n_level_idx:
			if idx in equals:
				continue
			n_ops=[]
			for idxx in idx:
				n_ops.append(effects[int(idxx)])
			n_op = n_ops[0]
			for o in range(1, len(n_ops)):
				n_op = n_op @ n_ops[o]
			#if i > 1 and arreqclose_in_list(n_op, ops):
				#print("EQUAL")
				#print(idx)
				#equals.append(idx)
				#pass
			#else:
			ops.append(n_op)


	#print(ops)
	#print(len(ops))

	Gamma = np.zeros([len(ops), len(ops)], dtype=complex)
	for i in range(len(ops)):
		for j in range(i, len(ops)):
			if i == j:
				Gamma[i,j] = np.trace(s1 @ ops[i].conj().T @ ops[j])
			else:
				value = np.trace(s1 @ ops[i].conj().T @ ops[j])
				Gamma[i,j] = value
				Gamma[j, i] = value.conj()

	#print(Gamma.real)
	Gamma = Gamma.real
	zero_matrix = np.zeros(Gamma.real.shape)
	#Next apply Gram-Schmidt process to obtain matrix basis for the space spanned by Gamma matrices
	if len(M_basis) == 0:
		M_basis.append(Gamma)
		T.append(Gamma)
		continue
	else:
		T.append(Gamma)
		Gamma = Gamma - sum([proj(Gamma, M_mat) for M_mat in M_basis])
		if not arreqclose(Gamma, zero_matrix):
			M_basis.append(Gamma)



print(len(M_basis))

for i in range(len(M_basis)):
	M_basis[i] = M_basis[i] / np.sqrt(np.trace(M_basis[i] @ M_basis[i]))


Tm = sum(T)/len(T)

#print(Tm)
#print(Tm.shape)
#print(la.eig(Tm)[0])
eigen_vectors = la.eig(Tm)[1]
#print(eigen_vectors)
U = eigen_vectors.T
UT = la.inv(U)
D = U @ Tm @ UT
#print(D)
print(np.linalg.matrix_rank(D))
for i in range(len(D)):
	for j in range(len(D)):
		if np.isclose(D[i][j], 0):
			D[i][j] = 0
#print(D)
print(np.linalg.matrix_rank(D))
print(np.count_nonzero(D - np.diag(np.diagonal(D))))

V = np.zeros(shape=(np.linalg.matrix_rank(D), Tm.shape[1]))
print(V.shape)
for i in range(np.linalg.matrix_rank(D)):
	if D[i][i] == 0:
		V[i][i] = 0
	else:
		V[i][i] = 1
#print(V)
#print(D)

Vm = V @ U
VmT = Vm.conj().T

#print(Vm)
#print(Vm @ VmT)
#print((Vm @ VmT).shape)



'''
for i in range(len(M_basis)):
	for j in range(len(M_basis)):
		if not np.isclose(ip(M_basis[i], M_basis[j]), 0):
			print("{}{}".format(i,j))
			#print(ip(M_basis[i], M_basis[j]))
			print(np.isclose(ip(M_basis[i], M_basis[j]), 0))
			print(ip(M_basis[i], M_basis[j]))
'''	

#Next define SDP where the GAMMA matrix is defined as real linear combination of M_basis elements
#GAMMA is positive definite and GAMMA_(id, id) = 1

Gx = []
x_vars = []
for i in range(X):
	xx = []
	for j in range(len(M_basis)):
		xx.append(pic.RealVariable("{}-{}".format(i, j)))
	x_vars.append(xx)
	G = np.zeros(shape=M_basis[0].shape)
	for j in range(len(M_basis)):
		G += xx[j] * M_basis[j]
	Gx.append(G)

#Obtain observed data
x_l = []
for x in range(X):
	y_l = []
	for y in range(m):
		k_l = []
		for kk in range(0, k-1):
			k_l.append(Gx[x][(0, 2*y + 1 + kk)])
		y_l.append(k_l)
	x_l.append(y_l)
Prob = x_l


P = pic.Problem()

for i in range(X):
	#P.add_constraint(Gx[i] >> 0)
	P.add_constraint(Vm * Gx[i] * VmT >> 0)
	P.add_constraint(Gx[i][(0,0)] == 1)
	#P.add_constraint(sum(x_vars[i]) == 1)

#Contextuality constraint
#P.add_constraint(Gx[0] + Gx[1] == Gx[2] + Gx[3])


S1 = (- Prob[0][2][0] - Prob[0][1][0] - Prob[0][0][0] - Prob[1][4][0] - Prob[1][3][0] + Prob[1][0][0]
	  - Prob[2][5][0] + Prob[2][3][0] + Prob[2][1][0] + Prob[3][5][0] + Prob[3][4][0] + Prob[3][2][0]).real
P.set_objective("max", S1)

P.solve(solver = "mosek")
print(S1.value)

print("No restrictions")
print("5.999999999959565")
print("Qubit with restriction")
print("4.8284271")
print("Qubit with no restriction")
print("4.8989794")
'''
print("rho1")
for j in range(m):
	print(Prob[0][j][0])
print("rho2")
for j in range(m):
	print(Prob[1][j][0])
print("rho3")
for j in range(m):
	print(Prob[2][j][0])
print("rho4")
for j in range(m):
	print(Prob[3][j][0])
'''


