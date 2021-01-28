import picos
import itertools as it
import more_itertools as mit
from pprint import pprint
from qutip import *
import numpy as np
import random



def random_state(A=2, rank_lmt=2):
    a_ket = rand_ket(N=2**A, density=1, dims=None, seed=None)

    a_dm = a_ket * a_ket.dag()

    dms = [a_dm]

    total_dim = 2**A
    for i in range(rank_lmt):
        a_ket = rand_ket(N=2**A, density=1, dims=None, seed=None)
        #print(a_ket)
        #print(np.linalg.norm(a_ket))
        #die
        a_dm = np.array(a_ket.data @ np.conj(a_ket.data).T)
        dms.append(a_dm)

    convex_weights = np.random.normal(size=len(dms))
    convex_weights = convex_weights / np.linalg.norm(convex_weights)
    convex_weights = np.array([x**2 for x in convex_weights])

    total_dm = sum([convex_weights[i] * dms[i] for i in range(len(dms))])

    return np.array(total_dm)


def GME(rho, n):
    N = 2**n
    R = picos.Constant(rho)
    I = picos.Constant(np.eye(N))
    W = picos.HermitianVariable("W", N)
    P_ops = [(picos.HermitianVariable("P_{},{}".format(a, b), N), a) for a, b in mit.set_partitions(range(n), k=2)]
    P = picos.Problem()
    P.set_objective("min", (W | R).real)
    #P.add_constraint(picos.trace(W) == 1)
    P.add_list_of_constraints([P_op >> 0 for P_op, a in P_ops])
    P.add_list_of_constraints([I >> P_op for P_op, a in P_ops])
    P.add_list_of_constraints([(W - P_op).partial_transpose(a) >> 0 for P_op, a in P_ops])
    P.add_list_of_constraints([I >> (W - P_op).partial_transpose(a) for P_op, a in P_ops])
    # Solve the problem.
    P.solve(solver = "cvxopt")
    return (W | R).value
'''
r_1 = picos.HermitianVariable("r1", 4)
r_2 = picos.HermitianVariable("r2", 4)
r_3 = picos.HermitianVariable("r3", 4)


rho_1 = random_state()
rho_2 = random_state()
rho_3 = random_state()
R1 = picos.Constant(rho_1)
#R1 = picos.HermitianVariable("rho1", 4)
R2 = picos.Constant(rho_2)
R3 = picos.Constant(rho_3)
I = picos.Constant(np.eye(4))
observable = [picos.HermitianVariable("M_{}".format(i), 4) for i in range(1, 4)]

P = picos.Problem()
P.set_objective("min", (R1 | observable[0]) + (R2 | observable[1]) + (R3 | observable[2]) )
P.add_list_of_constraints([o >> 0 for o in observable])
P.add_constraint(sum(observable) == I)
P.solve(solver = "cvxopt")

print(observable)

print(np.trace(observable[0].value @ rho_1))
print(np.trace(observable[1].value @ rho_2))
print(np.trace(observable[2].value @ rho_3))
print(np.trace(observable[0].value @ rho_1) + np.trace(observable[1].value @ rho_2) + np.trace(observable[2].value @ rho_3))

obs = [np.array(o.value) for o in observable]

print(obs)
print(sum(obs))
print("linalg")
print(np.linalg.eig(obs[0])[0])
'''



#Pauli matrices
#s = [sigmax(), sigmay(), sigmaz()]
s = {0: sigmax(), 1: sigmay(), 2: sigmaz()}

#General qubit state, input as list of Bloch vector components, i.e. r = [rx, ry, rz]
def rho(r):
    if np.linalg.norm(r) != 1:
        r = np.array(r)/np.linalg.norm(r)
    return np.array((qeye(2) + sum([r[i] * s[i] for i in range(3)])) / 2)

r1 = rho([1, 0, 0])
r2 = rho([-1, 0, 0])
r3 = rho([0, 1, 0])
r4 = rho([0, -1, 0])
r5 = rho([0, 0, 1])
r6 = rho([0, 0, -1])
print(r1)
print(r2)
print(r3)
print(r4)
print(r5)
print(r6)

x = [picos.Constant(r1), picos.Constant(r5), picos.Constant(r2), picos.Constant(r3), picos.Constant(r4), picos.Constant(r6)]
#perms = it.permutations(x)

#best = 0
#best_permutation = x

x[0] = picos.HermitianVariable("r1", 2)
x[1] = picos.HermitianVariable("r2", 2)
x[2] = picos.HermitianVariable("r3", 2)
x[3] = picos.HermitianVariable("r4", 2)
x[4] = picos.HermitianVariable("r5", 2)
x[5] = picos.HermitianVariable("r6", 2)

M1 = [picos.HermitianVariable("M_1_1", 2), picos.HermitianVariable("M_1_2", 2), picos.HermitianVariable("M_1_3", 2)]
M2 = [picos.HermitianVariable("M_2_1", 2), picos.HermitianVariable("M_2_2", 2), picos.HermitianVariable("M_2_3", 2)]
M3 = [picos.HermitianVariable("M_3_1", 2), picos.HermitianVariable("M_3_2", 2), picos.HermitianVariable("M_3_3", 2)]
M4 = [picos.HermitianVariable("M_4_1", 2), picos.HermitianVariable("M_4_2", 2), picos.HermitianVariable("M_4_3", 2)]
I = picos.Constant(np.eye(2))

#i = 0

#for perm in perms:
 #   i += 1
 #   if i % 25 == 0:
  #      print("{} out of 720".format(i))
   # x = [p for p in perm]
P = picos.Problem()
P.set_objective("max", (x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + 
                       (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) +
                       (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0]))

print((x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + 
                       (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) +
                       (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0]))
print(type((x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + 
                       (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) +
                       (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0])))



P.add_list_of_constraints([m >> 0 for m in M1])
P.add_list_of_constraints([m >> 0 for m in M2])
P.add_list_of_constraints([m >> 0 for m in M3])
P.add_list_of_constraints([m >> 0 for m in M4])
P.add_constraint(sum(M1) == I)
P.add_constraint(sum(M2) == I)
P.add_constraint(sum(M3) == I)
P.add_constraint(sum(M4) == I)

P.solve(solver = "cvxopt")

value = picos.value((x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) + (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0])).real
print(value)
#if value > best:
#    best = value
#    best_permutation = perm
#    print(best)
        

#print("best")
#print(best)
#print(best_permutation)

states = [np.array(picos.value(r)) for r in M1]
for r in states:
    print(r)
'''
print(r1)
print(r5)
print(r2)
print(r3)
print(r4)
print(r6)
'''