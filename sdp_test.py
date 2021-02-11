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

def random_rho(A=3, rank_lmt=3):
    a_ket = rand_ket(N=A, density=1, dims=None, seed=None)

    a_dm = a_ket * a_ket.dag()

    dms = [a_dm]

    total_dim = A
    for i in range(rank_lmt):
        a_ket = rand_ket(N=A, density=1, dims=None, seed=None)
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
'''
r1 = rho([1, 0, 0])
r2 = rho([-1, 0, 0])
r3 = rho([0, 1, 0])
r4 = rho([0, -1, 0])
r5 = rho([0, 0, 1])
r6 = rho([0, 0, -1])
'''



value1 = 1
value2 = 2

tol = 10e-10


#c_412
r1 = random_rho(A=3, rank_lmt=3)
r2 = random_rho(A=3, rank_lmt=3)
r3 = random_rho(A=3, rank_lmt=3)
r4 = random_rho(A=3, rank_lmt=3)

x = [r1, r2, r3, r4]
x = [picos.Constant(r) for r in x]
print(x)

while value2 - value1 >= tol:

    x = [picos.Constant(r.value) for r in x]

    M12 = [picos.HermitianVariable("M_12_1", 3), picos.HermitianVariable("M_12_2", 3)]
    M13 = [picos.HermitianVariable("M_13_1", 3), picos.HermitianVariable("M_13_2", 3)]
    M14 = [picos.HermitianVariable("M_14_1", 3), picos.HermitianVariable("M_14_2", 3)]
    M23 = [picos.HermitianVariable("M_23_1", 3), picos.HermitianVariable("M_23_2", 3)]
    M24 = [picos.HermitianVariable("M_24_1", 3), picos.HermitianVariable("M_24_2", 3)]
    M34 = [picos.HermitianVariable("M_34_1", 3), picos.HermitianVariable("M_34_2", 3)]
    M = [M12, M13, M14, M23, M24, M34]
    I = picos.Constant(np.eye(3))

    P = picos.Problem()
    P.set_objective("max",  (-x[0] | M14[0]) - (x[0] | M13[0]) - (x[0] | M12[0]) -
                            (x[1] | M24[0]) - (x[1] | M23[0]) + (x[1] | M12[0]) -
                            (x[2] | M34[0]) + (x[2] | M23[0]) + (x[2] | M13[0]) +
                            (x[3] | M34[0]) + (x[3] | M24[0]) + (x[3] | M14[0]))

    P.add_list_of_constraints([m >> 0 for m in M12])
    P.add_list_of_constraints([m >> 0 for m in M13])
    P.add_list_of_constraints([m >> 0 for m in M14])
    P.add_list_of_constraints([m >> 0 for m in M23])
    P.add_list_of_constraints([m >> 0 for m in M24])
    P.add_list_of_constraints([m >> 0 for m in M34])
    P.add_constraint(sum(M12) == I)
    P.add_constraint(sum(M13) == I)
    P.add_constraint(sum(M14) == I)
    P.add_constraint(sum(M23) == I)
    P.add_constraint(sum(M24) == I)
    P.add_constraint(sum(M34) == I)

    P.solve(solver = "cvxopt")

    value1 = picos.value((x[0] | M14[1]) + (x[0] | M13[1]) + (x[0] | M12[1]) + (x[1] | M24[1]) + (x[1] | M23[1]) + (x[1] | M12[0]) + (x[2] | M34[1]) + (x[2] | M23[0]) + (x[2] | M13[0]) + (x[3] | M34[0]) + (x[3] | M24[0]) + (x[3] | M14[0])).real
    print(value1)


    M12 = [picos.Constant(np.array(picos.value(m))) for m in M12]
    M13 = [picos.Constant(np.array(picos.value(m))) for m in M13]
    M14 = [picos.Constant(np.array(picos.value(m))) for m in M14]
    M23 = [picos.Constant(np.array(picos.value(m))) for m in M23]
    M24 = [picos.Constant(np.array(picos.value(m))) for m in M24]
    M34 = [picos.Constant(np.array(picos.value(m))) for m in M34]
    M = [M12, M13, M14, M23, M24, M34]
    x = [picos.HermitianVariable("r{}".format(i), 3) for i in range(1,5)] 

    P = picos.Problem()
    P.set_objective("max",  (-x[0] | M14[0]) - (x[0] | M13[0]) - (x[0] | M12[0]) -
                            (x[1] | M24[0]) - (x[1] | M23[0]) + (x[1] | M12[0]) -
                            (x[2] | M34[0]) + (x[2] | M23[0]) + (x[2] | M13[0]) +
                            (x[3] | M34[0]) + (x[3] | M24[0]) + (x[3] | M14[0]))


    P.add_list_of_constraints([r >> 0 for r in x])
    P.add_list_of_constraints([picos.trace(r) == 1 for r in x])
    #P.add_constraint(x[0] + x[1] == x[2] + x[3])

    P.solve(solver = "cvxopt")

    value2 = picos.value((x[0] | M14[1]) + (x[0] | M13[1]) + (x[0] | M12[1]) + (x[1] | M24[1]) + (x[1] | M23[1]) + (x[1] | M12[0]) + (x[2] | M34[1]) + (x[2] | M23[0]) + (x[2] | M13[0]) + (x[3] | M34[0]) + (x[3] | M24[0]) + (x[3] | M14[0])).real
    print(value2)

    print(value2 - value1)

print("VALUE")
print((-x[0] | M14[0]) - (x[0] | M13[0]) - (x[0] | M12[0]) -
                            (x[1] | M24[0]) - (x[1] | M23[0]) + (x[1] | M12[0]) -
                            (x[2] | M34[0]) + (x[2] | M23[0]) + (x[2] | M13[0]) +
                            (x[3] | M34[0]) + (x[3] | M24[0]) + (x[3] | M14[0]).value.real)

print("4.828427122593918")
print("5.999999999959565") #no equivalences
die

'''
#c_421

r1 = random_state(A=1, rank_lmt=2)
r2 = random_state(A=1, rank_lmt=2)
r3 = random_state(A=1, rank_lmt=2)
r4 = random_state(A=1, rank_lmt=2)
r5 = random_state(A=1, rank_lmt=2)
r6 = random_state(A=1, rank_lmt=2)



x = [r1, r2, r3, r4, r5, r6]
x = [picos.Constant(r) for r in x]

while value2 - value1 >= tol:

    x = [picos.Constant(r.value) for r in x]

    M1 = [picos.HermitianVariable("M_1_1", 2), picos.HermitianVariable("M_1_2", 2), picos.HermitianVariable("M_1_3", 2)]
    M2 = [picos.HermitianVariable("M_2_1", 2), picos.HermitianVariable("M_2_2", 2), picos.HermitianVariable("M_2_3", 2)]
    M3 = [picos.HermitianVariable("M_3_1", 2), picos.HermitianVariable("M_3_2", 2), picos.HermitianVariable("M_3_3", 2)]
    M4 = [picos.HermitianVariable("M_4_1", 2), picos.HermitianVariable("M_4_2", 2), picos.HermitianVariable("M_4_3", 2)]
    M = [M1, M2, M3, M4]
    I = picos.Constant(np.eye(2))

    P = picos.Problem()
    P.set_objective("max", (x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + 
                           (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) +
                           (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0]))

    P.add_list_of_constraints([m >> 0 for m in M1])
    P.add_list_of_constraints([m >> 0 for m in M2])
    P.add_list_of_constraints([m >> 0 for m in M3])
    P.add_list_of_constraints([m >> 0 for m in M4])
    P.add_list_of_constraints([sum(m) == I for m in M])

    P.solve(solver = "cvxopt")

    value1 = picos.value((x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) + (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0])).real
    print(value1)


    M1 = [picos.Constant(np.array(picos.value(m))) for m in M1]
    M2 = [picos.Constant(np.array(picos.value(m))) for m in M2]
    M3 = [picos.Constant(np.array(picos.value(m))) for m in M3]
    M4 = [picos.Constant(np.array(picos.value(m))) for m in M4]
    M = [M1, M2, M3, M4]

    x = [picos.HermitianVariable("r{}".format(i), 2) for i in range(1,7)] 

    P = picos.Problem()
    P.set_objective("max", (x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + 
                           (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) +
                           (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0]))


    P.add_list_of_constraints([r >> 0 for r in x])
    P.add_list_of_constraints([picos.trace(r) == 1 for r in x])
    #P.add_list_of_constraints([x[0] + x[1] == x[2] + x[3], x[0] + x[1] == x[4] + x[5]])


    P.solve(solver = "cvxopt")

    value2 = picos.value((x[0] | M3[2]) + (x[0] | M4[2]) + (x[1] | M2[2]) + (x[1] | M4[1]) + (x[2] | M2[1]) + (x[2] | M3[1]) + (x[3] | M1[2]) + (x[3] | M4[0]) + (x[4] | M1[1]) + (x[4] | M3[0]) + (x[5] | M1[0]) + (x[5] | M2[0])).real
    print(value2)

    print(value2 - value1)





for i in range(6):
    for j in range(i+1 , 6):
        print("{}{}".format(i,j))
        print(np.trace(np.array(x[i].value) @ np.array(x[j].value)))
'''
x_coord = []
y_coord = []
z_coord = []

states = [np.array(r.value) for r in x]


for r in x:
    rho = np.array(r.value)
    #print(rho)
    x = np.trace(s[0] @ rho)
    y = np.trace(s[1] @ rho)
    z = np.trace(s[2] @ rho)
    x_coord.append(x)
    y_coord.append(y)
    z_coord.append(z)
    print("{} {} {}".format(x, y, z))
    print(np.sqrt(x**2 + y**2 + z**2))

m_coords = []
for m in M:
    m_c = []
    for r in m:
        mes = np.array(r.value)
        x = np.trace(s[0] @ mes)
        y = np.trace(s[1] @ mes)
        z = np.trace(s[2] @ mes)
        vec = np.array([x, y, z])
        m_c.append(vec)
    m_coords.append(m_c)




import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D

# Create a sphere
r = 1
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0.0:pi:100j, 0.0:2.0*pi:100j]
x = r*sin(phi)*cos(theta)
y = r*sin(phi)*sin(theta)
z = r*cos(phi)

#Import data
xx = x_coord
yy = y_coord
zz = z_coord

#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.2, linewidth=0)

size=60

ax.scatter(xx[0],yy[0],zz[0],color="blue",s=size)
ax.scatter(xx[1],yy[1],zz[1],color="green",s=size)
ax.scatter(xx[2],yy[2],zz[2],color="yellow",s=size)
ax.scatter(xx[3],yy[3],zz[3],color="purple",s=size)
#ax.scatter(xx[4],yy[4],zz[4],color="orange",s=size)
#ax.scatter(xx[5],yy[5],zz[5],color="black",s=size)

ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("auto")
plt.tight_layout()
plt.show()

'''
print("M1")
#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.2, linewidth=0)
size=60
ax.scatter(xx[3],yy[3],zz[3],color="black",s=size)
ax.scatter(xx[4],yy[4],zz[4],color="black",s=size)
ax.scatter(xx[5],yy[5],zz[5],color="black",s=size)
ax.scatter(m_coords[0][0][0],m_coords[0][0][1],m_coords[0][0][2],color="green",s=size)
ax.scatter(m_coords[0][1][0],m_coords[0][1][1],m_coords[0][1][2],color="green",s=size)
ax.scatter(m_coords[0][2][0],m_coords[0][2][1],m_coords[0][2][2],color="green",s=size)
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("auto")
plt.tight_layout()
plt.show()

print("M2")
#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.2, linewidth=0)
size=60
ax.scatter(xx[1],yy[1],zz[1],color="black",s=size)
ax.scatter(xx[2],yy[2],zz[2],color="black",s=size)
ax.scatter(xx[5],yy[5],zz[5],color="black",s=size)
ax.scatter(m_coords[1][0][0],m_coords[1][0][1],m_coords[1][0][2],color="green",s=size)
ax.scatter(m_coords[1][1][0],m_coords[1][1][1],m_coords[1][1][2],color="green",s=size)
ax.scatter(m_coords[1][2][0],m_coords[1][2][1],m_coords[1][2][2],color="green",s=size)
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("auto")
plt.tight_layout()
plt.show()

print("M3")
#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.2, linewidth=0)
size=60
ax.scatter(xx[0],yy[0],zz[0],color="black",s=size)
ax.scatter(xx[2],yy[2],zz[2],color="black",s=size)
ax.scatter(xx[4],yy[4],zz[4],color="black",s=size)
ax.scatter(m_coords[2][0][0],m_coords[2][0][1],m_coords[2][0][2],color="green",s=size)
ax.scatter(m_coords[2][1][0],m_coords[2][1][1],m_coords[2][1][2],color="green",s=size)
ax.scatter(m_coords[2][2][0],m_coords[2][2][1],m_coords[2][2][2],color="green",s=size)
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("auto")
plt.tight_layout()
plt.show()

print("M4")
#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.2, linewidth=0)
size=60
ax.scatter(xx[0],yy[0],zz[0],color="black",s=size)
ax.scatter(xx[1],yy[1],zz[1],color="black",s=size)
ax.scatter(xx[3],yy[3],zz[3],color="black",s=size)
ax.scatter(m_coords[3][0][0],m_coords[3][0][1],m_coords[3][0][2],color="green",s=size)
ax.scatter(m_coords[3][1][0],m_coords[3][1][1],m_coords[3][1][2],color="green",s=size)
ax.scatter(m_coords[3][2][0],m_coords[3][2][1],m_coords[3][2][2],color="green",s=size)
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("auto")
plt.tight_layout()
plt.show()



'''


print("M13")
#Set colours and render
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(
    x, y, z,  rstride=1, cstride=1, color='c', alpha=0.2, linewidth=0)
size=60
ax.scatter(xx[0],yy[0],zz[0],color="black",s=size)
ax.scatter(xx[2],yy[2],zz[2],color="black",s=size)
ax.scatter(m_coords[1][0][0],m_coords[1][0][1],m_coords[1][0][2],color="green",s=size)
ax.scatter(m_coords[1][1][0],m_coords[1][1][1],m_coords[1][1][2],color="green",s=size)
ax.set_xlim([-1,1])
ax.set_ylim([-1,1])
ax.set_zlim([-1,1])
ax.set_aspect("auto")
plt.tight_layout()
plt.show()



for i in range(len(states)):
    for j in range(i, len(states)):
        print("{}{}".format(i,j))
        print(np.trace(states[i] @ states[j]))