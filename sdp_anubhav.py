import numpy as np
import picos
import time

start = time.time()

X = 4 	# six settings for the preparation device
Y = 6	# three settings for the measurement device
K = 2	# binary outcomes
O = 2*Y*(K-1)+1	# total number of operators in our list
conP = [] # an empty list for constraints on the moment matrix level
conM = [] # an empty list for constraints on the substrate level
Prob = [] # a cell for probabilities
S = [] # an empty list for upper bounds
G = [] # cell for X moment matrices

'''
#Positive semidefinite matrices for PSD rank
conPSD = []
A = []
B = []
psd_rank = 2
for y in range(Y):
	A.append([picos.HermitianVariable("A_{}_{}".format(y, i), psd_rank) for i in range(X)])
	B.append([picos.HermitianVariable("B_{}_{}".format(y, k), psd_rank) for k in range(K)])
'''

for x in range(X):
	G.append(picos.HermitianVariable("Moment_matrix_{}".format(x), (O, O))) # declare hermitian SDP variables
	conP.append(G[x] >> 0) # Semi-definiteness contraint
conP.append(G[0] + G[1] == G[2] + G[3]) # preparation equivalences
#conP.append(G[0] + G[1] == G[4] + G[5])


# function to return the position of the operators
def idx(y, k, u):
	return 2*(K-1)*y + 2*k + u + 2 - 1


for x in range(X):
	for y in range(Y):
		for k in range(0, K-1):
			conM.append(G[x][0*O + idx(y,k,0)] == G[x][O * idx(y,k,1) + 0])
			conM.append(G[x][O * idx(y,k,0) + 0] == G[x][0*O + idx(y,k,1)])
	for j in range(O):
		conM.append(G[x][O*j + j] == 1) # unitary constraints

#Measurement equivalences
'''
for x in range(X):
	for j in range(O):
		sum1 = 0
		sum2 = 0
		for y in range(Y):
			for k in range(0, K-1):
				sum1 += G[x][O*j + idx(y,k,0)] + G[x][O*j + idx(y,k,1)]
				sum2 += G[x][O*idx(y,k,0) + j] + G[x][O*idx(y,k,1) + j]
		conM.append(sum1 == 0)
		conM.append(sum2 == 0)
'''

#Obtain observed data
x_l = []
for x in range(X):
	y_l = []
	for y in range(Y):
		k_l = []
		for k in range(0, K-1):
			k_l.append(0.5 + 0.25 * (G[x][0*O + idx(y,k,0)] + G[x][0*O + idx(y,k,1)]))
		y_l.append(k_l)
	x_l.append(y_l)
Prob = x_l

'''
#Enforce psd rank on observed date
for x in range(X):
	for y in range(Y):
		outcome_sum = 0
		for k in range(0, K-1):
			#print(Prob[x][y][k])
			#print(A[y])
			#print(B[y])
			#conPSD.append(Prob[x][y][k] >> (A[y][x] | B[y][k]))
			conPSD.append((A[y][x] | B[y][k]) < Prob[x][y][k].real)
			conPSD.append((A[y][x] | B[y][k]) > Prob[x][y][k].real)
			outcome_sum += Prob[x][y][k]
		#conPSD.append((A[y][x] | B[y][K-1]) < (1 - outcome_sum).real)
		#conPSD.append((A[y][x] | B[y][K-1]) > (1 - outcome_sum).real)


print(conPSD[0])
print(conPSD[1])
'''


S1 = (- Prob[0][2][0] - Prob[0][1][0] - Prob[0][0][0] - Prob[1][4][0] - Prob[1][3][0] + Prob[1][0][0]
	  - Prob[2][5][0] + Prob[2][3][0] + Prob[2][1][0] + Prob[3][5][0] + Prob[3][4][0] + Prob[3][2][0]).real
P = picos.Problem()
P.set_objective("max", S1)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
#P.add_list_of_constraints(conPSD)
P.solve(solver = "cvxopt")
print(S1.value)
S.append(S1)

print(np.max(np.array(G[0].value)))
print(np.min(np.array(G[0].value)))

end = time.time()

print("Time spent: {} s".format(end-start))

'''
S = []

S1 = (Prob[0][0][0] + Prob[2][1][0] + Prob[4][2][0]).real
P = picos.Problem()
P.set_objective("max", S1)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S1.value)
S.append(S1)





S2 = (Prob[0][0][0] + Prob[1][1][0] + Prob[4][2][0]).real
P = picos.Problem()
P.set_objective("max", S2)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S2.value)
S.append(S2)

S3 = (Prob[0][0][0] - Prob[2][0][0] -2 * Prob[4][0][0] - 2 * Prob[1][1][0] + 2 * Prob[2][1][0] + 2 * Prob[4][2][0]).real
P = picos.Problem()
P.set_objective("max", S3)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S3.value)
S.append(S3)

S4 = (2 * Prob[0][0][0] - Prob[1][1][0] + 2 * Prob[2][1][0]).real
P = picos.Problem()
P.set_objective("max", S4)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S4.value)
S.append(S4)

S5 = (Prob[0][0][0] - Prob[4][0][0] + Prob[1][1][0] + Prob[2][1][0] + 2 * Prob[4][2][0]).real
P = picos.Problem()
P.set_objective("max", S5)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S5.value)
S.append(S5)

S6 = (Prob[0][0][0] - Prob[4][0][0] + 2 * Prob[1][1][0] + 2 * Prob[4][2][0]).real
P = picos.Problem()
P.set_objective("max", S6)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S6.value)
S.append(S6)

S7 = (Prob[0][0][0] - Prob[3][0][0] - 2 * Prob[4][0][0] - 2 * Prob[1][1][0] + 2 * Prob[2][1][0] + 2 * Prob[4][2][0]).real
P = picos.Problem()
P.set_objective("max", S7)
P.add_list_of_constraints(conM)
P.add_list_of_constraints(conP)
P.solve(solver = "cvxopt")
print(S7.value)
S.append(S7)
'''