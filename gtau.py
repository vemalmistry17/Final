import matplotlib.pyplot as plot
from random import *
import numpy
from math import *
file = open('gtau0.txt','w+') 
#Function to find exponential of matrix
def expmatrix(M):
	#Eigenvalues and vectors of the Hamiltonian matrix
	eigval, eigvec = numpy.linalg.eig(M)
	#Vector of exponentials of eigenvalues
	expvec = []	
	for i in range(0, len(n)):
		expvec.append(exp(eigval.item(i)))
	#Turning the exponential eigenvalues into a matrix
	diagmatrix = numpy.diag(numpy.array(expvec))
	#Unitary transformation on diagonal exponential matrix 
	ematrix = eigvec * diagmatrix *\
	numpy.transpose(numpy.conjugate(eigvec))
	return(ematrix)

#Function to calculate the weight 
def weight(y,spin):
	#B is a list of each of the 7 potential(U) matricies
	B = []
	for r in range(0, len(y)):
		#Calculate the potential matrix for each spin
		d = numpy.array([[exp(spin*lamda*y[r]),0,0,0],\
		[0,1,0,0],[0,0,1,0],[0,0,0,1]])
		Pot = numpy.matrix(d)
		jj = KE*Pot
		B.append(jj)
	z = numpy.zeros((4,4))
	I = numpy.eye(4)
	#7 by 7 block matrix with 4 by 4 blocks. 28 by 28 in total
	x = numpy.block([[I,z,z,z,z,z,B[6]],
	[-B[0],I,z,z,z,z,z],
	[z,-B[1],I,z,z,z,z],
	[z,z,-B[2],I,z,z,z],
	[z,z,z,-B[3],I,z,z],
	[z,z,z,z,-B[4],I,z],
	[z,z,z,z,z,-B[5],I],])
	#Determinant of entire x matrix
	o = numpy.linalg.det(x)
	#Inverse of entire x matrix
	inv = numpy.linalg.inv(x)
	return(o,inv)
	
#Parameters that we'll define
V = 0.5
t = 1
epsilon = 0
beta = 5
l = 7
deltatau = beta/l
P = 0
lamda = acosh(exp(deltatau*P/2))

xaxis = []
vecgreenup = []

# Hamiltonian matrix
n = numpy.array([[-epsilon,V,0,0],\
[V,0,t,0],[0,t,0,t],[0,0,t,0]])
H0 = deltatau*numpy.matrix(n) 
KE = expmatrix(H0)
a = [1,1,1,1,1,1,1]
b = [0,0,0,0,0,0,0]
greenup = numpy.zeros((4*len(a),4*len(a)))
greendown = numpy.zeros((4*len(a),4*len(a)))

# k is the number of iteration of the spin configuration
for k in range(0, 5000):
	for j in range(0, len(a)):
		for m in range(0, len(a)):
			#Each element of a is copied into b
			b[m]=a[m]
		#Single spin is flipped
		b[j] = -b[j]
		#Weight calculated for old and new 
		#spin configuration
		aweight = weight(a,1)[0]*weight(a,-1)[0] 
		bweight = weight(b,1)[0]*weight(b,-1)[0]
		greenup = greenup + weight(a,1)[1]
		greendown = greendown + weight(a,-1)[1]
		ran = random()
		#If ratio of weights is greater than some 
		#random number, accept b as a 
		if bweight/aweight > ran:
			for m in range(0, len(a)):
				a[m]=b[m]
	print(k)

#Average of each matrix
avggreenup = greenup/(len(a)*(k+1))
pos = 0
#Green's function
for val in range(pos,len(a)*len(KE),4):
	xaxis.append(val/28)
	vecgreenup.append(2*(avggreenup[val,0])-1)

#Write the result to a file	
for item in xaxis,vecgreenup: 
	file.write('%s\n' % item)