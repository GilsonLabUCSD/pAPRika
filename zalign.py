import numpy as np

m1 = np.asarray([-2.1,3.5,-6.2])
m2 = np.asarray([1.4,6.7,-2.6])

trans = m1*-1.0

m1 = m1 + trans
m2 = m2 + trans

z = np.asarray([0.0,0.0,1.0])

x = np.cross(m2,z)/np.linalg.norm(np.cross(m2,z))

th = np.arccos(np.dot(m2,z)/(np.linalg.norm(m2)*np.linalg.norm(z)))

I = np.identity(3)
A = np.asarray([ [0,-1.0*x[2],x[1]], [x[2],0,-1.0*x[0]], [-1.0*x[1],x[0],0] ])
R = I + np.dot(np.sin(th),A) + np.dot((1.0 - np.cos(th)),np.dot(A,A))

print np.dot(R,m2)

