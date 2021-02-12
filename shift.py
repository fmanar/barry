# construct a spline in barycentric coordinates
# then shift the basis to check the polynomial doesn't change
#
# poly is quadratic, bary condition solved by hand
# y = b0 * l0**0 * l1**0
#   + b1 * l0**1 * l1**0
#   + b2 * l0**0 * l1**1
#   + b3 * l0**2 * l1**0
#   + b4 * l0**1 * l1**1
#   + b5 * l0**0 * l1**2
import matplotlib.pyplot as plt
import numpy as np

import barry

powers = [[0, 0],
        [1, 0],
        [0, 1],
        [2, 0],
        [1, 1],
        [0, 2]]

def poly(l, b):
    s = 0
    for i, p in enumerate(powers):
        s += b[i] * l[0]**p[0] * l[1]**p[1]
    return s

# 1D basis
v0 = np.array([[1, 2]])
v1 = np.array([[0, 3]])
# v0 = np.array([[1, 3]])
# v1 = np.array([[0, 1]])
# v0 = np.array([[2, 0]])
# v1 = np.array([[1, 3]])

# BCs
# a list of [(lambda), derivative, value]
bc_set = [
    { 'l':(1, 0), 'd':0, 'y':1 },
    { 'l':(0.5, 0.5), 'd':0, 'y':0.5 },
    { 'l':(0, 1), 'd':0, 'y':1.5 },
    ]
lhs = np.zeros((6, 6))
rhs = np.zeros((6, 1))

# listed bcs
for i, bc in enumerate(bc_set):
    l = bc['l']
    d = bc['d']
    y = bc['y']
    for j, p in enumerate(powers):
        lhs[i,j] = l[0]**p[0] * l[1]**p[1]
        rhs[i,0] = y

# bary constraints
lhs[3,1] =  v0[0,1]
lhs[3,2] = -v0[0,0]
lhs[4,3] = -v0[0,0]
lhs[4,4] =  v0[0,1]
lhs[5,3] = -v0[0,1]**2
lhs[5,5] =  v0[0,0]**2

print(lhs)
print(rhs)
b = np.linalg.solve(lhs, rhs)
print(b)

x0 = []
y0 = []
x1 = []
y1 = []
for l0 in np.linspace(0, 1):
    l = barry.fill(np.array([l0, None]))
    x0.append(barry.get_pos(v0, l))
    y0.append(poly(l, b))
    x1.append(barry.get_pos(v1, l))
    y1.append(poly(l, b))

fig, ax = plt.subplots()
ax.plot(x0,y0)
ax.plot(x1,y1)
plt.show()

