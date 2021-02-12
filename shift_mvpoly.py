# construct a spline in barycentric coordinates
# then shift the basis to check the polynomial doesn't change
#
# poly is quadratic, bary condition solved by hand
# re-write to use mvpoly class
#
# old
# y = b0 * l0**0 * l1**0
#   + b1 * l0**1 * l1**0
#   + b2 * l0**0 * l1**1
#   + b3 * l0**2 * l1**0
#   + b4 * l0**1 * l1**1
#   + b5 * l0**0 * l1**2
# new
# y = b0 * l0**0 * l1**0
#   + b1 * l0**0 * l1**1
#   + b2 * l0**0 * l1**2
#   + b3 * l0**1 * l1**0
#   + b4 * l0**1 * l1**1
#   + b6 * l0**2 * l1**0
# map
# b0 -> b0
# b1 -> b3
# b2 -> b1
# b3 -> b6
# b4 -> b4
# b5 -> b2
import matplotlib.pyplot as plt
import numpy as np

import barry

D = 2
P = 2

# 1D basis
v0 = np.array([[1, 2]])
v1 = np.array([[0, 3]])
# v0 = np.array([[1, 3]])
# v1 = np.array([[0, 1]])
# v0 = np.array([[2, 0]])
# v1 = np.array([[1, 3]])

poly = barry.mvpoly.MVPoly(D=D, P=P)

# BCs
# a list of [(lambda), derivative, value]
bc_set = [
    { 'l':(1, 0), 'd':0, 'y':1 },
    { 'l':(0.5, 0.5), 'd':0, 'y':0.5 },
    { 'l':(0, 1), 'd':0, 'y':1.5 },
    ]
n = len(poly)
lhs = np.zeros((n, n))
rhs = np.zeros((n, 1))

# listed bcs
for i, bc in enumerate(bc_set):
    l = bc['l']
    d = bc['d']
    y = bc['y']
    rhs[i,0] = y
    for j in range(len(poly)):
        lhs[i,j] = poly.eval_term(j, l, 1)

# bary constraints
lhs[3,3] =  v0[0,1]
lhs[3,1] = -v0[0,0]
lhs[4,6] = -v0[0,0]
lhs[4,4] =  v0[0,1]
lhs[5,6] = -v0[0,1]**2
lhs[5,2] =  v0[0,0]**2

# remove extra terms (with total degree greater than P)
lhs[6,5] = 1
lhs[7,7] = 1
lhs[8,8] = 1

print(lhs)
print(rhs)
b = np.linalg.solve(lhs, rhs)
print(b)

for i, val in enumerate(b):
    poly.set_coef(i, val)

x0 = []
y0 = []
x1 = []
y1 = []
for l0 in np.linspace(0, 1):
    l = barry.fill(np.array([l0, None]))
    x0.append(barry.get_pos(v0, l))
    y0.append(poly.eval(l))
    x1.append(barry.get_pos(v1, l))
    y1.append(poly.eval(l))

fig, ax = plt.subplots()
ax.plot(x0,y0)
ax.plot(x1,y1)
plt.show()

