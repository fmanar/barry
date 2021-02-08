import numpy as np

import barry


l = np.array([0.25, None])
print(l)
barry.fill(l)
print(l)

basis = np.array([[0, 1, 1],
                  [1, 0, 1]])
pos = np.array([[0.75, 0.75]]).T
lam = barry.get_lam(basis, pos)
pos2 = barry.get_pos(basis, lam)
print(basis)
print(pos)
print(lam)
print(pos2)
