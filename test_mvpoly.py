import numpy as np

import barry.mvpoly as bmvp

D = 2
P = 2
coef = np.ones((P+1,)*D)

poly = bmvp.MVPoly(coef)
print(f'D {poly.D}')
print(f'P {poly.P}')
print(poly.coef)
print(poly)

x = (0.5,)*D
y = poly.eval(x)
print(f'f({x}) = {y}')