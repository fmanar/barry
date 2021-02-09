# builds the 1D basis spline for a given order
# assumes uniform segments
import matplotlib.pyplot as plt
import numpy as np

P = 3

n_unk = (P + 1)**2
n_seg = P + 1
n_coe = P + 1

knots = np.arange(n_seg + 1)

poly_base = np.polynomial.Polynomial([1,]*n_coe)
A = np.zeros((n_unk, n_unk))
b = np.zeros((n_unk, 1))

row = 0

# sum
for i_seg in range(n_seg):
    x = 0.5*(knots[i_seg] + knots[i_seg + 1])
    for p in range(n_coe):
        A[row, i_seg*n_coe + p] = x**p
b[row] = 1
row += 1

# edge conditions
poly_der = poly_base
for i_der in range(P):
    for p in range(poly_der.degree() + 1 ):
        j_left = i_der + p
        j_right = (n_seg - 1)*n_coe + i_der + p
        A[row, j_left] = poly_der.coef[p]*knots[0]**p
        A[row+1, j_right] = poly_der.coef[p]*knots[-1]**p
    row += 2
    poly_der = poly_der.deriv()

# continuity
for i_knt in range(1, len(knots) - 1):
    poly_der = poly_base
    for i_der in range(P):
        for p in range(poly_der.degree() + 1 ):
            j_left = (i_knt - 1)*n_coe + i_der + p 
            j_right = i_knt*n_coe + i_der + p
            A[row, j_left] = poly_der.coef[p]*knots[i_knt]**p
            A[row, j_right] = -poly_der.coef[p]*knots[i_knt]**p
        row += 1
        poly_der = poly_der.deriv()
coef = np.linalg.solve(A, b)

fig, ax = plt.subplots()
for i_seg in range(n_seg):
    coef_local = np.squeeze(coef[(i_seg*n_coe):((i_seg + 1)*n_coe)])
    domain = [knots[i_seg], knots[i_seg + 1]]

    print(f'seg {i_seg}: {coef_local}')

    p = np.polynomial.Polynomial(
        coef=coef_local,
        domain=domain,
        window=domain,
        )
    ax.plot(*p.linspace())
plt.show()
    

