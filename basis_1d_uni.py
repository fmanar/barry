# Builds the 1D basis spline for a given degree
# Field Manar, 2021-02-09
#
# Given a degree (maximum power) of P each segment has P + 1 coefficients
# to solve for. This requires P + 1 segments, giving (P + 1)**2 total unknowns.
#
# The spline is found by forcing continuity up to the (P - 1)th derivative at
# all segment boundaries (knots). There are (P + 2) knots so this supplies (P +
# 2)*P conditions; one less than the (P + 1)**2 unknowns.
#
# The final condition is to set the overall scale.  In this case by forcing all
# the segments to sum to one at the segment midpoint.  This should be sufficient
# to make them sum to one everywhere.  This constraint is actuall listed first
# because it will be the only one in the P = 0 case.
#
# I wrote this because I wanted to test how one could construct basis splines
# without using the recusive definition universally given.
#
import matplotlib.pyplot as plt
import numpy as np

def index(segment, power, max_power, derivative=0):
    return segment*(max_power + 1) + power + derivative

P = 4

n_unk = (P + 1)**2
n_seg = P + 1
n_coe = P + 1

knots = np.arange(n_seg + 1)

A = np.zeros((n_unk, n_unk))
b = np.zeros((n_unk, 1))

# sum
# first equation expresses magnitude constraint
row = 0
for seg in range(n_seg):
    x = 0.5*(knots[seg] + knots[seg + 1])
    for p in range(n_coe):
        col = index(seg, p, P)
        A[row, col] = x**p
b[row] = 1

# continuity
# consider left endpoint of each segment
# rhs stays all zeros
#   end points have 0 for all derivatives (including 0 derivative)
#   interior points are derivatives of p_left - p_right = 0
row = 1
# poly base has all unit coefficients for plugging into matrix
poly_base = np.polynomial.Polynomial([1,]*n_coe)
for seg in range(0, n_seg + 1):
    poly_der = poly_base
    # repeat for all applicable derivatives
    for der in range(P):
        for p in range(poly_der.degree() + 1 ):
            # add coefficient from previous segment at right side
            if seg != 0:
                col_left = index(seg-1, p, P, der)
                A[row, col_left] = poly_der.coef[p]*knots[seg]**p
            # add coefficient from current segment at left side
            if seg != n_seg:
                col_right = index(seg, p, P, der)
                A[row, col_right] = -poly_der.coef[p]*knots[seg]**p
        row += 1
        poly_der = poly_der.deriv()
coef = np.linalg.solve(A, b)

# plot
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
    

