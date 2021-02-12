# a class to hold multivariate polynomials
#
# Roche uses the terms:
#   - partial degree in x_i: highest exponent of variable x_i
#   - max degree: greatest partial degree (highest exponent anywhere)
#   - total degree: greatest sum of exponenets in a single term
#
# a polynomial has D input variables and total degree P
#
# multivariate polynomials are stored in one of three ways:
#   - dense: use a ndarray with D dimensions of length P
#   - recursive: collect terms recursively and store in a recursive list.
#       Effectively a tree.
#   - sparse: store only non-zero terms as a list of (coefficient, exponent)
#       tuples
#
# for the present purposes, we are primarily concerned with only evaluation and
# storage, not manipulation beyond differentiation. The polynomials considered
# are of low degree, total degree, and (probably?) dense. Since we are
# also all about solving coefficients an explicit linear numbering system for
# the terms is very useful.
#
# for us, D ~< 4, P ~< 5
#   - dense: 5**4 = 625 coefficients
#   - recusive: (~0.5)*(5**4) = 312 coefficients 
#   - sparse: ~= recusive
#
# The total degree requirement means that a dense representation will waste
# space. But the ordering requirement doesn't play well with the recursive
# requirement.  The density points away from sparse representation.
#
# store recusively, access densely?
#
# perhaps I'm getting ahead of myself.  the dense rep will waste space, and this
# will be a problem for large splines, but lets ignore it for now.
# indexing is really easy with ravel
import numpy as np

class MVPoly:
    def __init__(self, coef=None, D=None, P=None):
        if coef is None:
            self.D = D
            self.P = P
            self.shape = (P+1,)*D
            self.coef = np.zeros(self.shape)
        else:
            self.D = len(coef.shape)
            self.P = coef.shape[0] - 1
            self.shape = coef.shape
            self.coef = coef

    def get_exponent(self, i):
        return np.unravel_index(i, self.shape)

    def set_coef(self, i, val):
        ij = np.unravel_index(i, self.shape)
        self.coef[ij] = val

    def get_coef(self, i):
        ij = np.unravel_index(i, self.shape)
        return self.coef[ij]

    def eval(self, x):
        val = 0
        for i, c in enumerate(self.coef.flatten()):
            val += self.eval_term(i, x, c)
        return val

    def eval_term(self, i, x, c):
        if len(x) != self.D:
            raise ValueError("x is not of correct dimension")
        e = self.get_exponent(i)
        val = c
        for _x, _e in zip(x, e):
            val *= _x**_e
        return val

    def __len__(self):
        return self.coef.size

    def __str__(self):
        string = 'y = \n'
        for ic, c in enumerate(self.coef.flatten()):
            string += f'{ic}:   + {c}'
            for ix, e in enumerate(self.get_exponent(ic)):
                string += f' * x[{ix}]**{e}'
            string += '\n'    
        return string
