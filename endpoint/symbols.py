"""Rational Laurent-polynomial arithmetic for the finite-lag multiplier."""
from fractions import Fraction as F
from math import comb


def clean(p):
    return {k: F(v) for k, v in p.items() if v}


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, F(0)) + v
    return clean(out)


def scale(p, a):
    return clean({k: a*v for k, v in p.items()})


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i+j] = out.get(i+j, F(0)) + a*b
    return clean(out)


def star(p):
    return {-k: v for k, v in p.items()}


def real(p):
    return scale(add(p, star(p)), F(1, 2))


ONE = {0: F(1)}
D = {-1: F(-1), 0: F(2), 1: F(-1)}
d = {0: F(1), -1: F(-1)}


def support(weights, omega, block):
    p = {0: F(omega)}
    for (b, lag), val in weights.items():
        if b == block:
            p = add(p, {0: val, lag: -val})
    return clean(p)


def symbol(weights, omega, gamma, epsilon=F(0)):
    """Closed formula, entries (11,22,12)."""
    lx, ly = (support(weights, omega, b) for b in ('X', 'Y'))
    b = {0: F(1), -1: F(gamma)-1}
    aa = real(mul(star(d), lx))
    bb = add(real(mul(mul(star(b), d), lx)), scale(D, -epsilon))
    cc = scale(add(mul(D, lx), mul(b, add(star(lx), scale(star(ly), -1)))), F(1, 2))
    return aa, bb, cc


def original_symbol(weights, omega, gamma, epsilon=F(0)):
    """Independent construction from the original normal/point rows."""
    b = {0: F(1), -1: F(gamma)-1}
    normals = {'X': [scale(d, -1), scale(b, -1)], 'Y': [{}, b]}
    points = {'X': [ONE, d], 'Y': [ONE, {}]}
    mat = [[{}, {}], [{}, {}]]
    for block in ('X', 'Y'):
        ell = support(weights, omega, block)
        for i in range(2):
            for j in range(2):
                v = mul(mul(star(normals[block][i]), ell), points[block][j])
                adj = star(mul(mul(star(normals[block][j]), ell), points[block][i]))
                mat[i][j] = add(mat[i][j], scale(add(v, adj), F(-1, 2)))
    mat[1][1] = add(mat[1][1], scale(D, -epsilon))
    assert mat[1][0] == star(mat[0][1])
    return mat[0][0], mat[1][1], mat[0][1]


def divide_D(p):
    """Exact Laurent long division; nonzero remainder is an error."""
    if not p:
        return {}
    remainder = dict(p)
    lowest = min(p)
    quotient = {}
    while remainder and max(remainder) >= lowest + 2:
        top = max(remainder)
        coeff = -remainder[top]
        quotient[top-1] = coeff
        remainder = add(remainder, scale(mul({top-1: coeff}, D), -1))
    assert not remainder, f'nonzero quotient remainder: {remainder}'
    assert mul(quotient, D) == clean(p)
    return clean(quotient)


def poly_add(a, b):
    c = [F(0)]*max(len(a), len(b))
    for p in (a, b):
        for k, v in enumerate(p):
            c[k] += v
    return trim(c)


def trim(a):
    while len(a) > 1 and not a[-1]:
        a.pop()
    return a or [F(0)]


def poly_scale(a, s):
    return trim([s*v for v in a])


def poly_mul(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            c[i+j] += v*w
    return trim(c)


def chebyshev(n):
    ts = [[F(1)], [F(0), F(1)]]
    us = [[F(1)], [F(0), F(2)]]
    while len(ts) <= n:
        for seq in (ts, us):
            seq.append(poly_add([F(0)]+poly_scale(seq[-1], 2), poly_scale(seq[-2], -1)))
    return ts[:n+1], us[:n+1]


def real_imag_polynomials(p):
    n = max(map(abs, p), default=0)
    ts, us = chebyshev(n)
    c, e = [p.get(0, F(0))], [F(0)]
    for k in range(1, n+1):
        plus, minus = p.get(k, F(0)), p.get(-k, F(0))
        c = poly_add(c, poly_scale(ts[k], plus+minus))
        e = poly_add(e, poly_scale(us[k-1], plus-minus))
    return c, e


def positivity_polynomials(quotients, kappa):
    aa, bb, cc = quotients
    assert aa == star(aa) and bb == star(bb)
    a = real_imag_polynomials(add(aa, {0: -kappa}))[0]
    b = real_imag_polynomials(add(bb, {0: -kappa}))[0]
    c, e = real_imag_polynomials(cc)
    det = poly_add(poly_mul(a, b), poly_scale(poly_mul(c, c), -1))
    det = poly_add(det, poly_scale(poly_mul([F(1), F(0), F(-1)], poly_mul(e, e)), -1))
    return a, det, (b, c, e)


def eval_poly(p, x):
    out = F(0)
    for c in reversed(p):
        out = out*x + c
    return out


def bernstein_direct(p, left, right):
    n = len(p)-1
    width = right-left
    powers = [sum(p[k]*comb(k,j)*left**(k-j)*width**j for k in range(j,n+1))
              for j in range(n+1)]
    return [sum(powers[j]*F(comb(k,j),comb(n,j)) for j in range(k+1))
            for k in range(n+1)]


