"""Verify both sides of the refined endpoint enclosure with rational arithmetic.

Usage: python3 -I -S -B verify.py  (or: python3 -S -B verify.py).
No package installation, solver, network connection, or precomputed positivity
partitions are used. Assertions must be enabled.
"""
import sys
from pathlib import Path
from fractions import Fraction as F
from math import lcm
import json

if not __debug__:
    raise RuntimeError('Assertions must be enabled: run without -O.')
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import symbols as q
import gram

START = F(377, 200)
LOWER = F(75475298801, 40000000000)
UPPER = F(47172066094067, 25000000000000)


def positive_on_interval(polynomial):
    """Prove p(x)>0 on [-1,1] by an exhaustive rational Bernstein tree.

    All coefficients in a node share a positive denominator; only their integer
    numerators matter. On bisection the child denominator acquires a factor
    of 2**n. The kth intermediate row of de Casteljau has denominator 2**k
    times the parent denominator, explaining the binary shifts below.
    """
    coefficients = q.bernstein_direct(polynomial, F(-1), F(1))
    den = lcm(*(v.denominator for v in coefficients))
    stack = [([int(v * den) for v in coefficients], 0)]
    n = len(coefficients) - 1
    leaf_count = 0
    while stack:
        c, level = stack.pop()
        if min(c) > 0:
            leaf_count += 1
            continue
        assert min(c[0], c[-1]) > 0, 'Nonpositive endpoint'
        assert level < 28, 'Positivity unresolved at depth limit'
        left, right = [c[0] << n], [c[-1] << n]
        for k in range(1, n + 1):
            c = [a + b for a, b in zip(c, c[1:])]
            left.append(c[0] << (n - k))
            right.append(c[-1] << (n - k))
        stack.extend([(left, level+1), (list(reversed(right)), level+1)])
    return leaf_count


def positive_certificate():
    c = json.loads((HERE / 'positive.json').read_text())
    w = {(b,k): F(v) for b,k,v in c['weights']}
    omega, eps, kappa = (F(c[k]) for k in ('omega', 'epsilon', 'kappa'))
    assert c['lag'] == 256 and len(w) == len(c['weights'])
    assert all(b in ('X','Y') and isinstance(k,int) and 0<abs(k)<=c['lag'] and v>=0
               for (b,k),v in w.items())
    assert min(omega, eps, kappa)>0 and omega+sum(w.values())==1
    assert sum((k if b=='X' else -k)*v for (b,k),v in w.items())==0
    assert list(map(F, c['gammas'])) == [START, LOWER]
    def H(gamma):
        s = q.symbol(w, omega, gamma, eps)
        assert s == q.original_symbol(w, omega, gamma, eps)
        return tuple(q.divide_D(p) for p in s)
    for gamma in [START, LOWER]:
        h = H(gamma)
        a, det, _ = q.positivity_polynomials(h, kappa)
        counts = [positive_on_interval(p) for p in (a, det)]
        print(f'PASS: gamma={gamma}; Bernstein leaf counts A,det: {counts}', flush=True)
    # An explicit extension beyond LOWER, using the coefficient one-norm of
    # the affine derivative of the quotient as a uniform operator norm bound.
    derivative = [q.add(b,q.scale(a,-1)) for a,b in zip(h,H(LOWER+1))]
    m = 1 + sum((2 if i==2 else 1)*abs(v)
                for i,p in enumerate(derivative) for v in p.values())
    delta = min((2-LOWER)/2, kappa/(2*m))
    assert delta>0 and LOWER+delta<2 and delta*m<=kappa/2
    print('PASS: strict right continuation of the lower endpoint', flush=True)


def upper_certificate():
    c = json.loads((HERE / 'periodic.json').read_text())
    assert c['n']==1024 and F(c['gamma'])==UPPER
    r = gram.audit(c)
    assert r['status']=='PASS' and F(r['residual_bounds'][0])>0
    assert len(r['supports'])==2*c['n'] and r['full_support_count']==2*c['n']**2
    # Every original support is affine in gamma with this Gram matrix fixed.
    # One positive delta suffices for the finite set of strict inequalities.
    lower = min(F(s['lower']) for s in r['supports'])
    slope = max(abs(F(t)) for s in r['supports'] for t in s['slope_bounds'])
    delta = min((UPPER-LOWER)/2, lower/(2*(1+slope)))
    assert delta>0 and all(F(s['lower'])>delta*max(abs(F(t)) for t in s['slope_bounds'])
                           for s in r['supports'])
    print(f"PASS: period {c['n']}, {len(c['factors'])} Fourier factors, "
          f"{len(r['supports'])} strict cyclic representatives, positive residual; "
          'strict left continuation of the upper endpoint', flush=True)


if __name__ == '__main__':
    positive_certificate()
    upper_certificate()
    print(f'PASS: {LOWER} < Gamma_* < {UPPER}; width {UPPER-LOWER}')
