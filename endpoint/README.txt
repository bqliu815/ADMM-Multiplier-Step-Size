Exact verification of the refined endpoint enclosure
===================================================

Run from this directory with Python 3.9 or later:

    python3 -I -S -B verify.py

Nothing needs to be installed. The verifier reads two JSON files, writes
only to standard output, and exits on the first failed assertion. Do not
use Python's -O flag. A successful run ends in

    PASS: 75475298801/40000000000 < Gamma_* < 47172066094067/25000000000000; width 2171721/12500000000000

The algebraic checks do not search for a multiplier or solve an optimization
problem. The analytic implications for ADMM are proved in the article.

1. Uniform convergence
----------------------

positive.json supplies ONE multiplier for the entire interval from 377/200
to 75475298801/40000000000. The main finite-memory convergence theorem
covers the range up to 377/200. Each entry [B,k,v] is a support multiplier
of nonnegative rational weight v in block B=X or Y. Here k<0 means a past
point; for k>0 the comparison uses a future point, and is shifted backward
by k steps in the causal supply. The common origin coefficient is omega.
The weights obey the normalization and zero signed first moment:

  omega + sum v = 1,    sum_X k v - sum_Y k v = 0.

For a unit complex z put D=2-z-z^(-1), d=1-z^(-1), b=1+(gamma-1)z^(-1),
and ell_B=omega+sum_{k in B} v(1-z^k). The normal rows are (-d,-b),(0,b);
the corresponding point rows are (1,d),(1,0), in coordinates (eta,w).
Subtract epsilon*D*diag(0,1) from minus the Hermitian part of the sum of
normal-row adjoint times ell_B times point row. This gives the raw supply
symbol Q. symbols.py constructs Q in two algebraically equivalent ways;
verify.py checks agreement and divisibility of every entry by D in the
Laurent ring. This division defines H=Q/D also at z=1 by continuity.

Writing x=cos(theta), the first diagonal element and determinant of
H-kappa*I are real polynomials in x. A positive minimum Bernstein
coefficient is sufficient for strict positivity on the entire domain of
a Bernstein expansion. An exhaustive dyadic subdivision of [-1,1] proves
this at both endpoint parameters (Sylvester's criterion). In order of
increasing gamma, the numbers of terminal subintervals for the diagonal
and determinant polynomials are (7,9) and (7,24). No floating-point
samples or saved subdivision outcomes are used. The MATRIX is affine in
gamma, so positivity at both endpoints covers the entire interval. A
uniform bound on its derivative gives a positive extension to the right.
The bridge from a strict full-circle supply to convergence of actual
initialized ADMM trajectories is analytical; a polynomial positivity
check alone does not replace it.

2. Failure of uniform convergence
---------------------------------

periodic.json gives n=1024, the rational test factor gamma, a common
integer denominator d0, and 75 rational Fourier factors. An entry with
mode j and integers (a,b,c,e) represents the complex coordinate pair
u=(a+ib)/d0, v=(c+ie)/d0. Define one coordinate block of eta_t and w_t
from u exp(2 pi i j t/n) and v exp(2 pi i j t/n) by taking real/imaginary
parts and multiplying by sqrt(m_j/n), where m_j=2 except m_{n/2}=1.
The Nyquist factors are real. Concatenate blocks from all the factors.
This gives a real cyclic Gram matrix that is positive semidefinite by
construction. The vectors themselves need not have rational coordinates.

Set r_t=w_t-w_{t-1}, x_t=eta_t+r_t, a_t=eta_{t-1}-x_t-gamma*w_{t-1},
and c_t=r_t+gamma*w_{t-1}. The polytopes are the convex hulls of the
origin and the x_t, or of the origin and the eta_t, respectively.
For every nontrivial displacement k verify positivity of the inner
products <a_t,x_t-x_{t+k}> and <c_t,eta_t-eta_{t+k}>; the origin
comparisons are <a_t,x_t> and <c_t,eta_t>. Cyclic symmetry reduces
2 n^2 comparisons to 2 n=2048 representatives. gram.py computes strict
rational lower bounds for each representative and for ||r_0||^2. Its
outward-rounded half-angle bounds for roots of unity use only integer
arithmetic (intervals.py). The projection-cycle criterion in the paper
therefore realizes an actual nonconvergent ADMM orbit. Since the finite
support margins are positive and affine in gamma for this same tuple,
the failure persists on a left neighborhood of the upper test factor.

These two sufficient tests localize a uniform initial convergence
endpoint. They do not determine its exact value or endpoint attainment.

Manuscript conventions: the block labels X,Y in positive.json denote F,J
(the tilted functions are tilde F and J). The undamped quotient is H_gamma;
the positive matrix whose diagonal and determinant are tested is
T_gamma = H_gamma - kappa*I. The time discount alpha is one in these tests.
