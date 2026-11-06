Finite-memory convergence certificates and the period-48 counterexample
=======================================================================

From this directory, run:

    python3 -I -S -B verify.py
    python3 -I -S -B realizations.py

Python 3.9 or later is sufficient. Use the standard library only and do not
enable Python's -O or -OO mode. Both commands print PASS only after the
exact arithmetic checks have succeeded. To retain all rational values, use
the optional argument --output filename.json on either command.

positive.json supplies the four intervals, their common endpoint weights,
and the independent entries for the storage blocks in the appendix on
rational finite-memory certificates. In each case:

* Order the lag rows by distance d=1,...,ell. Their four integer columns
  are the backward and forward weights for F, then for J. Divide by
  lag_denominator.
* In each block, enumerate the flow indices (i,j) in lexicographic order,
  -1<=i<=ell-1, 0<=j<=ell-1, i!=j. Divide flow_numerators by
  flow_denominator. The missing flows (i,-1) follow from flow balance.
* Pack the upper triangle of each 2*ell square principal block row by
  row, including the diagonal. The two blocks, with their own denominators,
  refer to the left and right endpoints, respectively. The constant-history
  kernel identities determine the two omitted rows and columns of storage.

verify.py reconstructs both finite matrix inequalities, the kernel
identities and the interpolation weights. It checks 224 strict rational
LDL pivots at the eight interval endpoints. Coverage between endpoints
uses the parameter-affinity proof in the article.

cycle48.json contains the period-48 cyclic Gram input. p and q contain
lags 0,...,24 of symmetric circulants; c contains all 48 cross lags (the row c^K in the paper).
gamma0, sigma, and rank_one_coefficient specify the rank-one correction,
and gamma_hat and active_slope are the threshold and slope for this
particular cycle. These 98 independent integer entries are also given in
the article. The verifier uses fraction-free Bareiss elimination and checks
the 96 affine support representatives, equivalent by cyclic symmetry to
all 4608 support inequalities. Their slopes, their endpoint values and the
nonzero residual certify a nonconvergent orbit as explained in the text.

realizations.py independently reconstructs the strict cycle, the smooth
strongly convex interpolation data, and a geometrically growing trajectory
with finite piecewise-quadratic objectives. It prints a short summary by
default; --output includes strict support margins, interpolation constants,
and growth coefficients. The regularity and infinite-trajectory conclusions
follow from the analytic arguments in the article.

Notation: chi in the article is the geometric growth factor (kappa in the
realization script, reported as growth_factor). It is distinct from the
positivity margin kappa in the separate unit-circle certificate. The F
and J subscripts in the finite-memory data denote the two block labels.
