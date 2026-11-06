# The multiplier step size in classical two-block ADMM: convergence and counterexamples

Exact rational certificates and Python verification code for the paper by Benqi Liu, Qiming Dai, and Zaiwen Wen.

With a fixed penalty $\beta>0$, the factor $\gamma$ scales only the multiplier update $\lambda^{k+1}=\lambda^k+\gamma\beta r^{k+1}$, where $r^{k+1}$ is the constraint residual. The uniform convergence interval starts at zero and gives a guarantee for every well-defined exact trajectory in the problem class of Theorem 2.1. Its right endpoint is $\Gamma_\ast$, defined in (2.6). The proof of Corollary 2.3 uses the certificates below to establish the strict bounds

$$
1.886882470025 < \Gamma_\ast < 1.886882643762680.
$$

The exact value of $\Gamma_\ast$ and convergence at $\gamma=\Gamma_\ast$ remain open.

## Run the checks

Download and unpack the [version cited in the paper](https://github.com/bqliu815/ADMM-Multiplier-Step-Size/archive/dbc85ba3e257a389e95106c2b3c50431ce1ef9d1.zip). From the repository root, run the following commands with Python 3.9 or later. They use only the standard library and require assertions to be enabled.

```sh
python3 -I -S -B core/verify.py
python3 -I -S -B core/realizations.py
python3 -I -S -B endpoint/verify.py
```

| Program | What it checks |
| --- | --- |
| [`core/verify.py`](core/verify.py) | The four finite-memory certificates in Table 1, and the Gram matrix and support inequalities for the exact period-48 orbit on two fixed polytopes (Proposition 4.3). |
| [`core/realizations.py`](core/realizations.py) | The strict support and interpolation inequalities for the smooth strongly convex example, and the recurrences for geometric growth with piecewise-quadratic objectives (Theorem 2.2). |
| [`endpoint/verify.py`](endpoint/verify.py) | For the lower bound, the program checks polynomial positivity on $[-1,1]$ using Bernstein coefficients. For the upper bound, it checks strict support inequalities and a positive residual norm for a separate 1024-phase orbit (Appendix C). |

On success, the first two programs report `"status": "PASS"`. The last prints the two exact rational bounds on `Gamma_*` and their difference on its final `PASS` line.

The [finite-memory and period-48 guide](core/README.txt) and the [endpoint guide](endpoint/README.txt) describe the data formats and arithmetic checks. Formal proofs for the finite-memory analysis and period-48 construction are also available in [ReasBook (Lean 4.32.0)](https://optpku.github.io/ReasBook/versions/v4.32.0/papers/glowinski_admm_counterexample/).
