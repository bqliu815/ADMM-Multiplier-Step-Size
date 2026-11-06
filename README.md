# The multiplier step size in classical two-block ADMM: convergence and counterexamples

This repository contains exact numerical certificates and verification code for the paper by Benqi Liu, Qiming Dai, and Zaiwen Wen. The proofs and definitions of the convergence interval and multiplier update appear in the paper.

## Verification

From the repository root, run the following three commands with Python 3.9 or later. They use only the standard library. Keep assertions enabled.

```sh
python3 -I -S -B core/verify.py
python3 -I -S -B core/realizations.py
python3 -I -S -B endpoint/verify.py
```

- `core/verify.py` checks the four finite-memory certificates and the period-48 Gram matrix and support inequalities.
- `core/realizations.py` checks the strict support inequalities used for smooth interpolation and the geometrically growing example.
- `endpoint/verify.py` verifies the refined lower bound by rational polynomial inequalities and the upper bound by a 1024-phase periodic orbit.

For the last check, the certified interval is

```text
1.886882470025 < Gamma_* < 1.886882643762680
```

The [finite-memory and period-48 data guide](core/README.txt) and the [endpoint data guide](endpoint/README.txt) specify how the rational coefficients and Fourier factors are stored.
