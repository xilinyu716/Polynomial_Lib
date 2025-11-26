# Polynomial Lib

Industrial-grade polynomial library with a C++ core (including FFT multiplication) and a Python API. Supports addition, subtraction, multiplication, division, derivative, evaluation, GCD, Sturm sequence for real root counting, real root finding, and Bezout identity.

## Features
- Core operations: `+ - * /`, derivative, evaluate
- Algebraic algorithms: `gcd`, `bezout_identity`
- Real root analysis: `sturm_sequence`, `num_real_roots(a, b)`, `find_real_roots()`
- Performance: hybrid multiplication (`O(n^2)` for small degree, FFT `O(n log n)` above a runtime-configurable threshold)

## Install
```
python -m pip install .
```

## Build & Test (for contributors)
- Requirements: `g++`, `make` (optional: `pybind11`)
- Build shared library: `make lib`
- Build Python extension (optional): `python -m pip install --user pybind11 && make pyext`
- Run all tests: `make test`

## Usage
After installation, import the package and use the Python API:
```python
from polylib import Polynomial, set_fft_threshold

p1 = Polynomial([1, -3, 2])
p2 = Polynomial([1, -1])
print(p1 + p2)
q, r = p1 / p2
print(q, r)
set_fft_threshold(128)
```



