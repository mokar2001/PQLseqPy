# PQLseqPy

**PQLseqPy** is a fast Python implementation of Penalized Quasi-Likelihood for sequence count data inspired by **PQLseq** (Sun et al., 2019; PMID: 30020412), with added flexibility and significant performance improvements.

---

## 📦 Installation

> **Python** ≥ 3.8 recommended.

#### Conda

```bash
conda install -c bioconda PQLseqPy
```

#### Pip

```bash
pip install PQLseqPy
```

---

## 🧾 Parameters & Attributes 

```python
GLMM(X,Y, K, 
  fixed_tau=None,
  tau2_set_to_zero=False,
  verbose=False,
  starting_step_size=1,
  error_tolerance=1e-5,
  max_iter=200,
  regularization_factor=0)
```

#### Parameters
| Name                    | Type                            | Description                                                                                                                               |
| ----------------------- | ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| `X`                     | `np.ndarray`                    | Covariate matrix. **First column must be ones** (intercept).                                                                              |
| `Y`                     | `np.ndarray`                    | Binomial outcome: `[successes, failures]`.                                                                                                |
| `K`                     | `np.ndarray`                    | Covariance for the random effect.                                                                                                   |
| `fixed_tau`             | `tuple(float, float)` or `None` | If provided `(tau1, tau2)`, variance components are **fixed** during fitting.                                                             |
| `tau2_set_to_zero`      | `bool`                          | If `True` and `fixed_tau is None`, constrain `τ2 = 0` during inference.                                                                   |
| `verbose`               | `bool`                          | Print per-iteration diagnostics.                                                                                                          |
| `starting_step_size`    | `float`                         | Initial step size for Newton/AI updates. Auto-decays every 10 iters.                                                                      |
| `error_tolerance`       | `float`                         | Convergence tolerance (relative change in `β` and `τ`).                                                                                   |
| `max_iter`              | `int`                           | Maximum iterations for the AI loop.                                                                                                       |
| `regularization_factor` | `float`                         | Regularization in matrix inversions. |

#### Attributes (after `.fit()`)

| Name             | Type               | Description                                                             |
| ---------------- | ------------------ | ----------------------------------------------------------------------- |
| `beta`           | `np.ndarray (k,)`  | Estimated fixed-effect coefficients.                                    |
| `se_beta`        | `np.ndarray (k,)`  | Standard errors of `beta`.                                              |
| `z_beta`         | `np.ndarray (k,)`  | z-statistics (`beta / se_beta`).                                        |
| `p_beta`         | `np.ndarray (k,)`  | Two-sided p-values (from `chi2` on `z^2`, df=1).                        |
| `tau`            | `np.ndarray (2,)`  | Estimated `(τ1, τ2)`.                                                   |
| `sigma2`         | `float`            | `τ1 + τ2`.                                                              |
| `h2`             | `float` or `nan`   | Narrow-sense heritability proxy: `τ1 / (τ1 + τ2)` (nan if `σ²=0`).      |
| `variance_model` | `str`              | Diagnostic label for the variance-structure branch used that iteration. |
| `cov`            | `np.ndarray (k,k)` | Approximate covariance matrix of `beta`.                                |
| `converged`      | `bool`             | Convergence flag.                                                       |
| `iter`           | `int`              | Number of iterations used.                                              |
| `elapsed_time`   | `float`            | Seconds elapsed.                                                        |

#### Methods

* `fit() -> GLMM`: fits the model in-place and returns `self`.
* `summary() -> (pandas.Series, pandas.DataFrame)`: returns:

  * `param` (Series): `converged`, `variance_model`, `iter`, `elapsed_time`, `tau1`, `tau2`, `sigma2`, `h2`
  * `estimates` (DataFrame): columns `beta`, `se_beta`, `z_beta`, `p_beta` indexed by covariate (`x1`, `x2`, …)

---

## 📚 Example

#### Code
```python
import numpy as np
from PQLseqPy import GLMM

n = 100
rng = np.random.default_rng(5)
X = np.hstack((np.ones((n,1)), rng.standard_normal((n,3))))
Y = np.hstack((rng.integers(0, 10, (n, 1)), rng.integers(1, 10, (n, 1))))
K = rng.standard_normal((n, 150)); K = K @ K.T

res = GLMM(
    X, Y, K,
    verbose=False,            
    starting_step_size=1.0,  
    error_tolerance=1e-6,     
    max_iter=300              
).fit()

print(res.summary())
```

#### Output
```terminal
converged                   True
variance_model    tau1>0, tau2>0
iter                          18
elapsed_time            0.012182
tau1                    0.001796
tau2                    0.036274
sigma2                  0.038071
h2                      0.047186

GLMM      beta   se_beta    z_beta    p_beta
x1   -0.044621  0.086473 -0.516017  0.605842
x2    0.013953  0.087819  0.158882  0.873762
x3   -0.062288  0.086084 -0.723571  0.469329
x4    0.142272  0.101597  1.400358  0.161406
```

---

## 📖 Citation

If you use this software in academic work, please cite:

* **Akbari, et al. (2024)**. *Pervasive findings of directional selection realize the promise of ancient DNA to elucidate human adaptation.* (PMID: 39314480)

---
## 🛠 Documentation & Development
To view the full API documentation and contributor guide, visit our [Documentation Website](http://example.com).
