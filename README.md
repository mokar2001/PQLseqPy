# PQLseqPy

**PQLseqPy** is a fast Python implementation of Penalized Quasi-Likelihood for sequence count data inspired by **PQLseq** (Sun et al., 2019; PMID: 30020412), with added flexibility and significant performance improvements.

---

## 📦 Installation

> **Python** ≥ 3.8 recommended.

Conda:

```bash
conda install -c bioconda PQLseqPy
```

pip:

```bash
pip install PQLseqPy
```

---

## ✨ Quick Start

```python
import numpy as np
from PQLseqPy import GLMM

# Simulated data
n = 100
rng = np.random.default_rng(0)
X = np.hstack((np.ones((n, 1)), rng.standard_normal((n, 2))))  # intercept + 2 covariates
Y = np.hstack((rng.integers(0, 10, (n, 1)), rng.integers(1, 10, (n, 1))))  # successes, failures
G = rng.standard_normal((n, 500))  # genotypes or random features
K = G @ G.T  # covariance (n x n, PSD)

# Fit the model (default: infer tau1, tau2)
res = GLMM(X, Y, K).fit()

# Summaries
param, coef = res.summary()
print(param)  # model/variance params
print(coef)   # fixed effects table
```

---

## 🧠 Conceptual Model

* **Response**: `Y` is an `(n, 2)` matrix: first column = successes, second column = failures (binomial).
* **Mean**: `μ_i = lib_size_i * logistic(η_i)`, where `lib_size_i = Y[i,0] + Y[i,1]`
* **Linear predictor**: `η = Xβ + u`, with random effect `u ~ N(0, τ1 K + τ2 I)`
* **Variance components**:

  * `τ1` scales the provided covariance `K`
  * `τ2` is the iid (residual) Gaussian component

---

## 🧾 API Reference

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
| `K`                     | `np.ndarray`                    | Covariance (PSD) for the random effect.                                                                                                   |
| `fixed_tau`             | `tuple(float, float)` or `None` | If provided `(tau1, tau2)`, variance components are **fixed** during fitting.                                                             |
| `tau2_set_to_zero`      | `bool`                          | If `True` and `fixed_tau is None`, constrain `τ2 = 0` during inference.                                                                   |
| `verbose`               | `bool`                          | Print per-iteration diagnostics.                                                                                                          |
| `starting_step_size`    | `float`                         | Initial step size for Newton/AI updates. Auto-decays every 10 iters.                                                                      |
| `error_tolerance`       | `float`                         | Convergence tolerance (relative change in `β` and `τ`).                                                                                   |
| `max_iter`              | `int`                           | Maximum iterations for the AI loop.                                                                                                       |
| `regularization_factor` | `float`                         | Tikhonov-like regularization in matrix inversions: uses `(A + λI)^{-1}` with `λ = regularization_factor`. Helps ill-conditioned problems. |

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

## 📚 End-to-End Examples (All Arguments)

Below are small, focused examples showing how to use each set of options.

### 1) Default inference of `τ1` and `τ2`

```python
import numpy as np
from PQLseqPy import GLMM

n = 200
rng = np.random.default_rng(1)
X = np.hstack((np.ones((n, 1)), rng.standard_normal((n, 3))))
Y = np.hstack((rng.integers(0, 20, (n, 1)), rng.integers(1, 20, (n, 1))))
K = rng.standard_normal((n, 300)); K = K @ K.T

res = GLMM(X, Y, K).fit()
param, coef = res.summary()
print(param)
print(coef)
```

### 2) **Null model** (intercept-only fixed effects)

```python
import numpy as np
from PQLseqPy import GLMM

n = 150
rng = np.random.default_rng(2)
X = np.ones((n, 1))  # only intercept
Y = np.hstack((rng.integers(0, 15, (n, 1)), rng.integers(1, 15, (n, 1))))
K = rng.standard_normal((n, 200)); K = K @ K.T

res = GLMM(X, Y, K).fit()
param, coef = res.summary()
print(param)  # β is scalar (intercept only)
print(coef)
```

### 3) **Fixed variance components** (`fixed_tau=(τ1, τ2)`)

```python
import numpy as np
from PQLseqPy import GLMM

n = 120
rng = np.random.default_rng(3)
X = np.hstack((np.ones((n,1)), rng.standard_normal((n,2))))
Y = np.hstack((rng.integers(0, 12, (n, 1)), rng.integers(1, 12, (n, 1))))
K = rng.standard_normal((n, 100)); K = K @ K.T

# Force tau1=0.05, tau2=0.01
res = GLMM(X, Y, K, fixed_tau=(0.05, 0.01)).fit()
param, coef = res.summary()
print(param[['tau1','tau2','variance_model']])
print(coef)
```

### 4) **Constrain τ2 = 0** during inference

```python
import numpy as np
from PQLseqPy import GLMM

n = 180
rng = np.random.default_rng(4)
X = np.hstack((np.ones((n,1)), rng.standard_normal((n,4))))
Y = np.hstack((rng.integers(0, 18, (n, 1)), rng.integers(1, 18, (n, 1))))
K = rng.standard_normal((n, 400)); K = K @ K.T

res = GLMM(X, Y, K, tau2_set_to_zero=True).fit()
param, coef = res.summary()
print(param[['tau1','tau2','variance_model','h2']])
```

### 5) **Verbose** diagnostics + **adaptive step size** showcase

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
    verbose=True,             # print per-iteration info
    starting_step_size=1.0,   # initial step size
    error_tolerance=1e-6,     # tighter convergence
    max_iter=300              # allow more iterations if needed
).fit()
print(res.summary()[0])
```

### 6) **Regularization** for ill-conditioned matrices

```python
import numpy as np
from PQLseqPy import GLMM

n = 120
rng = np.random.default_rng(7)

# Collinear X to force ill-conditioning
x1 = rng.standard_normal(n)
X = np.column_stack([np.ones(n), x1, x1 + 1e-6 * rng.standard_normal(n)])

Y = np.hstack((rng.integers(0, 12, (n, 1)), rng.integers(1, 12, (n, 1))))
K = rng.standard_normal((n, 90)); K = K @ K.T

res = GLMM(X, Y, K, regularization_factor=1e-3).fit()
param, coef = res.summary()
print(param)
print(coef)
```

---

## 📖 Citation

If you use this software in academic work, please cite:

* **Akbari, et al. (2024)**. *Pervasive findings of directional selection realize the promise of ancient DNA to elucidate human adaptation.* (PMID: 39314480)

---