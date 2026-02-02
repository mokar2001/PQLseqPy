# PQLseqPy
[![GitHub release](https://img.shields.io/github/v/release/mokar2001/PQLseqPy)](ChangeLog.md)
[![PyPI version](https://img.shields.io/pypi/v/PQLseqPy.svg)](https://pypi.org/project/PQLseqPy/)
[![PyPI total downloads](https://static.pepy.tech/badge/pqlseqpy)](https://pepy.tech/project/pqlseqpy)
[![Bioconda version](https://anaconda.org/bioconda/PQLseqPy/badges/version.svg)](https://anaconda.org/bioconda/PQLseqPy)
[![Bioconda downloads](https://anaconda.org/bioconda/PQLseqPy/badges/downloads.svg)](https://anaconda.org/bioconda/PQLseqPy)
[![Bioconda release date](https://anaconda.org/bioconda/PQLseqPy/badges/latest_release_date.svg)](https://anaconda.org/bioconda/PQLseqPy)

**PQLseqPy** is a fast Python implementation of Penalized Quasi-Likelihood for **binomial count data**, designed for **selection analysis using allele counts** (Akbari et al. 2026).  
It is inspired by **PQLseq** (Sun et al. 2019), and adapted to population-genetic settings where observations are **alternative vs reference allele counts**, possibly **pooled across multiple individuals**.

The method implements a Generalized Linear Mixed Model (GLMM) with a **logit link**, supports flexible random-effect covariance structures, and provides substantial performance and numerical-stability improvements over the original PQLseq implementation.

---

## 📦 Installation

Python ≥ 3.8 recommended.

### Conda
```bash
conda install -c bioconda PQLseqPy
````

### Pip

```bash
pip install PQLseqPy
```

---

## 🚀 Quick start

```python
from PQLseqPy import GLMM

res = GLMM(X, Y, K).fit()
param, coef = res.summary()
```

This fits the model and returns fixed-effect estimates and variance components.
See the sections below for a detailed description of the model, assumptions, and parameters.

---

## 🧾 Model Overview

Each **row `i` corresponds to a single observational unit** in the analysis.
An observational unit represents a set of allele counts that share the same covariate values and random-effect structure.

In practice, a row may correspond to:

* a single individual at a locus, or
* allele counts pooled across multiple individuals at the same locus.

Each observation is modeled as a **binomial count**.

For each row `i`:

* `Y[i, 0]` = number of **alternative alleles**
* `Y[i, 1]` = number of **reference alleles**
* `n_i = Y[i, 0] + Y[i, 1]` is the **total allele count**

Counts may arise from **diploid or pseudo-haploid data**.

---

### Fixed effects and random effects

The model relates allele counts to covariates through:

* **Fixed effects (`β`)**: coefficients shared across all observations (for example, intercept, time, or other covariates).
* **Random effects (`u`)**: an unobserved additive term that captures correlated residual structure across observations not explained by the fixed effects.

Random effects are treated as a random vector, not as additional covariates.

The observation model is:

```
Y[i, 0] ~ Binomial(n_i, p_i)
logit(p_i) = X[i] β + u_i
```

The allele frequency `p_i` is **not observed directly**. It is an implicit quantity defined by the linear predictor `X[i] β + u_i` through the logistic function.

---

### Covariance structure of the random effects

Random effects are modeled as a **latent multivariate normal vector**:

```
u ~ N(0, V)
```

Here, `V` is the **covariance matrix** of the random effects, describing how residual deviations co-vary across observations.

In this implementation:

```
V = tau1 * K + tau2 * I
```

where:

* `K` is a user-supplied structure matrix (for example, a genetic relatedness matrix),
* `I` is the identity matrix,
* `tau1` and `tau2` are variance components.

The random effects `u` are **latent** and integrated out during model fitting.
The variance components `tau1` and `tau2` are estimated explicitly.

---

## 🧾 Parameters & Attributes

```python
GLMM(
  X, Y, K,
  fixed_tau=None,
  tau2_set_to_zero=False,
  verbose=False,
  starting_step_size=1,
  error_tolerance=1e-5,
  max_iter=200,
  regularization_factor=0
)
```

### Parameters

| Name                    | Type                       | Description                                                      |
| ----------------------- | -------------------------- | ---------------------------------------------------------------- |
| `X`                     | `np.ndarray (n, k)`        | Covariate matrix. **First column must be ones** (intercept).     |
| `Y`                     | `np.ndarray (n, 2)`        | Allele counts: `[alternative, reference]`.                       |
| `K`                     | `np.ndarray (n, n)`        | Structure matrix used to construct the random-effect covariance. |
| `fixed_tau`             | `(float, float)` or `None` | If provided, `(tau1, tau2)` are held fixed during fitting.       |
| `tau2_set_to_zero`      | `bool`                     | If `True` and `fixed_tau is None`, constrains `tau2 = 0`.        |
| `verbose`               | `bool`                     | Print per-iteration diagnostics.                                 |
| `starting_step_size`    | `float`                    | Initial step size for AI updates.                                |
| `error_tolerance`       | `float`                    | Relative convergence tolerance for `beta` and `tau`.             |
| `max_iter`              | `int`                      | Maximum number of iterations.                                    |
| `regularization_factor` | `float`                    | Diagonal regularization used in matrix inversions.               |

**Notes on behavior**

* When `fixed_tau` is provided, the covariance structure is fixed but **`beta` is still iteratively updated** until convergence.
* `tau2_set_to_zero=True` enforces a **hard constraint** `tau2 = 0` throughout fitting. The model is initialized with `(tau1 = 1, tau2 = 0)`, estimates `tau1` only, and assumes no independent (identity) variance component.
* Internally, `n_i` is computed as `Y.sum(axis=1)` (named `lib_size` in the code).

---

### Attributes (after `.fit()`)

| Name             | Type                | Description                                            |
| ---------------- | ------------------- | ------------------------------------------------------ |
| `beta`           | `np.ndarray (k,)`   | Fixed-effect estimates.                                |
| `se_beta`        | `np.ndarray (k,)`   | Standard errors of `beta`.                             |
| `z_beta`         | `np.ndarray (k,)`   | z-statistics.                                          |
| `p_beta`         | `np.ndarray (k,)`   | Two-sided p-values (`z² ~ chi-square(1)`).             |
| `tau`            | `np.ndarray (2,)`   | Estimated variance components `(tau1, tau2)`.          |
| `sigma2`         | `float`             | `tau1 + tau2`.                                         |
| `h2`             | `float` or `nan`    | `tau1 / (tau1 + tau2)` (nan if `sigma2 = 0`).          |
| `variance_model` | `str`               | Variance-structure branch used in the final iteration. |
| `cov`            | `np.ndarray (k, k)` | Approximate covariance of `beta`.                      |
| `converged`      | `bool`              | Convergence flag.                                      |
| `iter`           | `int`               | Number of iterations used.                             |
| `elapsed_time`   | `float`             | Runtime in seconds.                                    |

Possible values of `variance_model` include:

* `tau1>0, tau2>0`
* `tau1>0, tau2=0`
* `tau1=0, tau2=0 (GLM)`
* `Fixed tau`

---

### Methods

* `fit() -> GLMM`
  Fits the model in place and returns `self`.

* `summary() -> (pandas.Series, pandas.DataFrame)`
  Returns:

  * `param`: convergence status, variance model, runtime, `tau1`, `tau2`, `sigma2`, `h2`
  * `estimates`: `beta`, `se_beta`, `z_beta`, `p_beta` indexed by `x1`, `x2`, …

---

## 📚 Example

```python
import numpy as np
from PQLseqPy import GLMM

n = 100
rng = np.random.default_rng(5)

X = np.hstack((np.ones((n, 1)), rng.standard_normal((n, 2))))
Y = np.hstack((
    rng.integers(0, 10, (n, 1)),
    rng.integers(1, 10, (n, 1))
))

G = rng.standard_normal((n, 10000))
K = G @ G.T

res = GLMM(X, Y, K).fit()
param, coef = res.summary()

print(param)
print(coef)
```

---

## ⚠️ Notes

* The model operates on **allele counts**, not observed allele frequencies.
* Allele frequencies `p_i` and random effects `u_i` are implicit quantities defined by the model.
* Use `regularization_factor` if numerical instability occurs.

---

## 📖 Citation

If you use this software in academic work, please cite:

* **Akbari et al. (2026)**. *Ancient DNA reveals pervasive directional selection across West Eurasia.*
* **Sun et al. (2019)**. *Heritability estimation and differential analysis of count data with generalized linear mixed models in genomic sequencing studies.* PMID: 30020412
