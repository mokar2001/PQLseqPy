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

---

## 🧬 Real Example

In real analysis, your data is likely stored in separate files (phenotypes, covariates, and a kinship matrix). The following example demonstrates how to:
1. Load data from CSVs.
2. **Align samples** to ensure $X$, $Y$, and $K$ match.
3. Add the required intercept term.
4. Fit the model and save results.

### 1. Generate "Dummy" Files (for reproducibility)
*Run this once to create sample files in your working directory.*

```python
import pandas as pd
import numpy as np

# Create 50 samples with IDs 'sample_0' to 'sample_49'
ids = [f"sample_{i}" for i in range(50)]

# 1. Count Data (Y): Successes and Failures
df_counts = pd.DataFrame({
    'methylated': np.random.randint(0, 20, 50),
    'unmethylated': np.random.randint(5, 30, 50)
}, index=ids)
df_counts.to_csv("example_counts.csv")

# 2. Covariates (X): Age, Sex, Batch
df_cov = pd.DataFrame({
    'age': np.random.normal(30, 10, 50),
    'sex': np.random.randint(0, 2, 50), # 0=Male, 1=Female
    'batch': np.random.randint(1, 4, 50)
}, index=ids)
df_cov.to_csv("example_covariates.csv")

# 3. Kinship Matrix (K)
# Generate a random positive semi-definite matrix
G = np.random.randn(50, 100)
K = G @ G.T
df_kinship = pd.DataFrame(K, index=ids, columns=ids)
df_kinship.to_csv("example_kinship.csv")

print("Dummy files generated!")
```

### 2. Load and Analze
```python
import pandas as pd
import numpy as np
from PQLseqPy import GLMM

# Load Data
counts = pd.read_csv("example_counts.csv", index_col=0)
covs = pd.read_csv("example_covariates.csv", index_col=0)
kinship = pd.read_csv("example_kinship.csv", index_col=0)

# Ensure we only use samples present in ALL three files and in the same order
common_ids = counts.index.intersection(covs.index).intersection(kinship.index)

print(f"Analyzing {len(common_ids)} overlapping samples.")

Y_matched = counts.loc[common_ids]
X_matched = covs.loc[common_ids]
K_matched = kinship.loc[common_ids, common_ids]

# Y: Must be [Success, Failure]
Y_matrix = Y_matched[['methylated', 'unmethylated']].values

# X: Add Intercept column (First column must be 1s)
# We use the covariates: Age, Sex, Batch
X_matrix = X_matched.values
X_matrix = np.hstack((np.ones((X_matrix.shape[0], 1)), X_matrix))

# K: Kinship matrix
K_matrix = K_matched.values

# Fit Model
model = GLMM(
    X=X_matrix, 
    Y=Y_matrix, 
    K=K_matrix,
    verbose=True
).fit()

# View & Save Results
stats, estimates = model.summary()

# Map generic x1, x2... back to real names
# Note: x1 is Intercept
feature_names = ['Intercept'] + list(X_matched.columns)
estimates.index = feature_names

print(estimates)
estimates.to_csv("pqlseqpy_results.csv")
```


## 📖 Citation

If you use this software in academic work, please cite:

* **Akbari, et al. (2024)**. *Pervasive findings of directional selection realize the promise of ancient DNA to elucidate human adaptation.* (PMID: 39314480)