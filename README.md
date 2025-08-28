# qglmm

**qglmm** is a fast implementation of a Generalized Linear Mixed Model (GLMM) for count/binomial data with variance components.  
It is inspired by **PQLseq** (Sun et al. 2019; PMID: 30020412), with added flexibility and significant performance improvements.

## ✨ Features

- Supports **Binomial family with logit link**
- Order-of-magnitude faster than PQLseq
- Handles **variance components** (`tau1`, `tau2`) with options:
  - Fixed values
  - Inference from data
- Stable Newton-Raphson updates with adaptive step size
- Regularization for numerical stability
- Easy to use API similar to `statsmodels`

## 📦 Installation
You can easily install fastGLMM via Conda:
```bash
conda install -c bioconda fastglmm
```

## 🚀 Usage
```python
import numpy as np
from fastglmm import GLMM

# Simulated data
n = 100
np.random.seed(0)
X = np.hstack((np.ones((n, 1)), np.random.randn(n, 2)))
Y = np.hstack((np.random.randint(0, 10, (n, 1)), np.random.randint(1, 10, (n, 1))))
G = np.random.randn(n, 500)
K = G @ G.T

# Fit model
res = GLMM(X, Y, K).fit()

# Summary
param, coef = res.summary()
print(param)
print(coef)
```

