import numpy as np
import pandas as pd
import time

import statsmodels.api as sm
import statsmodels.formula.api as smf
import scipy.stats

class GLMM:
    """
    GLMM: Generalized Linear Mixed Model for count data
    -------------------------------------------------------------
    This class provides an implementation of a Generalized Linear Mixed Model (GLMM) 
    with a logit link function and binomial family distribution. It is based on the 
    Binomial mixed model from PQLseq (Sun et al. 2019; PMID: 30020412), with added 
    flexibility and significant performance improvements. Key features include:
    
    - Ability to calculate the null model (only intercept as a covariate).
    - Support for predefined variance components (fixed_tau).
    - Order-of-magnitude faster computation compared to the original PQLseq implementation.
    - Improved numerical stability.
    
    Parameters:
    -----------
    X : numpy.ndarray
        Covariate matrix of size (n, k), where n is the sample size and k is the 
        number of covariates (including the intercept). The first column must be 
        all ones to represent the intercept.
    
    Y : numpy.ndarray
        A matrix of size (n, 2), where the first column represents the count of 
        successes (allele 1), and the second column represents the count of failures 
        (allele 0).
    
    K : numpy.ndarray
        An (n, n) covariance matrix structure for the random effect.
    
    fixed_tau : tuple, optional (default=None)
        A tuple (tau1, tau2) specifying the variance components:
        - tau1: Variance component for K.
        - tau2: Variance component for the identity matrix I (Gaussian error).
        If None, the model infers tau1 and tau2.
    
    tau2_set_to_zero : bool, optional (default=False)
        If True, forces tau2 = 0 during inference (only applicable if fixed_tau is None).
    
    verbose : bool, optional (default=False)
        Enables detailed logging for debugging purposes.
    
    starting_step_size : float, optional (default=1)
        Initial step size for the Newton-Raphson updates in the AI algorithm.
    
    error_tolerance : float, optional (default=1e-5)
        Convergence threshold for parameter updates.
    
    max_iter : int, optional (default=200)
        Maximum number of iterations before the algorithm terminates.
    
    regularization_factor : float, optional (default=0)
        Regularization parameter for the matrix inversion, approximating the inverse 
        of A by (A + I * regularization_factor).
    
    Methods:
    --------
    - matrix_inv: Computes the regularized inverse of a matrix.
    - update_step_size: Dynamically adjusts the step size during iterations.
    - check_convergence: Checks convergence based on changes in parameters.
    - initialize_tau: Initializes the tau vector based on user input or default settings.
    - summary: Returns a summary of the fitted model, including parameters and statistics.
    - fit: Fits the GLMM using the AI algorithm and returns the model object.


    USAGE:
        import numpy as np
        n = 100
        np.random.seed(0)
        X = np.hstack((np.ones((n, 1)), np.random.randn(n, 2)))
        Y = np.hstack((np.random.randint(0, 10, (n,1)), np.random.randint(1, 10, (n, 1))))
        G = np.random.randn(n, 10000)
        K = G @ G.T
        res = GLMM(X, Y, K).fit()
        param, coef = res.summary()
        print('-'*50)
        print(param)
        print('-'*50)
        print(coef)
        print('-'*50)
        --------------------------------------------------
        converged                   True
        variance_model    tau1>0, tau2=0
        iter                          35
        elapsed_time            0.040148
        tau1                    0.000033
        tau2                         0.0
        sigma2                  0.000033
        h2                           1.0
        dtype: object
        --------------------------------------------------
        GLMM      beta   se_beta    z_beta    p_beta
        x1   -0.117111  0.089267 -1.311925  0.189545
        x2    0.078483  0.089778  0.874188  0.382016
        x3   -0.141575  0.089530 -1.581311  0.113807
        --------------------------------------------------
    """
    def __init__(self, X, Y, K, fixed_tau=None, tau2_set_to_zero=False, verbose=False,
                 starting_step_size=1, error_tolerance=1e-5, max_iter=200, regularization_factor=0):
        self.X = X
        self.Y = Y
        self.K = K
        self.fixed_tau = np.asarray(fixed_tau) if fixed_tau is not None else None
        self.tau2_set_to_zero = tau2_set_to_zero
        self.verbose = verbose
        self.step_size = starting_step_size
        self.error_tolerance = error_tolerance
        self.max_iter = max_iter
        self.regularization_factor = regularization_factor        
        
    def matrix_inv(self, X, regularization_factor=0):
        I = np.eye(X.shape[0])
        return np.linalg.solve(X + I * regularization_factor, I)

    def update_step_size(self, iter, step_size):
        if (iter + 1) % 10 == 0:
            step_size *= 0.9
        return step_size

    def check_convergence(self, beta, beta0, tau, tau0):
        diff1 = max(abs(beta - beta0) / (abs(beta) + abs(beta0) + self.error_tolerance))
        diff2 = max(abs(tau - tau0) / (abs(tau) + abs(tau0) + self.error_tolerance))
        return (2*max(diff1, diff2)) < self.error_tolerance

    def initialize_tau(self):
        if self.fixed_tau is not None:
            return self.fixed_tau
        elif self.tau2_set_to_zero:
            return np.asarray([1, 0])
        else:
            return np.asarray([1, 1])
            
    def summary(self):
        par = pd.Series(None)
        par['converged'] = self.converged
        par['variance_model'] = self.variance_model
        par['iter'] = self.iter
        par['elapsed_time'] = self.elapsed_time
        par['tau1'] = self.tau[0]
        par['tau2'] = self.tau[1]
        par['sigma2'] = self.sigma2
        par['h2'] = self.h2
        
        estimates = pd.DataFrame([self.beta, self.se_beta, self.z_beta, self.p_beta], index=['beta', 'se_beta', 'z_beta', 'p_beta']).T
        estimates.index = ['x%i'%i for i in range(1, estimates.shape[0]+1)]
        estimates.columns.name = 'GLMM'
        return par, estimates
    
    def fit(self):
        starting_time = time.time()
        self.glm = sm.GLM(self.Y, self.X, family=sm.families.Binomial()).fit()
        X = self.X
        K = self.K
        y = self.Y[:, 0]
        lib_size = self.Y.sum(1)
        beta = self.glm.params
        tau = self.initialize_tau()
        eta = X @ beta
        mu = None
        g_prime_mu = None
        Ytilde = None

        l = []
        for iter in range(self.max_iter):
            if self.verbose:
                print('--------------- iter:', iter+1, '---------------')
            l+=[eta]
        
            # update mu, g_prime_mu, Ytilde
            mu = 1 / (1 + np.exp(-eta))*lib_size
            g_prime_mu = 1/(mu * (lib_size - mu) / lib_size)
            Ytilde = eta + g_prime_mu * (y-mu)
        
            # update H
            D = np.diag(1/g_prime_mu)
            Dinv = np.diag(g_prime_mu)
            I = np.eye(K.shape[0])
            V = tau[0]*K + tau[1]*I
            H = Dinv+V
            Hinv = self.matrix_inv(H, self.regularization_factor)
        
            # update XtHinvX_inv
            HinvX = Hinv @ X
            XtHinvX = X.T @ HinvX
            XtHinvX_inv = self.matrix_inv(XtHinvX, self.regularization_factor)
        
            # update P
            P = Hinv - HinvX @ XtHinvX_inv @ HinvX.T
            P = (P+P.T)/2
        
        
            # update AI    
            PK = P @ K
            PI = P
            PYtilde = P @ Ytilde
            IPYtilde = PYtilde
            KPYtilde = K @ PYtilde
            AI00 = KPYtilde.T @ P @ KPYtilde
            AI01 = KPYtilde.T @ P @ IPYtilde
            # AI10 = IPYtilde.T @ P @ KPYtilde
            AI10 = AI01
            AI11 = IPYtilde.T @ P @ IPYtilde
            AI = np.asarray([[AI00, AI01], [AI10, AI11]])
            d_qlr_tau1 = (PYtilde.T @ KPYtilde - np.trace(PK)) # Drop /2 To be consistent with PQLseq2
            d_qlr_tau2 = (PYtilde.T @ IPYtilde - np.trace(PI)) # Drop /2 To be consistent with PQLseq2
            d_qlr_tau = np.asarray([d_qlr_tau1, d_qlr_tau2])
            if self.verbose:
                print('d_qlr_tau: ', d_qlr_tau)
            
            # update tau
            if self.fixed_tau is not None:
                tau0 = tau
                variance_model = 'Fixed tau'
            else:
                tau0 = tau
                if all(tau0>0):
                    variance_model = 'tau1>0, tau2>0'
                    AI_inv = np.linalg.solve(AI, d_qlr_tau)
                elif all(tau0==0):
                    variance_model = 'tau1=0, tau2=0 (GLM)'
                    AI_inv = np.asarray([0,0])
                elif tau0[0]==0:
                    variance_model = 'tau1=0, tau2>0 (restart with tau1=1 and tau2=0)'
                    AI_inv = np.asarray([0,0])
                    tau0 = np.asarray([1,0])
                elif tau0[1]==0:
                    variance_model = 'tau1>0, tau2=0'
                    AI_inv = np.append(d_qlr_tau[0]/AI[0,0], 0)


                    
                if self.verbose:
                    print('variance_model: ', variance_model)
                step_size_iter = self.step_size
                while True:
                    tau = tau0 + step_size_iter*AI_inv
                    if np.isnan(tau[0]) or np.isnan(tau[1]):
                        raise ValueError(f"NaN in tau at iteration {iter+1}: tau = {tau}.")

                    if self.verbose:
                        print('step_size_iter:', step_size_iter)
                        print('tau0:', tau0)
                        print('tau:', tau)
                    tau[(tau0 < self.error_tolerance) & (tau < self.error_tolerance)] = 0
                    if any(tau<0):
                        step_size_iter = step_size_iter/2
                    else:
                        break
        
        
            # update beta
            beta0 = beta
            beta = XtHinvX_inv @ HinvX.T @ Ytilde
        
            # update eta
            eta = Ytilde - Dinv @ Hinv @ (Ytilde-X @ beta)  
        
            # update step_size
            self.step_size = self.update_step_size(iter, self.step_size)
        
            # Check Convergence 
            self.converged = self.check_convergence(beta, beta0, tau, tau0)
            if self.converged:
                break
                
        se_beta = np.diag(XtHinvX_inv)**0.5
        z_beta = beta/se_beta
        p_beta = scipy.stats.chi2.sf(z_beta**2, 1)
        sigma2 = sum(tau)
        if sigma2==0:
            h2=np.nan
        else:
            h2 = tau[0]/sigma2
        
        self.beta = beta
        self.se_beta = se_beta
        self.z_beta = z_beta
        self.p_beta = p_beta
        self.tau = tau
        self.sigma2 = sigma2
        self.h2 = h2
        self.variance_model = variance_model
        self.cov = XtHinvX_inv
        self.iter = iter
        self.elapsed_time = time.time()-starting_time
        return self