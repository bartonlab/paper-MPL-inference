from __future__ import division, print_function

# Code relating to the multivariate normal is only present
# in later version of SciPy
import pkg_resources
pkg_resources.require("scipy>=0.14")

import numpy as np
from math import log, pi, exp, isnan, sqrt
import scipy.optimize
import scipy.misc
import scipy.stats
import scipy.integrate
import warnings
import logging
import psutil
import itertools
import collections
import os

# This must be built prior to running this script.
import pyexp

def logfac(x):
    return scipy.special.gammaln(x + 1)

def inv_qf(M, x):
    # if __debug__: print("inv_qf")
    # Inverse quadratic form: x^T(M^-1)x
    # More stable than computing M^-1
    v, W = np.linalg.eigh(M)
    s = sum(np.dot(W[:, i].T, x).item()**2 / v[i] for i in range(len(x)) if v[i] > 0)
    return s

def latent_mvnorm_ll_nquad(x, mu, sigma, sampling):
    K = len(x)
    counts = x * sampling
    mvn = scipy.stats.multivariate_normal(mu, sigma)
    a0 = logfac(sampling) - logfac(counts) - logfac(sampling - counts)
    def f(*args):
        Z = np.array(args)
        return exp((a0 + counts * np.log(Z) + (sampling - counts) * np.log(1. - Z)).sum() + mvn.logpdf(Z))
    ret = scipy.integrate.nquad(f, ((0.0, 1.0),) * K)
    print(ret)
    return ret[0]

def latent_mvnorm_ll_laplace(x, mu, sigma, sampling):
    K = len(x)
    sgn, logdetsigma = np.linalg.slogdet(sigma)
    if sgn != 1.0:
        evs = sorted(np.linalg.eigvals(sigma))
        cn = abs(evs[-1] / evs[0])
        warnings.warn("Sigma is not PSD; condition number: %g" % cn)

    # Importance sampler
    counts = x * sampling
    logfac = lambda x: scipy.special.gammaln(x + 1)
    # expo += mvn.logpdf(Z)
    # expo -= prop.logpdf(Z)
    # ret = scipy.misc.logsumexp(expo) - log(splx.sum())
    # Prior distribution
    mvn = scipy.stats.multivariate_normal(mu, sigma)
    sigmainv = np.linalg.inv(sigma)
    a0 = logfac(sampling) - logfac(counts) - logfac(sampling - counts)

    def f(Z):
        return (a0 + counts * np.log(Z) + (sampling - counts) * np.log(1. - Z)).sum()
    def g(Z): 
        Z0 = np.minimum(Z, 0)
        Z1 = np.maximum(Z - 1, 0)
        return -(f(Z) + mvn.logpdf(Z)) + 5 * 0.5 * (np.dot(Z0, Z0) + np.dot(Z1, Z1))
    def grad(Z):
        Z0 = np.minimum(Z, 0)
        Z1 = np.maximum(Z - 1, 0)
        return -((counts / Z) - (sampling - counts) / (1. - Z) - np.dot(sigmainv, Z - mu)) + 5 * (Z0 + Z1)
                 
    xopt, fopt, gopt, Bopt, _, _, _ = scipy.optimize.fmin_bfgs(g, x, grad, full_output=True, gtol=1e-5)
    print(xopt)
    print(mu)
    print(x)
    print(fopt, f(xopt), mvn.logpdf(xopt))
    return -fopt + (K / 2) * log(2 * pi) - (1 / 2) * np.linalg.slogdet(Bopt)[1]

# latent_mvnorm_ll = latent_mvnorm_ll_is

def mvnorm_ll(x, mu, sigma, mask=False):
    'Log-likelihood function of N(mu, sigma) rv.'
    # if __debug__: print("mvnorm_ll")
    # const * -(1/2) * (x - mu)^T (sigma^-1) (x - mu)
    if mask:
        bdry = np.logical_and(x >= 0.05, x <= 0.95)
        if not np.any(bdry):
            print("no data points remain")
            print(x)
            return -np.inf
        x = x[bdry]
        mu = mu[bdry]
        sigma = sigma[np.ix_(bdry, bdry)]
    alpha = 0.05
    # print("pre-shrunk", np.linalg.eigvals(sigma)[[1, -1]])
    sigma = (1. - alpha) * sigma + alpha * np.eye(sigma.shape[0])
    try:
        np.linalg.cholesky(sigma)
    except np.linalg.LinAlgError:
        return -np.inf
    # print("shrunk", np.linalg.eigvals(sigma)[[1, -1]])
    mvn = scipy.stats.multivariate_normal(mu, sigma)
    return mvn.logpdf(x)

class Likelihood:
    Z = [ 
            (True, False, False), # zL
            (False, True, False), # zS
            (False, False, True), # zR
            (True, True, False), # zLS
            (True, False, True), # zLR
            (False, True, True), # zRS
            (True, True, True) # zLRS
        ]

    def __init__(self, data, init_haps, times, positions, selected_position, sampling,
            verbose=False, reset=False):
        '''Initialize the model. 

            T is the number of sampled generations; L is the length of
            the region being considered; K is number of founder lines;
            R is number of replicates.
            
            data: T x L x R matrix of derived allele frequencies
            init_haps: K x L matrix of initial haplotypes (boolean)
            times: T x 1 vector of time periods in which sampling occurred
            positions: L x 1 vector of positions at which sampling occurred
            selected_position: \in positions; position which is under selection. 
                               All other positions assumed neutral.
            sampling: T x L x R matrix if coverage depth (or something broadcastable to
                      those dimensions.) Uses latent model with binomial sampling overlaid. (Slow.)
                      If None, then use the non-latent model.
            verbose: Emit extra debugging info.
            reset: Clear cached recursive function call evaluations
                   during each call to mean_cov(). Needed in some cases to
                   control memory consumption.
        '''
        self.verbose = verbose
        self.reset = reset
        self.selected_position = selected_position
        self.data = data
        self.sampling = sampling
        if sampling is not None:
            sampling = np.array(sampling)
            if sampling.size == 1:
                self.sampling = np.zeros(data.shape)
                self.sampling[:] = sampling
        self.times = np.array(times, dtype=np.int32)
        self.positions = np.array(positions, dtype=np.int32)
        self.dists = abs(np.subtract.outer(self.positions, self.positions))
        self.n_replicates = self.data.shape[-1]
        T = len(self.times)
        L = len(self.positions)
        self.M = np.zeros([self.n_replicates, T, L], dtype=np.double)
        self.C = np.zeros([self.n_replicates, T, L, T, L], dtype=np.double)

        # Construct initial haplotypes
        if np.array(init_haps).shape == (L, L, 7):
            self.init_haps = init_haps
        else:
            assert init_haps.shape[1] == len(self.positions)
            self.init_haps = np.zeros([L, L, 7], dtype=np.double)
            F = init_haps.shape[0]
            for (i1, l1), (i2, l2) in itertools.product(enumerate(self.positions), repeat=2):
                if selected_position:
                    sp = self.positions.tolist().index(selected_position)
                else:
                    sp = i1
                cc = collections.Counter(map(tuple, init_haps[:, (i1, sp, i2)]))
                self.init_haps[i1, i2] = [cc[k] / F for k in self.Z]

    def mean_cov(self, N, r, h, s):
        '''Update the mean and covariance for this model for given
        parameter values.

            N: Effective population size.
            r: Recombination rate.
            h: Dominance paremeter.
            s: Fitness parameter.
        '''
        # DEBUG: CHANGED psutil.phymem_usage().percent TO psutil.virtual_memory().percent
        # FOR COMPATIBILITY WITH CURRENT VERSION
        if self.reset or psutil.virtual_memory().percent > 90.0:
            pyexp.reset()
        # Update mean and covariance matrices in place. Returns None.
        for i in range(self.n_replicates):
            pyexp.mean_cov(self.M[i], self.C[i], self.positions, 
                    list(self.positions).index(self.selected_position), 
                    self.init_haps, self.times, N, r, h, s,
                    False)

    def print_data(self):
        for i in range(self.data.shape[-1]):
            print(self.data[..., i])

    def maxlik(self, **kwargs):
        '''Maximize likelihood over selection parameter.'''
        if self.verbose: 
            self.print_data()
        def negll(x):
            if self.verbose:
                print(x)
            N = kwargs.get('N', x)
            log10r = kwargs.get('log10r', x)
            h = kwargs.get('h', x)
            return -self.likelihood(N, 10**log10r, h, x)
        return scipy.optimize.minimize_scalar(negll, bounds=kwargs['bounds'], 
                method='bounded', options={'disp': True})

class GaussianLikelihood(Likelihood):
    def _latent_mvnorm_ll(self, x, mu, sigma, sampling):
        K = len(x)
        try:
            np.linalg.cholesky(sigma)
        except np.linalg.LinAlgError:
            print(sigma)
            return -np.inf
        # Importance sampler
        counts = x * sampling
        a0 = logfac(sampling) - logfac(counts) - logfac(sampling - counts)
        # Prior distribution
        mvn = scipy.stats.multivariate_normal(mu, sigma)
        NS = 20000
        Z = mvn.rvs(NS)
        if len(Z.shape) == 1:
            Z.shape = (Z.shape[0], 1)
        splx = np.all([Z > 0, Z < 1], axis=(0, 2))
        Z = Z[splx]
        s = (a0 + counts * np.log(Z) + (sampling - counts) * np.log(1. - Z)).sum(axis=1)
        val = scipy.misc.logsumexp(s) - log(splx.sum())
        return val

    def likelihood(self, N, r, h, s):
        self.mean_cov(N, r, h, s)
        ipos = list(range(len(self.positions)))
        itimes = list(range(len(self.times)))
        dinds = np.ix_(itimes, ipos)
        ddinds = np.ix_(itimes, ipos, itimes, ipos)
        M = self.M[0][dinds]
        C = self.C[0][ddinds]
        dd = C.shape[0] * C.shape[1]
        dtas = self.data[dinds]
        if self.sampling is not None:
            lmv = [self._latent_mvnorm_ll(dtas[..., i].reshape(dd), 
                   M.reshape(dd), C.reshape(dd, dd), 
                   self.sampling[..., i].reshape(dd))
                   for i in range(self.n_replicates)]
            return sum(lmv)
        else:
            mv = [mvnorm_ll(dtas[..., i].reshape(dd), M.reshape(dd), C.reshape(dd, dd))
                  for i in range(self.n_replicates)]
            return sum(mv)
