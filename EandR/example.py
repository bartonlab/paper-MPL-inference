#!/usr/bin/env python2.7
from __future__ import division
import math
import numpy as np
np.set_printoptions(linewidth=140, precision=5, suppress=True)

from likelihood import GaussianLikelihood as Likelihood

# This script gives an example usage of the file likelihood.py (and
# accompanying libraries) to perform inference. The files <data.npz>,
# <positions.txt> and <init_haps.txt> are a simulated E&R experiment
# over a loosely linked (r=1e-4) region of 1000 loci. The census
# population size was N=2000 and 200 founder lines were used. The locus
# at position 568 has positive fitness of s=0.05 and all other loci are
# neutral. Pooled sequencing was simulated at generations 10, 20, 30,
# 40, 50 by sampling binomially from the underlying populations. Three
# Independent replicates were used. A total of 24 sites were segregating.

# 200 x 24 array of initial haplotypes (1=derived; 0=ancestral)
init_haps = np.loadtxt("example_data/init_haps.txt")
# 5 x 24 x 3 array of derived allele frequencies
data = np.load("example_data/data.npz")['arr_0']
# 24 x 1 vector of positions corresponding to seg. sites
positions = np.loadtxt("example_data/positions.txt", dtype=int)
times = [10, 20, 30, 40, 50]
r = 1e-4
N = 2000
# True selected position is 568

liks = []
print("## Pass 1. Testing")
for i, pos in enumerate(positions):
    dta = data[:, (i,), :]
    if np.min(dta[-1]) < 0.02:
        # Skip sites which are are lost
        continue
    ih = init_haps[:, (i,)]
    lik = Likelihood(data[:, (i,), :], init_haps[:, (i,)], times, [pos], pos, None)
    fit = lik.maxlik(N=N, log10r=math.log10(r), h=0.5, bounds=(-.15, .15))
    s_hat = fit.x
    ml = -fit.fun
    ml0 = lik.likelihood(N, r, 0.5, 0.0)
    lr = 2 * (ml - ml0)
    if lr < 0:
        print(dta, )
    print("Site: %d\tLoglik ratio: %g" % (pos, lr))
    liks.append((pos, lr))

print("## Pass 2. Estimation")
maxpos = sorted(liks, key=lambda tup: -tup[1])[0][0]
print("Site with highest LR: %d" % maxpos)
print("Estimating a 5-locus model at this site...")
mpi = tuple(positions).index(maxpos)
inds = [mpi - 5, mpi - 3, mpi, mpi + 3, mpi + 5]
lik = Likelihood(data[:, inds, :], init_haps[:, inds], times, positions[inds], maxpos, None, True)
fit = lik.maxlik(N=N, log10r=math.log10(r), h=0.5, bounds=(-.15, .15))
s_hat = fit.x
print("Estimated selection coefficient: %f" % s_hat)
