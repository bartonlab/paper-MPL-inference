import math
import numpy as np
import argparse
import sys
np.set_printoptions(linewidth=140, precision=5, suppress=True)

from likelihood import GaussianLikelihood as Likelihood

# This file is modified from Terhorst et al.'s example.py

def usage():
    print("")

def main(verbose=False):
    """
    Load a saved trajectory, estimate selection coefficients,
    and record the results.
    """

    # Read in data from command line
    
    parser = argparse.ArgumentParser(description='Estimate selection coefficients from an input trajectory.')
    parser.add_argument('-i', type=str,               help='input file')
    parser.add_argument('-o', type=str,               help='output file')
    parser.add_argument('-N', type=int, default=1000, help='effective population size')
    parser.add_argument('-l', action='store_true',    help='enable estimates with linkage')
    
    arg_list = parser.parse_args(sys.argv[1:])
    
    in_str      = arg_list.i
    out_str     = arg_list.o
    N           = arg_list.N
    use_linkage = arg_list.l

    # Load data in Terhorst et al. format

    _data = np.loadtxt(in_str)
    _L    = len(_data[0][2:])
    t0    = _data[0][0]
    t1    = 0
    for i in range(len(_data)):
        if _data[i][0]!=t0:
            t1 = i
            break

    times     = np.unique(_data.T[0])[1:]   # drop first time point (t=0)
    positions = np.array(range(1, _L+1),int)
    init_haps = []
    for i in range(t1):
        for j in range(int(_data[i][1])):
            init_haps.append(_data[i][2:])
    init_haps = np.array(init_haps)

    r         = 1e-8    # true recombination rate is zero, but zero is not allowed
    data      = np.zeros((len(times),_L,1))
    for i in range(len(times)):
        _t_data = np.array([_d[2:] for _d in _data if _d[0]==times[i]])
        _t_num  = np.array([ _d[1] for _d in _data if _d[0]==times[i]])
        _t_sum  = np.einsum('i,ij->j', _t_num, _t_data) / float(np.sum(_t_num))
        for j in range(_L):
            data[i][j][0] = _t_sum[j]
            # bump frequencies away from boundaries to prevent pathological inference
            if data[i][j][0]==1: data[i][j][0] -= 1/float(N)
            if data[i][j][0]==0: data[i][j][0] += 1/float(N)

    # Compute selection coefficients and save output

    _s      = []
    all_idx = list(range(len(positions)))
    for i, pos in enumerate(positions):
        dta = data[:, (i,), :]
        #if np.min(dta[-1]) < 0.02:
        #    # Skip sites which are are lost
        #    continue
        ih    = init_haps[:, (i,)]
        lik   = 0
        if use_linkage: lik = Likelihood(data[:, all_idx, :], init_haps[:, all_idx], times, positions[all_idx], pos, None, verbose=False, reset=True)
        else:           lik = Likelihood(data[:,    (i,), :], init_haps[:,    (i,)], times,              [pos], pos, None, verbose=False, reset=True)
        fit   = lik.maxlik(N=N, log10r=math.log10(r), h=0.5, bounds=(-.25, .25))
        s_hat = fit.x
        _s.append(s_hat)

    # Save output

    f = open(out_str, 'w')
    for s in _s: f.write('%lf\n' % s)
    f.close()

if __name__ == '__main__': main()

