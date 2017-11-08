# Convert a Wright-Fisher trajectory stored as a numpy file into a format readable by wfpath

import sys
import argparse
import numpy as np                          # numerical tools
from timeit import default_timer as timer   # timer for performance


###### Main functions ######


def usage():
    print("")


def main(verbose=False):
    """ Convert a numpy trajectory into plain text save the results. """
    
    # Read in simulation parameters from command line
    
    parser = argparse.ArgumentParser(description='Convert Wright-Fisher from numpy to plain text.')
    parser.add_argument('-i',   type=str, default='data/trajectory', help='input/output file (without extension)')
    parser.add_argument('-t',   type=int, default=0,                 help='start generation')
    parser.add_argument('-T',   type=int, default=0,                 help='final generation (after start)')
    parser.add_argument('--ns', type=int, default=0,                 help='number of sequences to sample at each time point')
    parser.add_argument('--dt', type=int, default=1,                 help='spacing between generations')
    parser.add_argument('-s',   type=int, default=None,              help='seed for random number generator')
    
    arg_list = parser.parse_args(sys.argv[1:])
    
    io_str = arg_list.i
    t0     = arg_list.t
    T      = arg_list.T
    ns     = arg_list.ns
    dt     = arg_list.dt
    seed   = arg_list.s
    
    rng = np.random.RandomState()
    if seed!=None:
        rng = np.random.RandomState(seed)
    
    # Import data and (optionally) downsample
    
    t    = np.load(io_str+'.npz', encoding = 'latin1')
    nVec = t['nVec']
    sVec = t['sVec']
    
    if T==0: T = len(nVec)
    
    f = open(io_str+'_T%d_ns%d_dt%d.dat' % (T, ns, dt), 'w')
    for i in range(t0, T+1, dt):
        if ns==0 or ns>=np.sum(nVec[i]):
            for j in range(len(nVec[i])):
                f.write('%d\t%d\t%s\n' % (i, nVec[i][j], ' '.join([str(int(k)) for k in sVec[i][j]])))
        else:
            iVec = np.zeros(np.sum(nVec[i]))
            ct   = 0
            for j in range(len(nVec[i])):
                iVec[ct:ct+nVec[i][j]] = j
                ct += nVec[i][j]
            iSample = rng.choice(iVec, ns, replace=False)
            for j in range(len(nVec[i])):
                nSample = np.sum(iSample==j)
                if nSample>0:
                    f.write('%d\t%d\t%s\n' % (i, nSample, ' '.join([str(int(k)) for k in sVec[i][j]])))
    f.close()


if __name__ == '__main__': main()

