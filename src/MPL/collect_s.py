# Convert a Wright-Fisher trajectory stored as a numpy file into a format readable by wfpath

import os
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
    parser.add_argument('-i',   type=str, default='data/trajectory',  help='input/output file (without extension)')
    parser.add_argument('-n',   type=int, default=0,                  help='number of trajectories')
    parser.add_argument('-m',   type=str,            action='append', help='inference methods')
    parser.add_argument('-t',   type=int,            action='append', help='starting generation')
    parser.add_argument('-T',   type=int,            action='append', help='final generation')
    parser.add_argument('--ns', type=int,            action='append', help='number of sequences to sample at each time point')
    parser.add_argument('--dt', type=int,            action='append', help='spacing between generations')
    
    def s_str(i, n, m, T, ns, dt):
        return '_'.join([i, str(n), 'T'+str(T), 'ns'+str(ns), 'dt'+str(dt), m])
    
    arg_list = parser.parse_args(sys.argv[1:])
    
    io_str  = arg_list.i
    n_list  = range(arg_list.n)
    m_list  = arg_list.m
    t_list  = arg_list.t
    T_list  = arg_list.T
    ns_list = arg_list.ns
    dt_list = arg_list.dt
    
    # Process input parameters
    
    t_list  = np.sort(np.unique(np.array( t_list,int)))
    T_list  = np.sort(np.unique(np.array( T_list,int)))
    ns_list = np.sort(np.unique(np.array(ns_list,int)))
    dt_list = np.sort(np.unique(np.array(dt_list,int)))
    
    # Get sequence length for header
    
    d = np.loadtxt(s_str(io_str, 0, m_list[0], T_list[0], ns_list[0], dt_list[0])+'.dat')
    L = len(d)
    
    head = 'trajectory,method,t0,T,ns,deltat,runtime,' + (','.join(['s%d' % i for i in range(L)]))
    f    = open(io_str + '_collected.csv', 'w')
    f.write('%s\n' % head)
    
    # Import data and (optionally) downsample
    
    for n in n_list:
        for t in t_list:
            for T in T_list:
                for ns in ns_list:
                    for dt in dt_list:
                        for m in m_list:
                            d = np.loadtxt(s_str(io_str, n, m, T, ns, dt)+'.dat')
                            r = 1
                            if os.path.isfile(s_str(io_str, n, m, T, ns, dt)+'_time.dat'):
                                r = 1e-9 * float(np.loadtxt(s_str(io_str, n, m, T, ns, dt)+'_time.dat'))
                            f.write('%d,%s,%d,%d,%d,%d,%lf,' % (n, m, t, T, ns, dt, r))
                            f.write(','.join(['%lf' % i for i in d]))
                            f.write('\n')
    
    f.close()


if __name__ == '__main__': main()

