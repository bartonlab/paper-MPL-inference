# A simple Wright-Fisher simulation with an additive fitness model

import sys
import argparse
import numpy as np                          # numerical tools
from timeit import default_timer as timer   # timer for performance


###### Main functions ######


def usage():
    print("")


def main(verbose=False):
    """ Simulate Wright-Fisher evolution of a population and save the results. """
    
    # Read in simulation parameters from command line
    
    parser = argparse.ArgumentParser(description='Wright-Fisher evolutionary simulation.')
    parser.add_argument('-o',   type=str,   default='data/trajectory', help='output destination')
    parser.add_argument('-N',   type=int,   default=1000,              help='population size')
    parser.add_argument('-L',   type=int,   default=50,                help='sequence length')
    parser.add_argument('-T',   type=int,   default=1000,              help='number of generations in simulation')
    parser.add_argument('--mu', type=float, default=1.0e-3,            help='mutation rate')
    parser.add_argument('--nB', type=int,   default=10,                help='number of beneficial mutations')
    parser.add_argument('--fB', type=float, default=0.03,              help='fitness effect of beneficial mutations')
    parser.add_argument('--nD', type=int,   default=10,                help='number of deleterious mutations')
    parser.add_argument('--fD', type=float, default=-0.03,             help='fitness effect of deleterious mutations')
    parser.add_argument('--random', type=int, default=0,               help='number of random starting sequences in the population')
    
    arg_list = parser.parse_args(sys.argv[1:])
    
    out_str   = arg_list.o
    N         = arg_list.N
    L         = arg_list.L
    T         = arg_list.T
    mu        = arg_list.mu
    nB        = arg_list.nB
    fB        = arg_list.fB
    nD        = arg_list.nD
    fD        = arg_list.fD
    randomize = arg_list.random
    
    # Set selection coefficients
    
    h = np.zeros(L)
    
    for i in range(nB): h[i]      = fB
    for i in range(nD): h[-(i+1)] = fD
    

    # _ SPECIES CLASS _ #

    class Species:

        def __init__(self, n = 1, f = 1., **kwargs):
            """ Initialize clone/provirus-specific variables. """
            self.n = n   # number of members
            self.f = f   # fitness
        
            if 'sequence' in kwargs: self.sequence = kwargs['sequence']  # sequence identifider
            else:                    self.sequence = np.zeros(L)
    
        @classmethod
        def clone(cls, s):
            return cls(n = 1, f = s.f, sequence = [k for k in s.sequence]) # Return a new copy of the input Species
    
        def mutate(self):
            """ Mutate and return self + new sequences. """
        
            # get number of individuals that mutate
            newSpecies = []
            
            if self.n>0:
                nMut    = np.random.binomial(self.n, mu * L)
                self.n -= nMut
        
                # process mutations
                site = np.random.randint(L, size = nMut)
                for i in site:
                    s             = Species.clone(self)
                    s.sequence[i] = 1 - s.sequence[i]
                    s.f           = 1. + np.sum(h * s.sequence)
                    newSpecies.append(s)
    
            # return the result
            if (self.n>0): newSpecies.append(self)
            return newSpecies
    
    # ^ SPECIES CLASS ^


    # Trial length and recording frequency
    
    tStart = 0       # start generation
    tEnd   = T       # end generation
    record = 1       # record data every (record) generations
    start  = timer() # track running time

    # Create species and begin recording
    
    pop, sVec, nVec = [], [], []
    
    if randomize>0:
        if randomize>N: randomize = N
        n_seqs    = int(N/randomize)
        temp_sVec = []
        temp_nVec = []
        i         = 0
        while (i * n_seqs < N):
            temp_seq = np.random.randint(0, 2, size = L)
            pop.append(Species(n = n_seqs, sequence = temp_seq))
            temp_sVec.append(temp_seq)
            temp_nVec.append(n_seqs)
            i += 1
        temp_pop = np.sum(temp_nVec)
        if temp_pop<N:
            temp_seq = np.random.randint(0, 2, size = L)
            pop.append(Species(n = N-temp_pop, sequence = temp_seq))
            temp_sVec.append(temp_seq)
            temp_nVec.append(N-temp_pop)
        sVec.append(np.array(temp_sVec))
        nVec.append(np.array(temp_nVec))
    
    else:
        pop  = [Species(n = N)]             # current population
        sVec = [np.array([np.zeros(L)])]    # array of sequences at each time point
        nVec = [np.array([N])]              # array of sequence counts at each time point

    # Evolve the population
        
    for t in range(tStart, tEnd):
        
        printUpdate(t, tEnd)    # status check
        
        # Select species to replicate
        
        r = np.array([s.n * s.f for s in pop])
        p = r / np.sum(r)
        n = np.random.multinomial(N, pvals = p)
            
        # Update population size and mutate
        
        newPop = []
        for i in range(len(pop)):
            pop[i].n = n[i]
            p        = pop[i].mutate()
            for j in range(len(p)):
                unique = True
                for k in range(len(newPop)):
                    if np.array_equal(p[j].sequence, newPop[k].sequence):
                        unique       = False
                        newPop[k].n += p[j].n
                        break
                if unique:
                    newPop.append(p[j])
        pop = newPop
        
        # Update measurements
        
        if t%record==0:
            nVec.append(np.array([s.n        for s in pop]))
            sVec.append(np.array([s.sequence for s in pop]))

    # End and output total time
    
    f = open(out_str+'.npz', 'w')
    np.savez_compressed(f, nVec=nVec, sVec=sVec)
    f.close()
    
    #f = open('examples/figx2-N_1e3-mu_1e-3.dat', 'w')
    #for i in range(len(nVec)):
    #    for j in range(len(nVec[i])):
    #        f.write('%d\t%d\t%s\n' % (i*record, nVec[i][j], ' '.join([str(int(k)) for k in sVec[i][j]])))
    #f.close()

    end = timer()
    print('\nTotal time: %lfs, average per generation %lfs' % ((end - start),(end - start)/float(tEnd)))


def printUpdate(current, end, bar_length=20):
    """ Print an update of the simulation status. h/t Aravind Voggu on StackOverflow. """
    percent = float(current) / end
    dash    = ''.join(['-' for k in range(int(round(percent * bar_length)-1))]) + '>'
    space   = ''.join([' ' for k in range(bar_length - len(dash))])

    sys.stdout.write("\rSimulating: [{0}] {1}%".format(dash + space, int(round(percent * 100))))
    sys.stdout.flush()


if __name__ == '__main__': main()

