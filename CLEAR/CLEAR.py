'''
Copyleft Feb 22, 2017 Arya Iranmehr, PhD Student, Bafna Lab, UC San Diego,  Email: airanmehr@gmail.com
'''
import os,sys,optparse
import numpy as np;
import pandas as pd;
import seaborn as sns
import pylab as plt;
import matplotlib as mpl
import Utils.Util as utl
import Utils.Plots as pplt
import Libs.Markov as mkv
#except:
#    import CLEAR.Libs.Markov as mkv
sys.path.insert(1,os.getcwd())
np.set_printoptions(linewidth=200, precision=5, suppress=True)
pd.options.display.max_rows = 20;
pd.options.display.expand_frame_repr = False


parser = optparse.OptionParser()
parser.add_option( '--sync', action="store", dest="syncFile", help="path to synchronized file created by popoolation2")
parser.add_option( '--pandas', action="store", dest="pandasFile", help="path to pandas dataframe")
parser.add_option( '--vcfgz', action="store", dest="vcfgzFile", help="path to vcfgz file. Each sample in the VCF file represetn a population. Sample name should be in the format RXFY where X is the replicate id and Y is the generation.")
parser.add_option( '--N', default=0, action="store", dest="N", help="population size. otherwise estimates population size.")
parser.add_option( '--Nt', default=False,action="store_true", dest="Nt", help="estimtes N in time for every consucutive pairs of samples")
parser.add_option( '--Nc', default=False,action="store_true", dest="Nc", help="estimates N for each chromosome separately")
parser.add_option( '--Nr', default=False,action="store_true", dest="Nr", help="estimates N for each replicate separately")
parser.add_option( '--alpha', action="store", dest="alpha", help="FDR cutoff for single loci test")
parser.add_option( '--Alpha', action="store", dest="alpha", help="FDR cutoff for window based test")
parser.add_option( '--map', action="store", dest="gmap", help="Genetic map")
parser.add_option( '--O', default='d', action="store", dest="O", help="output format")
parser.add_option( '--plot',default=False, action="store_true", dest="plot", help="output format")
parser.add_option( '--out',default=None, action="store", dest="out", help="output dataframe path")

options, args = parser.parse_args()



if __name__ == '__main__':
    if options.pandasFile is not None:
        CD=pd.read_pickle(options.pandasFile)
    elif options.syncFile is not None:
        CD=utl.SynchronizedFile.load(options.syncFile)
    elif options.vcfgz is not None:
        CD=utl.VCF.loadCD(options.vcfgz)
    else:
        print('Invalid input')
        exit()
    n=200
    if options.N:
        N=int(options.N)
    else:
    	a= mkv.estimateN(CD,Nt=options.Nt,Nc=options.Nc, Nr=options.Nr)
    	N=a.idxmax()
    	a=a.reset_index();a.columns=['N','Likelihood'];print(a)
    	print('Maximum Likelihood of N=',N)
 	
    HMM=mkv.HMM(eps=1e-2,CD=CD,gridH=[0.5],N=N,n=n,saveCDE=False,loadCDE=False,verbose=1,maxS=None)
    a= HMM.fit(False)
    print(a)
    if options.out is not None:
        a.to_pickle(options.out)
        print('Output is saved in pandas dataframe in {}.'.format(options.out))

    if options.plot:
        f=lambda x: x.alt-x.null
        a=f(a[0.5])
        fig,axes=plt.subplots(2,1,sharex=True,dpi=200)
        pplt.Manhattan(a.rename('$H$'),top_k=10,axes=[axes[0]])
        pplt.Manhattan(utl.scanGenome(a).rename(r'$\mathcal{H}$'),top_k=3,axes=[axes[1]])
        plt.show()

