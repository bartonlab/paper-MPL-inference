'''
Copyleft Apr 11, 2016 Arya Iranmehr, PhD Student, Bafna Lab, UC San Diego,  Email: airanmehr@gmail.com
'''
import numpy as np;
from numba import guvectorize

import Utils.Util as utl
np.set_printoptions(linewidth=40, precision=5, suppress=True)
import pandas as pd;
pd.options.display.max_rows = 20;
pd.options.display.expand_frame_repr = False;
pd.options.display.max_columns = 20
import os;            home=os.path.expanduser('~') +'/'
import sys;sys.path.insert(1,'/home/arya/workspace/bio/')
import scipy as sc
EPS=1e-320
@guvectorize(['void(float64[:], float64[:], float64[:,:])'],'(n),()->(n,n)')
def computeTransition(nu_t,N,T):
    N=int(N[0])
    fact= np.append([0],np.log(np.arange(1,2*N+1)).cumsum())
    lognu_t=np.log(nu_t+EPS)
    lognu_tbar=np.log(1-nu_t+EPS)
    for i in range(T.shape[0]):
        for j in range(T.shape[0]):
            jj=int(np.round(nu_t[j]*2*N))
            T[i,j]= np.exp(fact[-1] - fact[jj] - fact[2*N-jj]+ lognu_t[i]*jj +  lognu_tbar[i]*(2.*N-jj))
        if not nu_t[i]: T[i, 0] = 1;
        if nu_t[i] == 1: T[i, -1] = 1;

class Markov:
    @staticmethod
    def computePower(T,n,takeLog=False):
        Tn=T.copy(True)
        for i in range(n-1):
            Tn=Tn.dot(T)
        if takeLog:
            return Tn.applymap(np.log)
        else:
            return Tn

    @staticmethod
    def normalize(a): return (a.T/a.sum(1)).T

    @staticmethod
    def computeTransition(s, N, n, h=0.5, TNeutral=None):
        def getRow(x):
            if x.low==x.up:
                return TNeutral.iloc[int(x.low)]
            else:
                theta=[x.up-x.Nt,x.Nt-x.low]
                assert sum(theta)==1
                return theta[0]*TNeutral.iloc[int(x.low)]+theta[1]*TNeutral.iloc[int(x.up)]
        nu0=np.arange(2*n+1)/float(2*n)
        if TNeutral is None:
            TNeutral=Markov.normalize(pd.DataFrame(computeTransition(nu0,N),index=nu0,columns=nu0)).fillna(0)
        if not s:
            T=TNeutral
        else:
            Nt=pd.Series(map(lambda x: 2*n*max(min(utl.fx(x, s, h=h), 1.), 0.), nu0),index=nu0).rename('Nt').round(int(np.ceil(np.log10(2*n))))
            a=pd.concat([Nt,Nt.apply(np.floor).rename('low').astype(int),Nt.apply(np.ceil).rename('up').astype(int)],1)
            T=a.groupby(level=0).apply(lambda x: getRow(x.loc[x.name]))
        return Markov.normalize(T)
    @staticmethod
    def power_recursive(T, n, powers_cached):
        if n not in powers_cached.index:
            if n % 2 == 0:
                TT = Markov.power_recursive(T, n / 2,powers_cached)
                powers_cached[n]= TT.dot(TT)
            else:
                powers_cached[n]= T .dot( Markov.power_recursive(T, n - 1,powers_cached))
        return powers_cached[n]
    @staticmethod
    def Powers(T,powers):
        powers_cached =pd.Series([np.eye(T.shape[0]),T],index=[0,1])
        for n in powers:
            Markov.power_recursive(T, n, powers_cached)
        return powers_cached.loc[powers]

    @staticmethod
    def computeProb(X,T):
        return sum([np.log(T.loc[X[t,r],X[t+1,r]]) for t in range(X.shape[0]-1) for r in range(X.shape[1])])
class Binomial:
    @staticmethod
    def computeTransition(N,nu0):
        lognu_t=np.log(nu0)
        lognu_tbar=np.log(1-nu0)
        nu_t=np.arange(2*N+1)/(2.*N)
        logrange= np.log(np.arange(1,2*N+1))
        lograngesum=logrange.sum()
        T=[]
        for j in range(nu_t.shape[0]):
            jj=int(np.round(nu_t[j]*2*N,int(np.ceil(np.log10(2*N)))))
            T+= [np.exp(lograngesum - logrange[:jj].sum() - logrange[:2*N-jj].sum()+ lognu_t*jj +  lognu_tbar*(2.*N-jj))]
        return pd.Series(T,index=nu_t)
    @staticmethod
    def computeTransitionExact(N,nu0):
        from scipy.stats import binom
        nu_t=np.arange(2*N+1)/(2.*N)
        rv = binom(2*N, nu0)
        return pd.Series(map(lambda x: rv.pmf(x*2*N),nu_t),index=nu_t)
    @staticmethod
    def computeTransitionMatrixExact(N):
        nu_t=np.arange(2*N+1)/(2.*N)
        return pd.concat(map(lambda x: Binomial.computeTransitionExact(N,x).rename(x),nu_t),1).T

    @staticmethod
    def sampling(N,n):
        T=np.zeros((2*N+1,2*n+1))
        nu_t=np.arange(2*N+1)/float(2*N)
        y_t=np.arange(2*n+1)/float(2*n)
        logrange= np.log(np.arange(1,2*n+1))
        lograngesum=logrange.sum()
        lognu_t=np.log(nu_t+EPS)
        lognu_tbar=np.log(1-nu_t+EPS)
        for i in range(T.shape[0]):
            for j in range(T.shape[1]):
                T[i,j]= np.exp(lograngesum - logrange[:j].sum() - logrange[:2*n-j].sum()+ lognu_t[i]*j +  lognu_tbar[i]*(2.*n-j))
            if not nu_t[i]: T[i, 0] = 1;
            if nu_t[i] == 1: T[i, -1] = 1;
        return pd.DataFrame(T,index=nu_t,columns=y_t)

    @staticmethod
    def likelihood(cd, nu):
        c, d = cd
        p = sc.misc.comb(d, c) * (nu ** c) * ((1 - nu) ** (d - c));
        return p

class HMM:
    def __init__(self, eps=1e-1, CD=None, CDfname=None, path=None, verbose=1,
                 N=None, n=None, Ns=None,
                 gridH=[0.5,5],stepS=0.05, nSteps=20,maxS=1,
                 loadCDE=False,saveCDE=False,transitionsPath=None,batchSize=int(2e5),
                 precomputeTransitions=False,filterOutlierReplicate=0):

        if path is not None:utl.mkdir(path)
        self.CDfname=CDfname
        if CDfname is not None: self._CD=pd.read_pickle(self.CDfname);
        else: self._CD=CD
        self.filterOutlierReplicate=filterOutlierReplicate
        self.batchSize=batchSize
        self.path,self.gridH,self.stepS,self.eps,self.verbose=path,gridH,stepS,eps,verbose
        self.n=n;self.N=N;self.transitionsPath=transitionsPath;self.Ns=Ns;self.maxS=maxS;self.nSteps=nSteps
        if self.n is None: self.n=self.N
        if self.Ns is None: self.Ns=self.n
        if self.CDfname is None:
            path=self.path
        else:
            path=self.CDfname.replace('.df','.')

        if self.Ns is not None and self.n is not None:
            print(self.Ns,self.n)
            self.CD,self.E=HMM.precomputeCDandEmissions(CD=self._CD, n=self.Ns, N=self.n,path=path,loadCDE=loadCDE,saveCDE=saveCDE ,verbose=self.verbose)
        R=self._CD.columns.get_level_values('REP').unique()
        self.powers = pd.Series([pd.Series(self._CD[r].columns.get_level_values('GEN').unique()).diff().values[1:] for r in R],index=R)
        self.setStepS()
        if self.transitionsPath is None:
            if self.path is not None:
                self.transitionsPath=self.path+'T/';utl.mkdir(self.transitionsPath)
        if precomputeTransitions:
            self.computeTransitions()

    def setStepS(self):
        if self.maxS is None:
            self.maxS=self.findMaxS()
            self.stepS=self.maxS/self.nSteps
        elif self.stepS is None:
            self.stepS=self.maxS/self.nSteps

    def findMaxS(self):
        freq=lambda x: x.xs('C',level='READ',axis=1).sum(1)/x.xs('D',level='READ',axis=1).sum(1)
        x=self._CD.groupby( level='GEN',axis=1).apply(lambda x: freq(x)).sort_index(1)
        x[(x==1)|(x==0)]=None; x=x.dropna()
        s= (2./(x.columns[-1]-x.columns[0])*(utl.logit(x.iloc[:,-1])-utl.logit(x.iloc[:,0])).abs().replace(np.inf,None).dropna()).max()
        s='{:e}'.format(s)
        try:
            s=(int(s.split('.')[0])+1)*10**(-int(s.split('-')[1]))
        except:
            s=(int(s.split('.')[0])+1)*10**(int(s.split('+')[1]))
        return s
    def likelihood(self,s,CD):
        """
        Args: (it's more convenient for multiprocessing)
            args: a list of [R,s,h].
            R: is a dataframe for which each row is a position and columns are allele frequencies.
                ColumnsLevels= [REP, TIME] , IndexLevels=[CHROM,POS]
            s: is selection strength
            h: is overdominance
        Returns:
            a series containing likelihood of timeseries for the specific values of s and h.
        """
        try:
            if not s: return self.likes_null
        except: pass
        if not CD.shape[0]: return pd.Series()
        if s==0:s=int(0)
        if self.verbose>0:print('Computing for {} SNPs for s={} h={}'.format(CD.shape[0], s, self.h));sys.stdout.flush()
        try:
            T = pd.read_pickle(self.transitionsPath + 'N{}.S{:E}.H{:E}.df'.format(self.N,s, self.h))
        except:
            T=HMM.precomputeTransitions(((self.CD, (s, self.h), self.N, self.n, self.transitionsPath, None,self.verbose)))
        args = map(lambda x: (x, self.E, T, self.powers), utl.batch(CD, self.batchSize));
        f=(HMM.likelihoodBatch,HMM.likelihoodBatchFilterRep)[self.filterOutlierReplicate>0]
        likes = pd.concat(map(f, args)).rename((s, self.h))
        if self.verbose>1: print(pd.DataFrame(likes))
        return likes

    def likelihoodN(self,N,n):
        if n>N: n=N
        if self.verbose>0:
            print('Computing for N={}, n={}'.format(N,n));
        try:
            T = pd.read_pickle(self.transitionsPath + 'N{}.df'.format(N))
            T=pd.Series([T],index=[10])
        except:
            T=HMM.precomputeTransitions(((self._CD, (0, 0.5), N, n, None, None,self.verbose)))
        self.CD,self.E=HMM.precomputeCDandEmissions(CD=self._CD, n=n, N=n,loadCDE=False,saveCDE=False,verbose=self.verbose)
        args = map(lambda x: (x, self.E, T, self.powers), utl.batch(self.CD, self.batchSize));
        likes = pd.concat(map(HMM.likelihoodBatch, args)).rename(N)
        return likes

    @staticmethod
    def likelihoodBatch(args):
        CD, E, T, powers = args
        likes = pd.Series(0, index=CD.index)
        n=0
        startGen=CD.columns.get_level_values('GEN').min()
        for rep, df in CD.T.groupby(level=0):
            alpha = E.iloc[df.loc[(rep, startGen)]].values
            for step, power in zip(range(1, df.shape[0]), powers[rep]):
                alpha = alpha.dot(T.loc[power].values) * E.values[df.loc[rep].iloc[step].values]
            likes += utl.vectorizedLog(alpha.mean(1)) #it should be here
            n+=1
        return likes/n

    @staticmethod
    def likelihoodBatchFilterRep(args):
        CD, E, T, powers,filterN = args
        dfl=[]
        for rep, df in CD.T.groupby(level=0):
            alpha = E.iloc[df.loc[(rep, 0)]].values
            for step, power in zip(range(1, df.shape[0]), powers[rep]):
                alpha = alpha.dot(T.loc[power].values) * E.values[df.loc[rep].iloc[step].values]
            dfl+=[utl.vectorizedLog(alpha.mean(1))]
        df= pd.DataFrame(dfl,columns=CD.index)
        return df.apply(lambda x: x.sort_values()[1:].mean())


    def likelihoods(self,rangeH=[0.5],rangeS=np.arange(-0.5,0.5,0.02)):
        res=[]
        for h in rangeH:
            for s in rangeS:
                self.h=h
                res+=[self.likelihood(s ,self.CD).rename((h,s))]
        return pd.concat(res,1)

    def fitOne(self,h):
        self.h=h
        likes_null = self.likelihood(0 ,self.CD).rename('null');
        dfn = self.linesearch(likes_null, False)
        dfp = self.linesearch(likes_null, True)
        I=dfp.lik>=dfn.lik
        df = pd.concat([dfp[I],dfn[~I]])
        return pd.concat([df.rename(columns={'lik':'alt'}),likes_null],1)



    def fitN(self,rangeN=np.arange(1,15,1)*100,n=1000):
        return pd.concat(map(lambda x: self.likelihoodN(x,n),rangeN),1,keys=rangeN)

    def fitNLineSearch(self,rangeN=np.arange(1,15,1)*100,n=1000):
        likes=pd.Series(None)
        prev=-1e10
        for N in rangeN:
            likes.loc[N]=self.likelihoodN(N,n).mean()
            if likes.loc[N]<prev:
                return likes
            prev=likes.loc[N]
        return likes


    def fit(self,save):
        df=pd.concat(map(self.fitOne,self.gridH),1,keys=self.gridH)
        df.columns.names=['h','stat']
        if save:
            if self.CDfname is None:
                fname=self.path+ 'HMM.df'
            else:
                fname=self.CDfname.replace('.df','.HMM.df')
            df.to_pickle( fname)
        return df

    def maximumLikelihood(self, scores,save):
        print('computing scores....')
        a = scores.groupby(level=0, axis=1).apply(lambda x:x[(x.name, 'alt')] - x[(x.name, 'null')])
        h = a.abs().apply(lambda x: x.idxmax(), axis=1)
        sh = scores.groupby(level=[0, 1]).apply(lambda x: x[h.loc[x.name]].s.values[0])
        lrdiff = a.max(1)-a[0.5]
        df = pd.concat([a[0.5].abs(), scores[(0.5, 's')], h, a.abs().max(1), sh, lrdiff], axis=1)
        df.columns = ['lrDrectional', 'sDirectional', 'h', 'lr', 's', 'lrdiff']
        if save is not None:df.to_pickle(self.path + 'scores.df')
        return df

    def linesearch(self,init,PositiveS):
        sgn=[-1,1][PositiveS]
        S = np.arange(0,  sgn*self.maxS+1e-10,  sgn*self.stepS)[1:]
        i = pd.Series(True, index=init.index).values;
        mlprev = init.values.copy(True);
        mlcurrent = init.values.copy(True)
        mle = np.zeros(mlcurrent.size)
        ml = init.values.copy(True)
        for s in S:
            mlprev[i] = mlcurrent[i]
            mlcurrent[i] = self.likelihood(s,self.CD[i])
            i = mlcurrent > mlprev + self.eps
            sys.stdout.flush()
            if i.sum() == 0: break
            mle[i] = s
            ml[i] = mlcurrent[i]
        return pd.DataFrame([ml, mle], index=['lik', 's'], columns=self.CD.index).T


    def SHgrid(self,gridH ):
        S = np.arange(-self.maxS, self.maxS+0.0001, self.stepS)
        SS,HH=np.meshgrid(S,gridH)
        SH=zip(SS.reshape(-1), HH.reshape(-1))
        return SH


    def computeTransitions(self,numProc=8):
        utl.mkdir(self.transitionsPath)
        print('Computing Transitions for Real Data...')
        args=map(lambda sh: (self.CD,sh,self.N,self.n,self.transitionsPath,None,self.verbose), self.SHgrid(gridH=self.gridH))
        # Pool(8).\
        map(precomputeHelper,args)


    @staticmethod
    def precomputeCDandEmissions(CD,n,N,path=None,saveCDE=False,loadCDE=None,verbose=1):
        """
        0- reads C read counts of reference  and D counts of depth
        1- computes alternate allele reads based on reference and depth
        2- saves CD
        3- saves state conditional distributions P(nu|(c,d)) aka emissions
        """
        if loadCDE:
            try:
                return pd.read_pickle(path + 'CDEidx.df'),pd.read_pickle(path + 'E.df')
            except:
                pass
        if verbose>0:
            print('Precomputing CD (C,D)=(Derived count,total Count) and corresponding emission probabilities...',CD.shape)
        nu=pd.Series(np.arange(0, 1.0000001, 1./(2*n)), index=np.arange(0, 1.00001, 1./(2*n)))
        c = CD.xs('C', level='READ', axis=1)
        d = CD.xs('D', level='READ', axis=1)
        cd = pd.concat([pd.Series(list(zip(c[i], d[i]))) for i in c.columns], axis=1);
        cd.columns = c.columns;
        cd.index = c.index
        allreads = pd.Series(cd.values.reshape(-1)).unique();
        allreads = pd.Series(allreads, index=pd.MultiIndex.from_tuples(allreads, names=['c', 'd'])).sort_index()
        E= allreads.apply(lambda x: Binomial.likelihood(x, nu)).sort_index().fillna(0)
        if n!=N:
            Y=Binomial.sampling(N=N,n=n)
            E=pd.DataFrame(E.values.dot(Y.T.values),index=E.index,columns=Y.index)
        index = pd.Series(range(E.shape[0]), E.index).to_dict()
        CDEidx = cd.applymap(lambda x: index[x])
        if saveCDE:
            E.to_pickle(path + 'E.df')
            CDEidx.to_pickle(path + 'CDEidx.df')
        return CDEidx,E


    @staticmethod
    def Powers(CD):
        if 'READ' in CD.columns.names:
            g=CD.xs('D',level=2,axis=1)
        else:
            g=CD
        powers=g.groupby(level=0,axis=1).apply(lambda x: pd.Series(x[x.name].columns).diff().values[1:].astype(int))
        return np.unique(np.concatenate(powers.tolist()))

    @staticmethod
    def precomputeTransitions(args):
        CD, sh, N, n, path,powers,verbose=args
        s,h=sh
        if CD is not None:powers=HMM.Powers(CD)

        Tn=Markov.Powers(Markov.computeTransition(s, N, n, h=h).fillna(0),powers)
        if verbose>0:
            print('Computing Transition for s={}, h={}, N={}, n={}'.format( s,h,N,n))
        if path is not None:
            Tn.to_pickle('{}N{}.S{:E}.H{:E}.df'.format(path, N,s, h))
        return Tn

def precomputeHelper(args):
    HMM.precomputeTransitions(args)


def likelihoodsN(CD, rangeN, n=200, minAF=0.01, k=2000):
    print('Performing grid search on N=',rangeN)
    cd=utl.polymorphic(CD,index=False,minAF=minAF)
    I=np.random.choice(cd.shape[0],k)
    a=HMM(CD=cd.iloc[I],gridH=[0.5],verbose=-1).fitN(rangeN=rangeN,n=n).sort_index(1).mean(0).rename('Likelihood')
    a.index.name='N'
    return a
def estimateN(CD,Nt=False,Nc=False,Nr=False,name='',rangeN=None):
    if Nc:
        return CD.groupby(level='CHROM').apply(lambda x: estimateN(x,rangeN=rangeN,Nt=Nt,Nr=Nr,name=name+'Chromosome {} '.format(x.name)))
    if Nr:
        return CD.groupby(level='REP',axis=1).apply(lambda x: estimateN(x,rangeN=rangeN,Nt=Nt,name=name+'Replicate {} '.format(x.name)))
    if Nt:
        gens=sorted(map(int,CD.columns.get_level_values('GEN').unique()))
        a=pd.concat([estimateN(CD.loc[:,pd.IndexSlice[:,[i,j]]],rangeN=rangeN, name=name+'Between F{}-F{} '.format(i,j)) for i,j in zip(gens[:-1],gens[1:])],1,keys=gens[1:])
        return pd.DataFrame(a)
    if name is not None: print('\nEstimating N for',name)
    if rangeN is None:
        a=likelihoodsN(CD, rangeN=10 ** np.arange(2, 7))
        print(a)
        N=a.idxmax()
        rangeN=np.append(np.linspace(N/10,N,10),np.linspace(N,N*10,10)[1:])
        b=likelihoodsN(CD, rangeN=rangeN);
        b.loc[a.idxmax()]=a.max();b.sort_index(inplace=True)
    else:
        b=likelihoodsN(CD, rangeN=rangeN);
    return b
