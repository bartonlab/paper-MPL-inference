'''
Copyleft Oct 10, 2015 Arya Iranmehr, PhD Student, Bafna's Lab, UC San Diego,  Email: airanmehr@gmail.com
'''

import sys

import numpy as np
import pandas as pd
import pylab as plt

try:
    from SFSelect import metaSVM
    sys.modules['metaSVM'] = metaSVM
    #Models are loaded once.
    from SFSelect.SFselect import sfselect,SVM,SVMXP
except:
    SVM,SVMXP= None,None

class Estimate:
    @staticmethod
    def pi(snp):
        n=snp.shape[0]
        return np.mean([sum(np.logical_xor(snp[i,:],snp[j,:])) for i in range(n-1) for j in range(i+1,n)])

    @staticmethod
    def getSAFS(x,bins=10, n=None,fixedRangeHist=True,removeFixedSites=False,fold=False,normed=False,noBinSmalln=False):
        afs=Estimate.getAFS(x, bins, n=n, fixedRangeHist=fixedRangeHist,removeFixedSites=removeFixedSites,normed=normed,fold=fold,noBinSmalln=noBinSmalln)
        if afs is not None:
            return afs*afs.index.values
        else:
            return  None
    
    @staticmethod
    def getAFS(freqs,bins=10,n=None,fixedRangeHist=True,removeFixedSites=True,name=None,normed=False,fold=False,noBinSmalln=False):
        if removeFixedSites:
            x=freqs[(freqs>0) & (freqs<1)].copy(True)
        else:
            x=freqs.copy(True)
        if noBinSmalln:
            nn=1./x.value_counts().sort_index().index[0]
            if nn <2*bins:bins=-1
            if n is None: n=nn
        if not x.size:
            if n is not None:
                return pd.Series(0,index=range(1,n))
        if fold:
            x[x>0.5]=1-x[x>0.5]
            if n/2 <2*bins: bins=-1
        if bins>0:
            if fixedRangeHist:
                rang=[0,1]
                if x.max()<=0.5: rang=[0,0.5]
                counts,limits=np.histogram(x, bins,range=rang) # counts=xi
            else:
                counts,limits=np.histogram(x, bins) # counts=xi
            centers = 0.5*(limits[1:]+limits[:-1])  # center = i
            if n is not None:
                centers=centers*n
            afs= pd.Series(counts, index=centers)  #  site frequency spectrum
        else:
            afs=pd.Series(x.round(3)).value_counts()
            if not len(afs.index): return None
            if n is not None:
                afs.index=np.round(afs.index*n)

            afs=afs.sort_index()
            if 0 in  afs.index:
                if 1 in  afs.index:
                    afs.loc[1]+=afs.loc[0]
                else:
                    afs.loc[1]=afs.loc[0]
                afs=afs.drop(0).sort_index()
        if name is not None:
            afs.name=name

        if bins==-1:
            afs=(afs+pd.Series(0,range(1,(n,int(n/2)+1)[fold]))).fillna(0)
            afs.index=map(int,afs.index)
        if normed:
            return afs/float(afs.sum())
        return  afs
    @staticmethod
    def plotSAFS(x,bins=10,n=None, fixedRangeHist=True,removeFixedSites=False):
        Estimate.getSAFS(x, bins=bins, n=n, fixedRangeHist=fixedRangeHist,removeFixedSites=removeFixedSites).plot(kind='bar',width=1,grid=True);plt.xlim([-0.5,bins-0.5])
    
    @staticmethod
    def getWeights(safs,method,n):
        i=safs.index.values*1.0
        if method is 'watterson':
            w=1/i/(1./np.arange(1,n)).sum()
        elif method is 'pi':
            w=(n-i) /((n*(n-1)/2.))
        elif method is 'faywu':
            w=i / (n*(n-1)/2.)
        return w
    @staticmethod
    def getEstimate(x=None,x_XP_pop=None, n=None, snp=None, method='all', bins=-1, normalizeTajimaD=True, averageResults=False,
                    svm_model_sfselect=SVM, fixedRangeHist=True, removeFixedSites=True, selectionPredictor=False):
        """
        Compute different estimates either based on AF (x and n should be given) or SNP matrix (only SNP matrix suffices) 
        watterson: watterson's theta 
        pi: tajima's pi (avg heterozygosity) 
        faywu: FayWu estimate of theta
        tajimaD: pi-watterson ***********************FOR THIS CASE -TAJIMAD IS RETURNED TO BE A PREDICTOR OF SELECTION, (NEG TAJIMAD MEANS SELECTION)
        H: pi-faywu           ***********************FOR THIS CASE -H IS RETURNED TO BE A PREDICTOR OF SELECTION, (NEG H MEANS SELECTION)
        
        Parametes:
        x: can be either tensor, matrix or vector of AFS with dimension T x L x R, or L x R or L
        n: number of samples
        snp: snp matrix in dataframe format
        """
        if x is not None:   # if AF is given
            if len(x.shape)==3: # if x is tensor T x L x R
                return pd.DataFrame([Estimate.getEstimate(x=xt, method=method, bins=bins, n=n, normalizeTajimaD=normalizeTajimaD, averageResults=averageResults, svm_model_sfselect=svm_model_sfselect, fixedRangeHist=fixedRangeHist,removeFixedSites=removeFixedSites, selectionPredictor=selectionPredictor) for xt in x])
            elif len(x.shape)==2: # if x is matrix L x R
                if averageResults:
                    return np.mean([Estimate.getEstimate(x=xr, method=method, bins=bins, n=n, normalizeTajimaD=normalizeTajimaD, averageResults=averageResults, svm_model_sfselect=svm_model_sfselect, fixedRangeHist=fixedRangeHist,removeFixedSites=removeFixedSites, selectionPredictor=selectionPredictor) for xr in x.T])
                else:
                    return ([Estimate.getEstimate(x=xr, method=method, bins=bins, n=n, normalizeTajimaD=normalizeTajimaD, averageResults=averageResults, svm_model_sfselect=svm_model_sfselect, fixedRangeHist=fixedRangeHist,removeFixedSites=removeFixedSites, selectionPredictor=selectionPredictor) for xr in x.T])
            elif len(x.shape)==1: # if x is L-dim vector of AF
                if method=='SFSelect':
                    if x_XP_pop is not None:
                        svm_model_sfselect = SVMXP
                    return sfselect(x,neut_pop_freqs=x_XP_pop, svm=svm_model_sfselect,removeFixedSites=removeFixedSites)['score']
                safs=Estimate.getSAFS(x=x, bins=bins, n=n, fixedRangeHist=fixedRangeHist,removeFixedSites=removeFixedSites)
        elif snp is not None:  # if SNP is given
            safs=snp.sum(0).value_counts().sort_index()
            safs=safs[safs.index!=0]
            safs = safs* safs.index.values
            n=snp.shape[0]
        else:
            return None

        if method=='all':
            return ({
                'w': Estimate.getEstimate(x=x, n=n, snp=snp, method='watterson',
                                          bins=bins, normalizeTajimaD=True,
                                          averageResults=averageResults,
                                          svm_model_sfselect=svm_model_sfselect,
                                          fixedRangeHist=fixedRangeHist,
                                          removeFixedSites=removeFixedSites,
                                          selectionPredictor=selectionPredictor),
                'pi': Estimate.getEstimate(x=x, n=n, snp=snp, method='pi',
                                           bins=bins, normalizeTajimaD=True,
                                           averageResults=averageResults,
                                           svm_model_sfselect=svm_model_sfselect,
                                           fixedRangeHist=fixedRangeHist,
                                           removeFixedSites=removeFixedSites,
                                           selectionPredictor=selectionPredictor),

                '{}D'.format(('', '-')[selectionPredictor]): Estimate.getEstimate(x=x, n=n, snp=snp, method='tajimaD',
                                                                                  bins=bins, normalizeTajimaD=True,
                                                                                  averageResults=averageResults,
                                                                                  svm_model_sfselect=svm_model_sfselect,
                                                                                  fixedRangeHist=fixedRangeHist,
                                                                                  removeFixedSites=removeFixedSites,
                                                                                  selectionPredictor=selectionPredictor),
                '{}H'.format(('', '-')[selectionPredictor]): Estimate.getEstimate(x=x, n=n, snp=snp, method='H',
                                                                                  bins=bins, normalizeTajimaD=True,
                                                                                  averageResults=averageResults,
                                                                                  svm_model_sfselect=svm_model_sfselect,
                                                                                  fixedRangeHist=fixedRangeHist,
                                                                                  removeFixedSites=removeFixedSites,
                                                                                  selectionPredictor=selectionPredictor),
                'SFSelect': Estimate.getEstimate(x=x, n=n, snp=snp, method='SFSelect', bins=bins, normalizeTajimaD=True,
                                                 averageResults=averageResults, svm_model_sfselect=svm_model_sfselect,
                                                 fixedRangeHist=fixedRangeHist, removeFixedSites=removeFixedSites,
                                                 selectionPredictor=selectionPredictor)})

        if safs is None:return  None
        if method is 'tajimaD':
            w = Estimate.getWeights(safs, 'pi', n) - Estimate.getWeights(safs, 'watterson', n)
            if normalizeTajimaD:
                if x is not None:
                    m = len(x)
                if snp is not None:
                    m = snp.shape[1]
                sig = Estimate.tajimaDstd(n=n, m=m)
                w /= sig
            if selectionPredictor:  w*=-1 
        elif method is 'H':
            w = Estimate.getWeights(safs, 'pi', n) - Estimate.getWeights(safs, 'faywu', n)
            if selectionPredictor:  w*=-1
        else:
            w=Estimate.getWeights(safs, method, n)
        return safs.dot(w)
        
    @staticmethod
    def tajimaDstd(n,m):
        a1=(1./np.arange(1,n)).sum();       a2=(1./(np.arange(1,n)**2).sum())
        b1=(n+1.0)/(3*n-3);                 b2=2.0*(n*n+n+3)/(9*n*n-9*n)
        c1=b1-1.0/a1;                       c2=b2 -(n+2)/(a1*n) +a2/(a1*a1)
        e1=c1/a1;                           e2=c2/(a1*a1+a2)
        return np.sqrt(e1*m+e2*m*m-e2*m)
    
    @staticmethod
    def watterson(snp=None,n=None,m=None):
        """
        computes watterson estimate of theta
        snp: is m x n numpy matrix
        n: number of samples
        m: number of segregating sites
        """
        if n is None or m is None:
            n,m=snp.shape
        return m/(1./np.arange(1,n)).sum()
    
    @staticmethod
    def rho(snp):
        return snp.shape[1]/np.log(snp.shape[0])
    
    @staticmethod
    def getAllEstimates(snp):
        thetaw=Estimate.watterson(snp)
        pi=Estimate.pi(snp)
        tajD=(pi-thetaw)#/Estimate.tajimaDstd(n=snp.shape[0],m=snp.shape[1])
        n=snp.shape[0]
        x=snp.mean(0)
        fay=2.*n/(n-1)*np.linalg.norm(x)**2
        sfsel = Estimate.getEstimate(x=x, n=snp.shape[0], method='SFSelect')
        return pd.DataFrame([('Theta', thetaw),('Pi',pi),('TajimaD',tajD),('m',int(snp.shape[1])), ('FayWu', fay),('SFSelect', sfsel)],columns=['method','estimate'])

    @staticmethod
    def getAllEstimatesX(X, n=200, bins=-1, method=None):
        x=X[(X!=0) & (X!=1)]
        m = x.size
        if not m:
            i,v=zip(*[('Theta', None), ('Pi', None), ('TajimaD', None), ('m', None), ('FayWu', None), ('SFSelect', None)])
            all = pd.Series(v,index=i,name=X.name)
            all.index.name='method'
            if method is None:
                return all
            else:
                return None
        thetaw = Estimate.watterson(n=n, m=m)
        pi = Estimate.getEstimate(x=x, n=n, method='pi', bins=bins)
        tajD = (pi - thetaw) #/ Estimate.tajimaDstd(n=n, m=m)
        fay = Estimate.getEstimate(x=x, n=n, method='H', bins=bins)
        sfsel = Estimate.getEstimate(x=x, n=n, method='SFSelect')
        i,v=zip(*[('Theta', thetaw), ('Pi', pi), ('TajimaD', tajD), ('m', int(m)), ('FayWu', fay), ('SFSelect', sfsel)])
        all = pd.Series(v,index=i,name=X.name)
        all.index.name='method'
        if method is None:
            return all
        else:
            return all.loc[method]
    @staticmethod
    def mu(S=6188,l=21647181,T=5300000):
        """S/2 = mu*l*T"""
        return S/(2.*l*T)
    @staticmethod
    def Ne(mu,S=25742,n=12,l=21647181):
        """S = mu*N*Ttot*l"""
        Ttot= sum(2./np.arange(1,n))
        print 'Ttot' , Ttot
        print 'Ttotmu' , (mu*Ttot)
        return S/(mu*Ttot*l) 
    
    
    @staticmethod
    def LDold(SNP,site=None,sites=None,positions=None,measure='DPrime'):
        """
        Computes All Measures of LD between all sites of a SNP matrix to a site 
        SNP: a pandas dataframe which its columns contains position 
        Site: index of the column which LD is computed for
        Sites: index of the columns which pairwise LD is computed for
        Positions: positions which pairwise LD is computed for
        measure: {'all', 'D','DPrime','Rho','RhoPrime', 'Rho2', 'DPrime2'}
        NOTE THAT RhoPrime=DPrime
        """
        if site is None:
            if sites is not None:
                D=pd.DataFrame(map( lambda x: Estimate.LD(SNP.iloc[:,sites],x,measure),range(sites.shape[0])));
            elif positions is not None:
                D=pd.DataFrame(map( lambda x: Estimate.LD(SNP.loc[:,positions],x,measure),range(positions.shape[0])));
            else:
                D=pd.DataFrame(map( lambda x: Estimate.LD(SNP,x,measure),range(SNP.shape[1])));
            D.index=D.columns
            return D

        LD=[]
        p0=(SNP.iloc[:,site]==0).mean()
        p1=(SNP.iloc[:,site]==1).mean()
        for i in range(SNP.shape[1]):
            q0=np.mean(SNP.iloc[:,i]==0)
            q1=np.mean(SNP.iloc[:,i]==1)
            x00=( ( SNP.iloc[:,i]==0) & (SNP.iloc[:,site]==0 ) ).mean()
            x11 = ((SNP.iloc[:, i] == 1) & (SNP.iloc[:, site] == 1)).mean()
            D = (x00 + x11) / 2. - p0 * q0
            if D<0:
                Dmax=min(p0*q0,p1*q1)
            else:
                Dmax=min(p0*q1,p1*q0)
            Dprime=D/Dmax
            denom=np.sqrt(p0*p1*q0*q1)
            if denom:
                rho=D/denom
                LD.append((D,Dprime, rho, rho/(Dmax/denom), {'p':[p0,p1],'q':[q0,q1],'x00':x00},Dmax))
            else:
                LD.append((D, Dprime, None, None) )
        LD= pd.DataFrame(LD,index=SNP.columns,columns=['D','DPrime','Rho','RhoPrime','Freq','Dmax'])
        LD['Rho2']=LD.Rho**2
        LD['DPrime2']=LD.DPrime**2
        if measure=='all':
            return LD
        else:
            return LD[measure]


    @staticmethod
    def LD(SNP,site=None,sites=None,positions=None,measure='DPrime'):
            """
            Computes All Measures of LD between all sites of a SNP matrix to a site
            SNP: a pandas dataframe which its columns contains position
            Site: index of the column which LD is computed for
            Sites: index of the columns which pairwise LD is computed for
            Positions: positions which pairwise LD is computed for
            measure: {'D','DPrime','Rho','RhoPrime'}
            NOTE THAT RhoPrime=DPrime
            """
            measures=np.array(['D','DPrime','Rho','RhoPrime'])
            measureid=np.where(measures==measure)[0][0]
            if site is None:
                if sites is not None:
                    D=pd.DataFrame(map( lambda x: LDvectorized(SNP.iloc[:,sites].values,x,measureid),range(sites.shape[0])),columns=SNP.columns[sites]);
                elif positions is not None:
                    D=pd.DataFrame(map( lambda x: LDvectorized(SNP.loc[:,positions].values,x,measureid),range(positions.shape[0])),columns=positions);
                else:
                    D=pd.DataFrame(map( lambda x: LDvectorized(SNP.values,x,measureid),range(SNP.shape[1])),columns=SNP.columns);
                D.index=D.columns
            else:
                D=pd.Series(LDvectorized(SNP.values,site,measureid),index=SNP.columns)
            return D

    @staticmethod
    def LD_usingR(H0):
        import rpy2.robjects as robjects
        from rpy2.robjects import pandas2ri
        robjects.r('options(warn=-1)')
        robjects.r('library(genetics)')
        genotype = H0.applymap(lambda x: ('A', 'G')[x == 0])
        genotype = pd.concat(
                [pd.concat([genotype.iloc[i], genotype.iloc[i + 1]], axis=1).apply(lambda x: '/'.join(x), axis=1) for i
                 in np.arange(0, H0.shape[0], 2)], axis=1).T
        c = robjects.r['LD'](robjects.r['makeGenotypes'](pandas2ri.py2ri_pandasdataframe(genotype.astype('str'))))
        c = pd.Series(map(lambda x: pd.DataFrame(pandas2ri.ri2py(x)), c[1:]), index=list(c.names[1:])).apply(
            lambda x: x.fillna(0))
        for x in c:
            x += x.T
            x.index = H0.columns;
            x.columns = H0.columns;
        c.apply(lambda x: np.fill_diagonal(x.values, None))
        return c

from numba import guvectorize
@guvectorize(['void(float64[:,:],int64[:],int64[:],float64[:])'],'(M,N),(),()->(N)')
def LDvectorized(SNP,site,measure,LD):
    """
    Computes All Measures of LD between all sites of a SNP matrix to a site
    SNP: a pandas dataframe which its columns contains position
    Site: index of the column which LD is computed for
    Sites: index of the columns which pairwise LD is computed for
    Positions: positions which pairwise LD is computed for
    measure: {'all', 'D','DPrime','Rho','RhoPrime', 'Rho2', 'DPrime2'}
    NOTE THAT RhoPrime=DPrime
    """
    # LD=np.zeros((SNP.shape[1],4) )
    site=site[0]
    p0=(SNP[:,site]==0).mean()
    p1=(SNP[:,site]==1).mean()
    for i in range(SNP.shape[1]):
        q0=np.mean(SNP[:,i]==0)
        q1=np.mean(SNP[:,i]==1)
        x00=( ( SNP[:,i]==0) & (SNP[:,site]==0 ) ).mean()

        D = (x00) - p0 * q0
        if D<0:
            Dmax=min(p0*q0,p1*q1)
        else:
            Dmax=min(p0*q1,p1*q0)
        Dprime=D/Dmax
        denom=np.sqrt(p0*p1*q0*q1)
        if denom:
            rho=D/denom
            LD[i]=[D,Dprime, rho, rho/(Dmax/denom)][measure[0]]
        else:
            LD[i]=[D, Dprime, None, None][measure[0]]