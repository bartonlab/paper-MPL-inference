'''
Copyleft Dec 17, 2015 Arya Iranmehr, PhD Student, Bafna Lab, UC San Diego,  Email: airanmehr@gmail.com
'''
import numba
import pandas as pd;
import numpy as np;



np.set_printoptions(linewidth=140, precision=5, suppress=True)
import os

try:
    import readline
    import rpy2.robjects as robjects
except:
    pass

def mkdir(path):
    if not os.path.exists(path): os.makedirs(path)
parentdir=lambda path:os.path.abspath(os.path.join(path, os.pardir))
home = os.path.expanduser('~') + '/'
paperPath = home + 'workspace/timeseries_paper/'
dataPath=home+'storage/Data/'
dataPathDmel=dataPath+'Dmelanogaster/'
dataPath1000GP=dataPath+'Human/20130502/ALL/'
paperFiguresPath = paperPath + 'figures/'
plotpath=home+'out/plots/';
outpath=home+'out/';
simoutpath=dataPath+'SimulationOutFiles/'
stdoutpath = outpath + 'std/'
mkdir(outpath)
mkdir(simoutpath)
mkdir(outpath)
mkdir(plotpath)
mkdir(stdoutpath)
from subprocess import Popen, PIPE, STDOUT
comaleName = r'\sc{Clear}'
def googleDocURL(name,url):return '=HYPERLINK("{}","{}")'.format(url,name)
import scipy.stats as st
def freq(X,fold=False):
        x=X.C/X.D
        if fold:x[x>0.5]=1-x[x>0.5]
        return x
class pval:
    #c=utl.scanGenome(a,f=lambda x: utl.chi2SampleMeanPval(x,1))
    @staticmethod
    def zscore(x): return ((x-x.mean())/x.std()).astype(float)
    @staticmethod
    def zscoreChr(x): return x.groupby(level='CHROM').apply(pval.zscore).astype(float)
    @staticmethod
    def chi2(x,df=1):return - np.log10(np.exp(st.chi2.logsf(x,df)))
    @staticmethod
    def chi2SampleMean(x,df,dftot=None):
        if not x.size: return None
        df=df*x.size
        mu=x.mean()
        if dftot is not None: df=dftot
        return -np.log10(np.exp(st.gamma.logsf(mu,df/2.,scale=2./x.size)))

    @staticmethod
    def zpval(z):return -pd.Series(1-st.norm.cdf(pval.zscore(z).values)+ 1e-16,index=z.index).apply(np.log)
    @staticmethod
    def zgenome(x): return -pd.Series(1-st.norm.cdf(pval.zscoreChr(x).values)+ 1e-16,index=x.index).apply(np.log)
    zgenome2tail= lambda x: -pd.Series(1-st.norm.cdf(pval.zscoreChr(x).abs().values)+ 1e-16,index=x.index).apply(np.log)

    @staticmethod
    def gammachi2Test(x,df):return -st.chi2.logsf(x,df), -st.gamma.logsf(x,df/2.,scale=2.),-st.gamma.logsf(x/df,df/2.,scale=2./df)

    @staticmethod
    def empirical(A,Z,positiveStatistic=True):#Z is null scores
        if positiveStatistic:
            a=A[A>0].sort_values()
            z=Z[Z>0].sort_values().values
        else:
            a=A.sort_values()
            z=Z.sort_values().values
        p=np.zeros(a.size)
        j=0
        N=z.size
        for i in range(a.size):
            while j<N:
                if a.iloc[i] <= z[j]:
                    p[i]=N-j +1
                    break
                else:
                    j+=1
            if j==N: p[i]=1
        return -pd.concat([pd.Series(p,index=a.index).sort_index()/(Z.size+1),A[A==0]+1]).sort_index().apply(np.log10)

def DataframetolaTexTable(df, alignment=None, fname=None):
    """
    Args:
        df: pandas dataframe
        alignment: python list of allignment of columns; default is ['c',..]; use ['c', 'p{4in}', 'c', 'c'] for wrapping
        fname: path to save latex table

    Returns:
        object: 
    Returns: latex table
    """
    if alignment is None: alignment = list('c' * (df.shape[1]))
    csv = r'\centering \begin{tabular}{' + '|'.join(alignment) + '}\n' + df.to_csv(sep='\t', index=False).replace('\t',
                                                                                                                  '\t&').replace(
        '\n', '\\\\\n').replace('\\\\\n', '\\\\\hline\n', 1) + r'\end{tabular}'
    csv = csv.replace("%", '\%').replace('inf',r'$\infty$')
    if fname is not None:
        with open(fname, 'w') as f:  print >> f, csv
    return csv

def files(mypath):
    return [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]

def mergeIntervals(data):
    csv = data.to_csv(header=None, sep='\t', index=None)
    csv = \
    Popen(['bedtools', 'merge', '-scores', 'max', '-i'], stdout=PIPE, stdin=PIPE, stderr=STDOUT).communicate(input=csv)[
        0]
    df = pd.DataFrame(map(lambda x: x.split(), csv.split('\n')),
                      columns=['CHROM', 'start', 'end', data.columns[-1]]).dropna()
    df.iloc[:, -1] = df.iloc[:, -1].astype(float)
    df[['start', 'end']] = df[['start', 'end']].astype(int)
    return df


def topK(x, k=2000):
    return x.sort_values(ascending=False).iloc[:k]
def batch(iterable, n=10000000):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]
def getQantilePvalues(X,kde=None,quantiles=np.arange(0,1.001,0.01)):
    if kde is None:   kde=getDensity(X)
    return getPvalKDE(getQuantiles(X, quantiles=quantiles), kde)

# def getDensity(X,width='silverman'):# not width can be string {scott or silverman} or positive real
def getDensity(X,width='scott'):# not width can be string {scott or silverman} or positive real
    from scipy.stats import gaussian_kde
    return gaussian_kde(X, bw_method=width)


def getPvalKDE(x, kde=None):
    if kde is None:   kde=getDensity(x)
    pval=x.apply(lambda y: kde.integrate_box_1d(y,np.inf))
    return -pval.apply(np.log10).sort_index()

def getQuantiles(X,quantiles):
    return X.quantile(quantiles,interpolation='nearest')

@numba.vectorize
def vectorizedLog(x):
    return float(np.log(x))
def numbaLog(df):
    return  pd.DataFrame(vectorizedLog(df.values),columns=df.index,index=df.index).astype(float)

@numba.vectorize
def vectorizedExp(x):
    return float(np.exp(x))
def numbaExp(df):
    return  pd.DataFrame(vectorizedExp(df.values),columns=df.index,index=df.index).astype(float)


def fx(x, s=0.0, h=0.5):
    Z=(1 + s) * x ** 2 + 2 * (1 + h * s) * x * (1 - x) + (1 - x) ** 2
    if Z>0:
        return ((1 + s) * x ** 2 + (1 + h * s) * x * (1 - x)) / (Z)
    else:
        return 0

def sig(x): return 1. / (1 + np.exp(-x))


def logit(p): return np.log(p) - np.log(1 - p)


# def logit_(p): return T.log(p) - T.log(1 - p)


# def sig_(x): return 1. / (1 + T.exp(-x))


def Nu(s, t, nu0, theta, n=2000): return Z(sig(t * s / 2 + logit(nu0)), n, theta)


def forward(x0=0.005,h=0.5,s=1,t=100):
    def f(x,h=0.5,s=1): return ((1+s)*x*x + (1+h*s)*x*(1-x) )/((1+s)*x*x + 2*(1+h*s)*x*(1-x)  +(1-x)**2)
    x=[x0]
    for i in range(t):
        x+=[f(x[-1],h,s)]
    return pd.Series(x)

floatX = 'float64'


def Z(nu, n, theta): return theta * (
nu * ((nu + 1) / 2. - 1. / ((1 - nu) * n + 1)) + (1 - nu) * ((n + 1.) / (2 * n) - 1. / ((1 - nu) * n + 1)))


def roundto(x, base=50000):
    return int(base * np.round(float(x)/base))

def ceilto(x, base=50000):
    return int(base * np.ceil(float(x)/base))

def puComment(fig, comment):
    if comment is not None:
        fig.text(.05, .05, 'Comment: ' + comment, fontsize=26, color='red')


def scanGenome(genome, f=lambda x: x.mean(), uf=None,winSize=50000, step=None, nsteps=5, minSize=None):
    """
    Args:
        genome: scans genome, a series which CHROM and POS are its indices
        windowSize:
        step:
        f: is a SCALAR function or dict of SCALAR fucntions e.g. f= {'Mean' : np.mean, 'Max' : np.max, 'Custom' : np.min}
        Only good for scanning a series with dictionary of scalar fucntions
        uf: is a universal function which returns a dataframe e.g. uf=lambda x: pd.DataFrame(np.random.rand(2,3))
        good for scanning a dataframe (which each column to be scanned) with a scalar or  univesal fucntions
    Returns:
    """
    if len(genome.shape)>1:
        # if genome.shape[1]>1:
        return genome.apply(lambda x: scanGenome(x,f=f,uf=uf,winSize=winSize,step=step,nsteps=nsteps,minSize=minSize))
    if step is None:step=winSize/nsteps
    df = genome.groupby(level='CHROM').apply(lambda ch: scanChromosome(ch.loc[ch.name],f,uf,winSize,step,minSize))
    if minSize is not None:
        df=df[scanGenome(genome,f=len,winSize=winSize,step=step)>=minSize]
    return df

def scanChromosome(x,f,uf,winSize,step,minSize):
    """
    Args:
        chrom: dataframe containing chromosome, positions are index and the index name should be set
        windowSize: winsize
        step: steps in sliding widnow
        f: is a SCALAR function or dict of SCALAR fucntions e.g. f= {'Mean' : np.mean, 'Max' : np.max, 'Custom' : np.min}
        uf: is a universal function which returns a dataframe e.g. uf=lambda x: pd.DataFrame(np.random.rand(2,3))
    Returns:
    """
    print('Chromosome',x.name)
    POS=x.index.get_level_values('POS')
    res=[]
    Bins=np.arange(max(0,roundto(POS.min()-winSize,base=step)), roundto(POS.max(),base=step),winSize)
    for i in range(int(winSize/step)):
        bins=i*step +Bins
        windows=pd.cut( POS, bins,labels=(bins[:-1] + winSize/2).astype(int))
        if uf is None:
            tmp=x.groupby(windows).agg(f)
            tmp.index=tmp.index.astype(int);
            tmp.index.name='POS'

        else:
            tmp=x.groupby(windows).apply(uf)
            tmp=tmp.reset_index()
            tmp.iloc[:,0]=tmp.iloc[:,0].astype(int)
            tmp.columns=['POS']+tmp.columns[1:].tolist()
            tmp= tmp.set_index(tmp.columns[:-1].tolist()).iloc[:,0]
        res+=[tmp]
    df=pd.concat(res).sort_index().dropna()
    # if minSize is not None:
    #     df[df.COUNT < minSize] = None
    #     df = df.loc[:, df.columns != 'COUNT'].dropna()
    return df



def scanGenomeSNP(genome, f=np.mean, winSize=300,skipFromFirst=0,step=None):
    if step is None:
        step=int(winSize/5)
    return  genome.groupby(level=0).apply(lambda x: scanChromosomeSNP(x.iloc[skipFromFirst:],f,winSize,step))

def scanChromosomeSNP(x,f,winSize,step):
    """
    Args:
        chrom: dataframe containing chromosome, positions are index and the index name should be set
        windowSize: winsize
        step: steps in sliding widnow
        f: is a function or dict of fucntions e.g. f= {'Mean' : np.mean, 'Max' : np.max, 'Custom' : np.min}
    Returns:
    """
    BinsStart=pd.Series(np.arange(0, roundto(x.size,base=step),winSize),name='start')
    def createBins(i):
        bins=pd.DataFrame(i*step +BinsStart)
        bins['end'] = bins.start+ winSize
        bins.index=((bins.start+bins.end)/2).astype(int)
        return bins
    bins=pd.concat(map(createBins,range(int(winSize/step)))).sort_index()
    bins[bins>x.size]=None
    bins=bins.dropna().astype(int)
    bins=bins.apply(lambda bin: f(x.iloc[range(bin.start,bin.end)]),axis=1)
    bins.index=x.index[bins.index]
    if bins.shape[0]:return bins.loc[x.name]

def fisher(a):
    return robjects.r(
        'fisher.test(rbind(c({},{}),c({},{})), alternative="less")$p.value'.format(a[0, 0], a[0, 1], a[1, 0], a[1, 1]))[
        0]


def CMH(x, num_rep=3):
    r = robjects.r
    response_robj = robjects.IntVector(x.reshape(-1))
    dim_robj = robjects.IntVector([2, 2, num_rep])
    response_rar = robjects.r['array'](response_robj, dim=dim_robj)
    testres = r['mantelhaen.test'](response_rar);
    pvalue = testres[2][0];
    return pvalue
def CMHcd(cd,DisCoverage=True,eps=1e-20,negLog10=True,damp=1):
    name='CMH '+'-'.join(cd.columns.get_level_values('GEN').unique().values.astype(str))
    a = cd+damp
    num_rep=cd.shape[1]/(2*cd.columns.get_level_values('GEN').unique().size)
    if DisCoverage:
        a.loc[:, pd.IndexSlice[:, :, 'D']] = (a.xs('D', level=2, axis=1) - a.xs('C', level=2, axis=1)).values
    a=a.apply(lambda x: CMH(x.values.reshape(num_rep, 2, 2)), axis=1).rename(name) +eps
    if negLog10: a= -a.apply(np.log10)
    return a

def getPvalFisher(AllGenes, putativeList, myList):
    cont = pd.DataFrame(getContingencyTable(AllGenes=AllGenes, putativeList=putativeList, myList=myList));
    pval = -np.log10(1 - fisher(cont.values))
    return pval, cont


def getContingencyTable(AllGenes, putativeList, myList):
    """
                     |COMALE|Other |
    -------------------------------------
    Putative(Knouwn) |  a   |  b   |  A
    Other            |  c   |  d   |  B
    -------------------------------------
                     |  C   |  D   |  N
    """
    N = AllGenes.size;
    A = putativeList.size;
    C = myList.size
    a = np.intersect1d(putativeList, myList).size
    b = A - a
    c = C - a
    d = (N - C) - b
    df = pd.DataFrame(np.array([[a, b], [c, d]]), index=['Putative', 'Other'], columns=['myList', 'Other'])
    return df



class BED:
    @staticmethod
    def drop_duplicates(file_a,file_b,outfile=None):
        a=pd.concat([pd.read_csv(file_a,sep='\t',header=None),pd.read_csv(file_b,sep='\t',header=None)]).drop_duplicates().sort_values([0,1])
        if outfile is None: return a
        a.to_csv(outfile,index=None,header=None,sep='\t')

    @staticmethod
    def save(dff, fname=None, intervalName=None, intervalScore=None, fhandle=None):
        if 'start' not in dff.columns:
            df = dff.copy(True)
            df['end'] = df.POS
            df.rename(columns={'POS': 'start'}, inplace=True)
        else:
            df = dff
        cols = ['start', 'end']
        if intervalName is not None:
            cols += [intervalName]
            if intervalScore is not None:
                cols += [intervalScore]
        bed = df.reset_index().sort_values(['CHROM', 'start']).set_index('CHROM')[cols]
        if fhandle is not None:
            print >> fhandle, bed.to_csv(header=False, sep='\t')
            fhandle.flush()
        else:
            bed = bed.to_csv(header=False, sep='\t', path_or_buf=fname)
            return bed
    @staticmethod
    def getIntervals(regions, padding=0):
        def get_interval(df, padding, merge=False):
            df = df.sort_index()
            df = pd.DataFrame([df.values, df.index.get_level_values('POS').values - padding,
                               df.index.get_level_values('POS').values + padding], index=['score', 'start', 'end']).T
            df.start = df.start.apply(lambda x: (x, 0)[x < 0])
            df['len'] = df.end - df.start
            return df.set_index('start')
        if len(regions.shape)==1:
            df = regions.groupby(level=0).apply(lambda x: get_interval(x, padding)).reset_index().set_index('CHROM')
        elif 'start' not in regions.columns:
            df = regions.groupby(level=0).apply(lambda x: get_interval(x, padding)).reset_index().set_index('CHROM')
        else:
            df=regions.copy(True)
            df.start-=padding;df.end+=padding
            df.loc[df.start<0,'start']=0
        df=df.reset_index().sort_values(['CHROM','start']).set_index('CHROM')
        df['name'] = range(df.shape[0])
        df.score=(df.score*1000).round()
        df[['start','end','score']]=df[['start','end','score']].applymap(int)
        csv = df[['start', 'end', 'name', 'score']].to_csv(header=False, sep='\t')
        csv = \
        Popen(['bedtools', 'merge', '-scores', 'max', '-i'], stdout=PIPE, stdin=PIPE, stderr=STDOUT).communicate(input=csv)[
            0]
        df = pd.DataFrame(map(lambda x: x.split(), csv.split('\n'))).dropna()
        df.columns=['CHROM', 'start', 'end', 'score']
        df.CHROM=df.CHROM.apply(INT)
        df=df.dropna().set_index('CHROM').astype(int)
        df.score /= 1000
        df['len'] = df.end - df.start
        return df

    @staticmethod
    def intersection(dfa, dfb, dfa_interval_name='Gene_ID',dfb_interval_name='len'):
        dfb.start = dfb.start.astype(int)
        dfb.end = dfb.end.astype(int)
        import tempfile
        with tempfile.NamedTemporaryFile()as f1, tempfile.NamedTemporaryFile() as f2:
            if 'POS' in dfa.index.names:
                BED.save(dfa.reset_index()[['CHROM', 'POS', dfa_interval_name]].drop_duplicates(),intervalName=dfa_interval_name, fhandle=f1)
            else:
                BED.save(dfa, intervalName=dfa_interval_name, fhandle=f1)
            BED.save(dfb, intervalName=dfb_interval_name, fhandle=f2)
            csv = Popen(['bedtools', 'intersect', '-wb', '-wa', '-a', f1.name, '-b', f2.name], stdout=PIPE, stdin=PIPE,
                        stderr=STDOUT).communicate()[0]
        df = pd.DataFrame(map(lambda x: x.split(), csv.split('\n')),).dropna()
        try:
            df=df.iloc[:,[0,1,2,3,7]]
            df.columns=['CHROM', 'start', 'end', dfa_interval_name, dfb_interval_name]
            df=df.set_index('CHROM')
            df[['start', 'end']] = df[['start', 'end']].astype(int)
            return df
        except:
            return None

    @staticmethod
    def saveBEDGraph(a,name,fout_name,color='255,0,0',browser_pos='chrX:1-1000',chromLen=None,URL=None):
        if len(a.shape)==1:
            b=a.dropna().rename('score').reset_index()
        else:
            b=a.dropna().reset_index()
        if 'start' not in b.columns:
            b['start']=b.POS-1;b['end']=b.POS
        if not b.CHROM.astype(str).apply(lambda x:'chr' in x).sum():
            b.CHROM=b.CHROM.apply(lambda x: 'chr{}'.format(x))
        if URL is not None:
            b['UCSC']=None
            for i,row in b.iterrows():
                b.loc[i,'UCSC'] ='=HYPERLINK("{}","GB")'.format(URL.format(row.CHROM,int(row.start),int(row.end)))
            b.to_csv(fout_name+'.tsv',sep='\t',index=None)
        with open(fout_name,'w') as fout:
            print >> fout,"browser position",browser_pos
            print >> fout,"browser hide all"
            if type(b.score.iloc[0])==str :
                print >> fout,'track name="{}" color={} '.format(name,color)
            else:
                print >> fout,'track type=bedGraph name="{}" autoScale=off  visibility=full color={} viewLimits={}:{} priority=20'.format(name,color,min(0,b.score.min()), np.ceil(b.score.max()*10)/10)
            print >>fout,b[['CHROM','start','end','score']].to_csv(sep=' ',index=None,header=None)
            fout.flush()
        import subprocess
        subprocess.call('rm {}'.format(fout_name+'.gz'),shell=True)
        subprocess.call('bgzip {}'.format(fout_name),shell=True)

    @staticmethod
    def saveBEDGraphDF(DF,fout_path,color='255,0,0',colors=None,browser_pos='chrX:10000-12000',ver=None,winSize=None):
        df=DF.reset_index()
        tracks=DF.columns
        colormap={'r':'255,0,0','g':'0,255,0','b':'0,0,255','k':'0,0,0'}
        if winSize is None:
            df['start']=df.POS-1;df['end']=df.POS
            fout_path+='.snps.bedgraph'
        else:
            df['start']=df.POS-winSize/5;df['end']=df.POS+winSize/5-1
            fout_path+='.intervals.bedgraph'
        if ver is not None:
            chromLen=getChromLen(ver)
            df=df.groupby('CHROM').apply(lambda x:x[(x.end<chromLen[x.name]).values & (x.start>=0) ])
        df.set_index(['CHROM','POS']).groupby(level=0).apply(lambda x:x.end)
        df.CHROM=df.CHROM.apply(lambda x: 'chr{}'.format(x))
        with open(fout_path,'w') as fout:
            print >> fout,"browser position",browser_pos
            print >> fout,"browser hide all"
            for i,track in enumerate(tracks):
                if colors is not None:color=colormap[colors[i]]
                print >> fout,'track type=bedGraph name="{}" autoScale=off  visibility=dense color={} viewLimits={}:{} priority=20'.format(track,color,min(0,df[track].min()), np.ceil(df[track].max()*10)/10)
                print >>fout,df[['CHROM','start','end',track]].dropna().to_csv(sep=' ',index=None,header=None)
            fout.flush()
        import subprocess
        subprocess.call('rm -f {}'.format(fout_path+'.gz'),shell=True)
        subprocess.call('bgzip {}'.format(fout_path),shell=True)

    @staticmethod
    def xmap_bed(Interval,hgFrom=19, hgTo=38):
        """
        Args:
            hgFrom: (int) assembly version eg: 19
            hgTo: int) assembly version eg: 38
            interval: dataframe with CHROM, start, end
        Returns:
            out: dataframe with CHROM, start, end
        """
        interval=Interval.copy(True)
        hasChr=False
        if 'chr' in interval.CHROM.iloc[0]: hasChr=True
        if not interval.CHROM.astype(str).apply(lambda x:'chr' in x).sum():
            interval.CHROM='chr'+interval.CHROM.apply(convertToIntStr)
        interval.start=interval.start.astype(int)
        interval.end=interval.end.astype(int)
        chainfile = "/home/arya/CrossMap-0.2.5/chains/hg{}ToHg{}.over.chain.gz".format(hgFrom, hgTo)
        in_file=home+'xmap.in.tmp'
        out_file=home+'xmap.out.tmp'
        import subprocess
        with open(in_file ,'w') as f1:
            BED.save(interval.reset_index()[['CHROM','start','end','index']], fhandle=f1,intervalName='index')
            cmd = "/home/arya/anaconda2/bin/CrossMap.py bed  {} {} {}".format(chainfile, in_file, out_file)
            subprocess.call(cmd,shell=True)
        maped=pd.DataFrame(map(lambda x: x.split(), open(out_file).readlines()),columns=['CHROM','start','end','ID']).dropna()
        maped.ID=maped.ID.astype('int')
        maped=maped.set_index('ID').sort_index()
        maped=pd.concat([interval,maped.loc[:,maped.columns!='CHROM']],1,keys=[hgFrom,hgTo])
        if not hasChr:
            maped[(hgFrom,'CHROM')]=maped[(hgFrom,'CHROM')].apply(lambda x: x[3:])
        maped.sort_values([(hgFrom,'CHROM'),(hgFrom,'start')])
        maped=maped.set_index((hgFrom,'CHROM'))
        maped.index.name='CHROM'
        os.remove(in_file)
        os.remove(out_file)
        os.remove(out_file+'.unmap')
        return maped

def createAnnotation(vcf ,db='BDGP5.75'):
    snps=loadSNPID()
    import subprocess
    fname=vcf.replace('.vcf','.SNPEFF.vcf')
    assert fname!=vcf
    cmd='java -Xmx4g -jar ~/bin/snpEff/snpEff.jar  -ud 0 -s snpeff.html {} {} | cut -f1-8 > {}'.format(db,vcf,fname)
    print(cmd)
    subprocess.call(cmd,shell=True)
    import vcf
    def saveAnnDataframe(x='ANN'):
        print(x)
        ffields = lambda x: x.strip().replace("'", '').replace('"', '').replace(' >', '')
        vcf_reader = vcf.Reader(open(fname, 'r'))
        csv=fname.replace('SNPEFF.vcf',x+'.csv')
        with open(csv,'w') as fout:
            print >>fout,'\t'.join(['CHROM','POS','REF']+map(ffields,vcf_reader.infos[x].desc.split(':')[1].split('|')))
            for rec in  vcf_reader:
                if x in rec.INFO:
                    for line in map(lambda y:('\t'.join(map(str,[rec.CHROM,rec.POS,rec.REF]+y))),map(lambda ann: ann.split('|') ,rec.INFO[x])):
                        # print line
                        print >>fout, line
        uscols=[range(10),range(6)][x=='LOF']
        df=pd.read_csv(csv,sep='\t',usecols=uscols).set_index(['CHROM','POS']).apply(lambda x: x.astype('category'))
        df.join(snps,rsuffix='_flybaseVCF').to_pickle(csv.replace('.csv','.df'))
    saveAnnDataframe('ANN')
    saveAnnDataframe('LOF')


def getEuChromatin(scores):
    def filter(x):
        try:
            return scores.loc[x.name][(x.loc[x.name].start <= scores.loc[x.name].index.values) & (x.loc[x.name].end >= scores.loc[x.name].index.values)]
        except:
            pass
    return getEuChromatinCoordinates().groupby(level=0).apply(filter)

def getEuChromatinCoordinates():
    """http://petrov.stanford.edu/cgi-bin/recombination-rates_updateR5.pl"""
    a = pd.Series(
    """X : 1.22 .. 21.21
    2L : 0.53 .. 18.87
    2R : 1.87 .. 20.86
    3L : 0.75 .. 19.02
    3R : 2.58 .. 27.44 """.split('\n'))
    return pd.DataFrame(a.apply(lambda x: [x.split(':')[0]] + x.split(':')[1].split('..')).tolist(),
                        columns=['CHROM', 'start', 'end']).applymap(str.strip).set_index('CHROM').astype(float) * 1e6

def getChromLen(ver=5):
    return pd.read_csv(home + 'storage/Data/Dmelanogaster/dmel{}.fasta.fai'.format(ver), sep='\t', header=None).replace({'dmel_mitochondrion_genome':'M'}).rename(columns={0:'CHROM'}).set_index(['CHROM'])[1].rename("length")
loadSNPID= lambda : pd.read_csv('/home/arya/storage/Data/Dmelanogaster/dm5.vcf',sep='\t',usecols=range(5),header=None,comment='#',names=['CHROM','POS','ID','REF','ALT']).set_index(['CHROM','POS'])


def smooth(a,winsize): return scan3way(a/a.sum(),winsize,np.mean)
def scan3way(a,winsize,f):
    return pd.concat([a.rolling(window=winsize).apply(f),
                      a.rolling(window=winsize,center=True).apply(f),
                      a.iloc[::-1].rolling(window=winsize).apply(f).iloc[::-1]],
                     axis=1).apply(lambda x: np.mean(x),axis=1)


def localOutliers(a, q=0.99,winSize = 2e6):
    def f(xx):
        window = int(winSize / 10000)
        th =scan3way(xx,window,f=lambda x: pd.Series(x).quantile(q))
        return xx[xx >= th].loc[xx.name]
    return a.loc[a.groupby(level=0).apply(f).index]

def renameColumns(df,suffix,pre=True):
    if pre:
        df.columns=map(lambda x:'{}{}'.format(suffix,x),df.columns)
    else:
        df.columns=map(lambda x:'{}{}'.format(x,suffix),df.columns)
    return df

def makeCategory(df,field):
    try:
        df[field]=df[field].astype('category')
    except:
        pass
    return df

import numbers
def isNumber(x):
    return isinstance(x, numbers.Number)
def convertToIntStr(x):
    if isNumber(x):
        return '{:.0f}'.format(x)
    else:
        return x


def GOEA(bg,study,assoc,alpha=0.05,propagate=False):
    """
    Args:
        bg: list, Background Gene set (possibly the polymorphic genes in the experiment), eg: ['FBgn111','FBgn112',...]
        study: list, a  subset of bg, eg: ['FBgn111']
        assoc: series which index is gene id and value is a set of go terms, eg.
        pd.Series([{GO:0005615, GO:0007566}],index=[ 'FBgn111'])
        alpha:significance level
        propagate: if propagate counts in GO hierarchy
    Returns:

    """
    from goatools.go_enrichment import GOEnrichmentStudy
    from goatools.obo_parser import GODag
    obodag = GODag(dataPath+"go-basic.obo")
    goea= GOEnrichmentStudy(bg,assoc.to_dict(),obodag,propagate_counts = propagate,alpha = alpha, methods = ['fdr_bh'])
    goea_results_all = goea.run_study(study)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < alpha]
    import tempfile
    print(goea_results_sig)
    with tempfile.NamedTemporaryFile()as f:
        goea.wr_tsv(f.name, goea_results_sig)
        df=pd.read_csv(f.name,sep='\t')
    return df

def loadAssociations(species='fly'):
    taxid={'fly':7227, 'human':9606,'mouse':10090,'rat':10116}
    from goatools.associations import read_ncbi_gene2go
    aa=pd.Series(read_ncbi_gene2go(dataPath+"gene2go", taxids=[taxid[species]]))
    if species == 'fly':
        bb=pd.read_pickle(dataPath+'fruitfly.mygene.df')
        bb.index=map(int,bb.index)
        aa=bb.join(aa.rename('GO')).set_index("FLYBASE")['GO']
    return aa

def getGeneName(geneIDs,species='fruitfly'):
    try:
        return pd.read_pickle(dataPath+'{}.mygene.df'.format(species))
    except:
        import mygene
        names=mygene.MyGeneInfo().querymany(geneIDs, scopes="entrezgene,flybase",  species=species, as_dataframe=True,fields='all')
        names.to_pickle(dataPath+'{}.mygene.df'.format(species))
        return names

class Dmel:
    @staticmethod
    def geneAssociations(relation=False):
        """
        Args:
            relation: if True returns 1 to N relation otherwise returns 1 to a set relation
        Returns:
        """
        # ftp://ftp.flybase.net/releases/current/precomputed_files/go/gene_association.fb.gz
        a=pd.read_csv(dataPathDmel+'gene_association.fb',sep='\t',comment='!',header=None).set_index(1)[4].rename('GO')
        a.index.name=None
        if relation:
            return a
        else:
            return a.groupby(level=0).apply(lambda x: set(x.tolist()))

    @staticmethod
    def geneCoordinates(assembly=5,allOrganisms=False):
        # ftp://ftp.flybase.net/releases/FB2014_03/precomputed_files/genes/gene_map_table_fb_2014_03.tsv.gz
        # ftp://ftp.flybase.net/releases/FB2016_04/precomputed_files/genes/gene_map_table_fb_2016_04.tsv.gz
        fname=dataPathDmel+('gene_map_table_fb_2016_04.tsv','gene_map_table_fb_2014_03.tsv')[assembly==5]
        import re
        try:
            return pd.read_pickle(fname+'.{}.df'.format(('dmel','all')[allOrganisms]))
        except:
            col=4;# sequence_loc column
            a=pd.read_csv(fname, sep='\t', comment='#',header=None,na_values='',keep_default_na=False).dropna(subset=[col])
            df=pd.DataFrame([re.split('\..|:|\(',r[col])[:-1]+r[:2].tolist()+['\\' not in r[0]] for _,r in a.iterrows()])
            df.columns=['CHROM','start','end','name','FBgn','dmel']
            df[['start', 'end']] = df[['start', 'end']].astype(int)
            df = df.sort_values(['CHROM', 'start'])
            df.index=range(df.shape[0])
            df.to_pickle(fname+'.all.df')
            df=df[df.dmel]
            df.index=range(df.shape[0])
            df.to_pickle(fname+'.dmel.df')
            return df[df.dmel]

def getFDR(a,T=[0.05,0.01,0.0025,0.001,0.0001]):
    b=pd.DataFrame([(t, a.size*t,a[a>=-np.log10(t)].size) for t in T],columns=['t','mt','discoveries']);
    b['fdr']=b.mt/b.discoveries
    return b

def to_hdf5(filename, df, metadf=None, **kwargs):
    store = pd.HDFStore(filename)
    store.put('data', df)
    if metadf is not None:
        store.put('meta', metadf)
    store.get_storer('data').attrs.metadata = kwargs
    store.close()

def read_hdf5(filename):
    with pd.HDFStore(filename) as store:
        data = store['data']
        metadf=None
        try:
            metadf = store['meta']
        except:
            pass
        metadata = store.get_storer('data').attrs.metadata
    return data, metadf, metadata

def loadGowinda(path='/home/arya/out/real/gowinda/',fname='cands.final.out.snp.tsv'):
    gowinda = pd.read_csv(path+fname, sep='\t', header=None)[[0, 4, 5, 6, 7, 8, 9]]
    gowinda.columns = ['GO ID', '-log($p$-value)', 'Hits', 'Num of Genes', 'Total Genes', 'GO Term', 'Genes']
    return gowinda
def saveGowinda(cands,all,path=outpath + 'real/gowinda/', fname='cands.final.txt'):
    cands.sort_index().reset_index().drop_duplicates().dropna().to_csv(path+fname,sep='\t', header=None, index=False)
    all.sort_index().reset_index().drop_duplicates().dropna().to_csv(path+'allsnps.txt',sep='\t', header=None, index=False)

def get_gene_coordinates(gene_name, hg="hg19"):
    """ Return the coordinates of genes
    :param gene_name: Gene Name
    :param hg: assembly version, must be in {hg16, hg17, hg18, hg19, hg38}
    :return:  [Gene Name, Chromosome, Start Position, End Position]
    """
    """ Return the coordinates of gene "g" in hg19/CRCH37
    :param g: Gene Name
    :return:  [Gene Name, Chromosome, Start Position, End Position]
    """
    if hg not in ["hg16", "hg17", "hg18", "hg19", "hg38"]:
        print("'%s' is not a valid assembly name! it must be in {hg16, hg17, hg18, hg19, hg38}."%hg)
        return None
    gene_name = gene_name.upper()
    GeneCordinateTable="/media/alek/DATA/DB/GeneCordinatesTable/%s.GeneCordinatesTable.csv"%hg
    df=pd.read_csv(GeneCordinateTable, sep='\t').drop_duplicates("name2")
    return list(df.loc[df["name2"] == gene_name][["name2", "chrom", "cdsStart", "cdsEnd"]].values.squeeze())

def get_overlapping_genes(chrom,region_st,region_end, hg="hg19"):
    """ Return the coordinates of genes
    :param gene_name: Gene Name
    :param hg: assembly version, must be in {hg16, hg17, hg18, hg19, hg38}
    :return:  [Gene Name, Chromosome, Start Position, End Position]
    """
    """ Return the coordinates of gene "g" in hg19/CRCH37
    :param g: Gene Name
    :return:  [Gene Name, Chromosome, Start Position, End Position]
    """
    if hg not in ["hg16", "hg17", "hg18", "hg19", "hg38"]:
        print("'%s' is not a valid assembly name! it must be in {hg16, hg17, hg18, hg19, hg38}."%hg)
        return None
    hg='hg19'

    GeneCordinateTable=dataPath+"/Human/GeneCordinatesTable/%s.GeneCordinatesTable.csv"%hg
    df=pd.read_csv(GeneCordinateTable, sep='\t').drop_duplicates("name2")
    c1 = df["chrom"] == "chr%s" % chrom
    c2 = ~ ((df["cdsStart"]>region_end)|(df["cdsEnd"]<region_st))
    return df.loc[c1&c2]['name2'].tolist()

class SynchronizedFile:
    @staticmethod
    def processSyncFileLine(x,dialellic=True):
        z = x.apply(lambda xx: pd.Series(xx.split(':'), index=['A', 'T', 'C', 'G', 'N', 'del'])).astype(float).iloc[:, :4]
        ref = x.name[-1]
        alt = z.sum().sort_values()[-2:]
        alt = alt[(alt.index != ref)].index[0]
        if dialellic:   ## Alternate allele is everthing except reference
            return pd.concat([z[ref].astype(int).rename('C'), (z.sum(1)).rename('D')], axis=1).stack()
        else:           ## Alternate allele is the allele with the most reads
            return pd.concat([z[ref].astype(int).rename('C'), (z[ref] + z[alt]).rename('D')], axis=1).stack()

    @staticmethod
    def load(fname = './sample_data/popoolation2/F37.sync'):
        print('loading',fname)
        cols=pd.read_csv(fname+'.pops', sep='\t', header=None, comment='#').iloc[0].apply(lambda x: map(int,x.split(','))).tolist()
        data=pd.read_csv(fname, sep='\t', header=None).set_index(range(3))
        data.columns=pd.MultiIndex.from_tuples(cols)
        data.index.names= ['CHROM', 'POS', 'REF']
        data=data.sort_index().reorder_levels([1,0],axis=1).sort_index(axis=1)
        data=data.apply(SynchronizedFile.processSyncFileLine,axis=1)
        data.columns.names=['REP','GEN','READ']
        data=SynchronizedFile.changeCtoAlternateAndDampZeroReads(data)
        data.index=data.index.droplevel('REF')
        return data

    @staticmethod
    def changeCtoAlternateAndDampZeroReads(a):
        C = a.xs('C', level=2, axis=1).sort_index().sort_index(axis=1)
        D = a.xs('D', level=2, axis=1).sort_index().sort_index(axis=1)
        C = D - C
        if (D == 0).sum().sum():
            C[D == 0] += 1
            D[D == 0] += 2
        C.columns = pd.MultiIndex.from_tuples([x + ('C',) for x in C.columns], names=C.columns.names + ['READ'])
        D.columns = pd.MultiIndex.from_tuples([x + ('D',) for x in D.columns], names=D.columns.names + ['READ'])
        return pd.concat([C, D], axis=1).sort_index(axis=1).sort_index()

def INT(x):
    try: return int(x)
    except: return x

def intIndex(df):
    names=df.index.names
    df=df.reset_index()
    df[names]=df[names].applymap(INT)
    return df.set_index(names).sort_index()

def uniqIndex(df,keep=False,subset=['CHROM','POS']): #keep can be first,last,None
    names=df.index.names
    if subset is None: subset=names
    return df.reset_index().drop_duplicates(subset=subset,keep=keep).set_index(names).sort_index()


def mask(genome,CHROM=None,start=None,end=None,interval=None,pad=0,returnIndex=False):
    if interval is not None: CHROM, start, end = interval.CHROM, interval.start, interval.end
    start-=pad;end+=pad
    if returnIndex:
        return (genome.index.get_level_values('CHROM')==CHROM) & (genome.index.get_level_values('POS')>=start)&(genome.index.get_level_values('POS')<=end)
    else:
        tmp=genome.loc[CHROM]
        return tmp[(tmp.index.get_level_values('POS')>=start)&(tmp.index.get_level_values('POS')<=end)]


def getRegionPrameter(CHROM,start,end):
    if start is not None and end is not None:CHROM='{}:{}-{}'.format(CHROM,start,end)
    elif start is None and end is not None:CHROM='{}:-{}'.format(CHROM,end)
    elif start is not None and end is None :CHROM='{}:{}-'.format(CHROM,start)
    return CHROM


class VCF:
    @staticmethod
    def getField(fname,field='POS'):
        fields={'CHROM':1,'POS':2,'ID':3}
        cmd="zgrep -v '#' {} | cut -f{}".format(fname,fields[field])
        return pd.Series(Popen([cmd], stdout=PIPE, stdin=PIPE, stderr=STDOUT,shell=True).communicate()[0].strip().split('\n')).astype(int)

    @staticmethod
    def header(fname):
        cmd="zgrep -w '^#CHROM' -m1 {}".format(fname)
        return Popen([cmd], stdout=PIPE, stdin=PIPE, stderr=STDOUT,shell=True).communicate()[0].split('\n')[0].split()
    @staticmethod
    def headerSamples(fname):
        return map(INT,VCF.header(fname)[9:])

    @staticmethod
    def loadPanel(fname=dataPath1000GP+'integrated_call_samples_v3.20130502.ALL.panel'):
        return  pd.read_table(fname,sep='\t').dropna(axis=1)

    @staticmethod
    def getDataframeColumns(fin,panel=None):
        cols=[]
        if panel is not None:
            panel=VCF.loadPanel(panel).set_index('sample')[['super_pop','pop']]
            f=lambda x: tuple(panel.loc[x].tolist())
            for x in VCF.headerSamples(fin):cols+=[f(x)+( x,'A'),f(x)+(x,'B')]
            cols=pd.MultiIndex.from_tuples(cols,names=['SPOP', 'POP', 'ID','HAP'])
        else:
            for x in VCF.headerSamples(fin):cols+=[( x,'A'),(x,'B')]
            cols=pd.MultiIndex.from_tuples(cols,names=[ 'ID','HAP'])
        return cols

    @staticmethod
    def getDataframe(CHROM,start=None,end=None,
                     fin=dataPath1000GP+'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
                     bcftools="/home/arya/bin/bcftools/bcftools",
                     panel=dataPath1000GP+'integrated_call_samples_v3.20130502.ALL.panel'):
        reg=getRegionPrameter(CHROM,start,end)
        fin=fin.format(CHROM)
        cmd="{} filter {} -i \"N_ALT=1 & TYPE='snp'\" -r {} | {} annotate -x INFO,FORMAT,FILTER,QUAL,FORMAT | grep -v '#' | tr '|' '\\t'|  tr '/' '\\t' | cut -f1-5,10-".format(bcftools,fin,reg,bcftools)
        csv=Popen([cmd], stdout=PIPE, stdin=PIPE, stderr=STDOUT,shell=True).communicate()[0].split('\n')
        df = pd.DataFrame(map(lambda x: x.split('\t'),csv)).dropna().set_index(range(5)).astype(int)
        df.index.names=['CHROM','POS', 'ID', 'REF', 'ALT']
        df.columns=VCF.getDataframeColumns(fin,panel)
        return df

    @staticmethod
    def computeFreqs(CHROM,start=None,end=None,
                     fin=dataPath1000GP+'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
                     panel=dataPath1000GP+'integrated_call_samples_v3.20130502.ALL.panel',hap=False):
        try:
            a=VCF.getDataframe(CHROM,start,end,fin=fin,panel=panel)
            if hap: return a
            if panel is not None:
                return pd.concat([a.mean(1).rename('ALL'),a.groupby(level=0,axis=1).mean(),a.groupby(level=1,axis=1).mean()],1)
            else:
                return a.mean(1).rename('ALL')
        except:
            return None
    @staticmethod
    def computeFreqsChromosome(CHROM,
                       fin=dataPath1000GP+'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz',
                       panel=dataPath1000GP+'integrated_call_samples_v3.20130502.ALL.panel',winSize=100000,hap=False):
        CHROM=INT(CHROM)
        #if CHROM == 'Y': fin=dataPath1000GP+'ALL.chr{}.phase3_integrated_v1b.20130502.genotypes.vcf.gz'
        #if CHROM == 'X': fin=dataPath1000GP+'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'
        import vcf
        try:
            L=vcf.Reader(open(fin.format(CHROM), 'r')).contigs[str(CHROM)].length
        except:
            L=vcf.Reader(open(fin.format(CHROM), 'r')).contigs['chr{}'.format(CHROM)].length
        print(CHROM,int(L/1e6),'Mbp')
        a=[VCF.computeFreqs(CHROM,start,end=start+winSize-1,fin=fin,panel=panel,hap=hap) for start in xrange(0,ceilto(L,winSize),winSize)]
        print(a)
        return intIndex(uniqIndex(pd.concat([x  for x in a if x is not None]),subset=['CHROM','POS']))
    @staticmethod
    def createGeneticMap(VCFin, chrom,gmpath=dataPath+'Human/map/GRCh37/plink.chr{}.GRCh37.map'):
        gm = pd.read_csv(gmpath.format(chrom), sep='\t', header=None,names=['CHROM','ID','GMAP','POS'])
        df = pd.DataFrame(VCF.getField(VCFin).rename('POS'))
        df['GMAP'] = np.interp(df['POS'].tolist(), gm['POS'].tolist(),gm['GMAP'].tolist())
        df['CHROM']=chrom
        df['ID']='.'
        df[['CHROM','ID','GMAP','POS']].to_csv(VCFin+'.map',sep='\t',header=None,index=None)

    @staticmethod
    def subset(VCFin, pop,panel,chrom):
        bcf='/home/arya/bin/bcftools/bcftools'
        assert len(pop)
        if pop=='ALL' or pop is None:return VCFin,None
        print('Creating a vcf.gz file for individuals of {} population'.format(pop))
        fileSamples='{}.{}.chr{}'.format(panel,pop,chrom)
        fileVCF=VCFin.replace('.vcf.gz','.{}.vcf.gz'.format(pop))
        os.system('grep {} {} | cut -f1 >{}'.format(pop,panel,fileSamples))
        cmd="{} view -S {} {} | {} filter -i \"N_ALT=1 & TYPE='snp'\" -O z -o {}".format(bcf,fileSamples,VCFin,bcf,fileVCF)
        os.system(cmd)
        return fileVCF

    @staticmethod
    def loadDP(fname):
        a= pd.read_csv(fname,sep='\t',na_values='.').set_index(['CHROM','POS'])
        a.columns=pd.MultiIndex.from_tuples(map(lambda x:(int(x.split('R')[1].split('F')[0]),int(x.split('F')[1])),a.columns))
        return  a

    @staticmethod
    def loadCD(vcfgz,vcftools='~/bin/vcftools_0.1.13/bin/vcftools'):
        """
            vcfgz: vcf file where samples are in the format of RXXFXXX
        """
        vcf=os.path.basename(vcfgz)
        path=vcfgz.split(vcf)[0]
        os.system('cd {0} && {1} --gzvcf {2} --extract-FORMAT-info DP && {1} --gzvcf {2} --extract-FORMAT-info AD'.format(path,vcftools,vcf))
        fname='out.{}.FORMAT'
        a=map(lambda x: VCF.loadDP(path +fname.format(x)) ,['AD','DP'])
        a=pd.concat(a,keys=['C','D'],axis=1).reorder_levels([1,2,0],1).sort_index(1)
        a.columns.names=['REP','GEN','READ']
        return a

def polymorphic(data, minAF=0.01,mincoverage=10,index=True):
    def poly(x):return (x>=minAF)&(x<=1-minAF)
    C,D=data.xs('C',level='READ',axis=1),data.xs('D',level='READ',axis=1)
    I=(C.sum(1)/D.sum(1)).apply(lambda x:poly(x)) & ((D>=mincoverage).mean(1)==1)
    if index:
        return I
    return data[I]
def loadGenes(Intervals=True):
    a=pd.read_csv(dataPath+'Human/WNG_1000GP_Phase3/gene_info.csv')[['chrom','pop','gene','POS_hg19']].rename(columns={'chrom':'CHROM','POS_hg19':'POS'})
    a.CHROM=a.CHROM.apply(lambda x: INT(x[3:]))
    a=a.set_index('pop')
    if Intervals:
        a['start']=a.POS-2e6
        a['end']=a.POS+2e6
        a['name']=a.gene
    return a

def normalizeIHS(a,field=None):
    if field is None:
        field=a.columns[-1]
    m=a.set_index('x')[field].groupby(level=0).mean().sort_index()
    s=a.set_index('x')[field].groupby(level=0).std().sort_index()
    return (a[field]-m.loc[a.x].values)/s.loc[a.x].values


def filterGap(a,assempbly=38,pad=50000):
    gap=pd.read_csv(dataPath+'Human/hg{}.gap'.format(assempbly),sep='\t')[['chrom','chromStart','chromEnd']].rename(columns={'chrom':'CHROM','chromStart':'start','chromEnd':'end'}).reset_index()
    gap.start-=pad;gap.end+=pad;gap.start[gap.start<0]=0
    gap.CHROM=gap.CHROM.apply(lambda x: x[3:])
    gap=BED.intersection(a,gap,dfa_interval_name=a.name,dfb_interval_name='index').rename(columns={'start':'POS'})
    gap.index=map(INT,gap.index);gap.index.name='CHROM'
    gap=gap.set_index('POS',append=True)[a.name].sort_index()
    a.loc[gap.index]=None;a=a.dropna()
    return a

def get_gap(assempbly=38,pad=50000):
    gap=pd.read_csv(dataPath+'Human/hg{}.gap'.format(assempbly),sep='\t')[['chrom','chromStart','chromEnd']].rename(columns={'chrom':'CHROM','chromStart':'start','chromEnd':'end'}).reset_index()
    gap.start-=pad;gap.end+=pad;gap.start[gap.start<0]=0
    gap.CHROM=gap.CHROM.apply(lambda x: x[3:])
    return gap
def loadFst(fname):
    a=pd.read_csv(fname,sep='\t');
    a.CHROM=a.CHROM.replace({'22':22});a=a.set_index(['CHROM','POS']).sort_index().iloc[:,0].rename('Fst')
    a=a.loc[range(1,23)]
    return a[a>0]
