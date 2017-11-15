'''
Copyleft May 11, 2016 Arya Iranmehr, PhD Student, Bafna Lab, UC San Diego,  Email: airanmehr@gmail.com
'''
import matplotlib as mpl
import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns

import Utils.Util as utl


def setStyle(style="darkgrid", lw=2, fontscale=1, fontsize=10):
    sns.axes_style(style)
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': fontsize});
    mpl.rc('text', usetex=True)
    sns.set_context(font_scale=fontscale, rc={"lines.linewidth": lw})


class PLOS:
    max_width = 7.5
    min_width = 2.6
    max_height = 8.75
    dpi = 300
    extention = 'tiff'

    @staticmethod
    def get_figsize(width=None, height=None):
        if width is not None:
            width = min(width, PLOS.max_width)
            return (width, 2. / 3 * width)
        else:
            return (6, 4)


def get_axis_limits(ax, upper=True):
    return ax.get_xlim()[(0, 1)[upper]], ax.get_ylim()[(0, 1)[upper]]


def annotate(comment, loc=1, fontsize=26, xpad=0.05, ypad=0.05, ax=None, axtoplot=None):
    """
    Args:
        comment: text
    """
    if ax is None: ax = plt.gca()
    if axtoplot is None: axtoplot = ax
    xrang = getAxRange(ax, 0)
    yrang = getAxRange(ax, 1)
    xy = get_axis_limits(ax, upper=False)[0] + xpad * xrang, get_axis_limits(ax)[1] - ypad * yrang
    axtoplot.annotate(comment, xy=xy, xycoords='data', size=fontsize, horizontalalignment='left',
                       verticalalignment='top')


def getAxRange(ax, axi=0):
    return get_axis_limits(ax, upper=True)[axi] - get_axis_limits(ax, upper=False)[axi]

def getColorMap(n):
    colors = ['darkblue', 'r', 'darkviolet', 'green', 'k', 'darkorange', 'olive', 'darkgrey', 'chocolate', 'rosybrown',
              'gold', 'aqua']
    if n == 1: return [colors[0]]
    if n <= len(colors):
        return colors[:n]
    return [mpl.cm.jet(1. * i / n) for i in range(n)]


def getMarker(n, addDashed=True):
    markers = np.array(['o', '^', 's', 'v', '<', '>', 'D', 'd', 'h', '*', 'p', '3', '2', '4', 'H', '8'])[:n]
    if addDashed: markers = map(lambda x: '--' + x, markers)
    return markers

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':30}) ;
mpl.rc('text', usetex=True)

def addGlobalPOSIndex(df,chroms):
    if df is not None:
        df['gpos'] = df.POS + chroms.offset.loc[df.CHROM].values
        df.set_index('gpos', inplace=True);
        df.sort_index(inplace=True)


def GenomeChromosomewise(df, candSNPs=None, genes=None, axes=None,outliers=None):
    markerSize = 6
    fontsize = 6
    chrsize = df.reset_index().groupby('CHROM').POS.max()
    if axes is None:
        if chrsize.shape[0]>1:
            _, axes = plt.subplots(int(np.ceil(chrsize.shape[0] / 2.)), 2, sharey=True, dpi=200, figsize=(10, 6));
            ax = axes.reshape(-1)
        else:
            ax = [plt.subplots(1,1, sharey=True, dpi=200, figsize=(10, 6))[1]]

    for j, (chrom, a) in enumerate(df.groupby(level=0)):
        if candSNPs is not None:
            try:
                candSNPs.loc[chrom]
                for pos in candSNPs.loc[chrom].index.values:
                    ax[j].axvline(pos, color='r', linewidth=0.5, alpha=0.5)
                    ax[j].annotate(
                        '{:.0f},{:.2f}'.format(candSNPs['rank'].loc[(chrom, pos)], candSNPs.nu0.loc[(chrom, pos)]),
                        xy=(pos, a.max()), xytext=(pos, a.max()), fontsize=fontsize - 2)

            except:
                pass

        if genes is not None:
            try:
                X = genes.loc[chrom]
                if len(genes.loc[chrom].shape) == 1:
                    X = pd.DataFrame(X).T
                for _, row in X.iterrows():
                    ax[j].fill_between([row.start, row.end], a.min(), a.max(), color='r')
                    ax[j].annotate(row['name'], xy=(row.start, a.max()), xytext=(row.start, a.max()),
                                   fontsize=fontsize - 2)


            except:
                pass

        ax[j].scatter(a.loc[chrom].index, a.loc[chrom], s=markerSize, alpha=0.8, edgecolors='none')

        if outliers is not None:
            try:
                ax[j].scatter(outliers.loc[chrom].index, outliers.loc[chrom], s=markerSize, c='r', alpha=0.8, edgecolors='none')
            except:
                pass

        setSize(ax[j], fontsize)
        ax[j].set_xlim([-1000, chrsize[chrom] + 1000])
        # ax[j].set_title(chrom, fontsize=fontsize+2)
        annotate(chrom, ax=ax[j],fontsize=fontsize+4)
        ax[j].locator_params(axis='x', nbins=10)
    plt.tight_layout(pad=0.1)
    plt.gcf().subplots_adjust(bottom=0.1)


def Manhattan(data, columns=None, names=None, fname=None, colors=['black', 'gray'], markerSize=20, ylim=None, show=True,
              std_th=None, top_k=None, cutoff=None, common=None, Outliers=None, shade=None, fig=None, ticksize=12,
              sortedAlready=False,lw=1,axes=None):
    def reset_index(x):
        if x is None: return None
        if 'CHROM' not in x.columns.values:
            return x.reset_index()
        else:
            return x
    if type(data) == pd.Series:
        DF = pd.DataFrame(data)
    else:
        DF = data

    if columns is None: columns=DF.columns
    if names is None:names=columns

    df = reset_index(DF)
    Outliers = reset_index(Outliers)
    if not sortedAlready: df = df.sort_index()
    if not show:
        plt.ioff()
    from itertools import cycle
    def plotOne(b, d, name, chroms,common,shade,ax):
        a = b.dropna()
        c = d.loc[a.index]
        if ax is None:
            ax=plt.gca()
        if shade is not None:
            for _ ,  row in shade.iterrows():
                ax.fill_between([row.gstart, row.gend], a.min(), a.max(), color='b', alpha=0.4)

        ax.scatter(a.index, a, s=markerSize, c=c, alpha=0.8, edgecolors='none')

        outliers=None
        if Outliers is not None:
            outliers=Outliers[name].dropna()
        if cutoff is not None:
            outliers = a[a >= cutoff[name]]
        elif top_k is not None:
            outliers = a.sort_values(ascending=False).iloc[:top_k]
        elif std_th is not None:
            outliers = a[a > a.mean() + std_th * a.std()]
        if outliers is not None:
            if len(outliers):
                ax.scatter(outliers.index, outliers, s=markerSize, c='r', alpha=0.8, edgecolors='none')
                ax.axhline(outliers.min(), color='k', ls='--',lw=lw)


        if common is not None:
            for ii in common.index: plt.axvline(ii,c='g',alpha=0.5)

        ax.axis('tight');
        ax.set_xlim(max(0,a.index[0]-10000), a.index[-1]);
        setSize(ax,ticksize)
        ax.set_ylabel(name, fontsize=ticksize * 1.5)
        if chroms.shape[0]>1:
            plt.xticks([x for x in chroms.mid], [str(x) for x in chroms.index], rotation=-90, fontsize=ticksize * 1.5)
        plt.setp(plt.gca().get_xticklabels(), visible=False)
        plt.locator_params(axis='y', nbins=4)
        mpl.rc('ytick', labelsize=ticksize)
        if ylim is not None:    plt.ylim(ymin=ylim)
    chroms = pd.DataFrame(df.groupby('CHROM').POS.apply(lambda x:x.max()-x.min()).rename('len').loc[df.reset_index().CHROM.unique()] + 1000)
    chroms = pd.DataFrame(df.groupby('CHROM').POS.apply(lambda x:x.max()).rename('len').loc[df.reset_index().CHROM.unique()] + 1000)
    chroms['offset'] = np.append([0], chroms.len.cumsum().iloc[:-1].values)
    chroms['color'] = [c for (_, c) in zip(range(chroms.shape[0]), cycle(colors))]
    chroms['start']=df.groupby('CHROM').POS.min()
    chroms['mid'] = [x + y / 2 for x, y in zip(chroms.offset+chroms.start, chroms.len)]
    chroms['mid'] = [x + y / 2 for x, y in zip(chroms.offset+chroms.start, chroms.len)]
    df['color'] = chroms.color.loc[df.CHROM].values
    df['gpos'] = df.POS + chroms.offset.loc[df.CHROM].values
    df['color'] = chroms.color.loc[df.CHROM].values
    df.set_index('gpos', inplace=True);

    if shade is not None:
        shade['gstart']=shade.start #
        shade['gend']=shade.end #
        if chroms.shape[0]>1:
            shade['gstart']+= chroms.offset.loc[shade.CHROM].values
            shade['gend']+=+ chroms.offset.loc[shade.CHROM].values
        if 'name' in shade.columns:
            shade.sort_values('gstart',ascending=False,inplace=True)
            shade['ID']=range(1,shade.shape[0]+1)
    addGlobalPOSIndex(common, chroms);
    addGlobalPOSIndex(Outliers, chroms)
    if fig is None and axes is None:
        fig,axes=plt.subplots(columns.size, 1, sharex=True,figsize=(20, columns.size * 4));
        if columns.size==1:
            axes=[axes]
    elif axes is None:
        axes=fig.axes
    for i in range(columns.size):
        plotOne(df[columns[i]], df.color, names[i], chroms,common, shade,axes[i])
    plt.setp(plt.gca().get_xticklabels(), visible=True)
    xlabel='Chromosome'
    if chroms.shape[0]==1:xlabel+=' {}'.format(chroms.index[0])
    plt.xlabel(xlabel, size=ticksize * 1.5)
    plt.gcf().subplots_adjust(bottom=0.2)
    if fname is not None:
        print('saving ', fname)
        plt.savefig(fname)
    if not show:
        plt.ion()

    return fig




def TimeSeries(data, methodColumn=None, ax=None, fname=None, color='r', ci=1,shade=[0,50],samplingTimes=None):
    """
    Args:
        data: a dataframe containing mu and st fields,
        methodColumn: when method column is given, it plots together
        ax:
        fname:
    Returns:
    """
    if ax is None: fig=plt.figure(figsize=(12,4), dpi=200)
    if methodColumn is None:
        dfs=[('aa',data)]
    else:
        dfs=data.groupby(methodColumn)
    for name,df in dfs:
        if 'color' in df.columns:color=df.color.unique()[0]
        df.mu.plot(linewidth=1, color=color, label=name, ax=ax)
        # plt.gca().fill_between(df.index,  (df.mu+df.st).apply(lambda x:min(x,1)), (df.mu-df.st).apply(lambda x:max(x,0)), color=color, alpha=0.25)
        ax.fill_between(df.index.values.astype(int), (df.mu + ci * df.st), (df.mu - ci * df.st), color=color,
                        alpha=0.25)
    if shade is not None:
        ax.axvspan(shade[0], shade[1], alpha=0.25, color='black')
        ax.set_xticks(np.append([50], plt.xticks()[0]))
    if samplingTimes is  not None:
        for t in samplingTimes:ax.axvline(t,color='k',ls='-',lw=0.5,alpha=0.5)

    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':26}) ;
    mpl.rc('text', usetex=True)
    if fname is not None:
        plt.savefig(fname)



def QQPval(a,z,nq=20, s=40, alpha=0.8, fname=None):
    """pplt.QQPval(exp10(logpa),exp10(logpz))"""
    def getQuantilesLog2():
        q=[1]
        for i in range(nq):q+=[q[-1]/2.]
        q=pd.Series(q,index=q).iloc[1:]
        return q
    q=getQuantilesLog2()
    qq=pd.DataFrame(q.apply(lambda x: [abs((x)),z.quantile(x),a.quantile(x)]).sort_index().tolist(),index=q,columns=['expected','null','data']).applymap(lambda x: -np.log10(x))
    plt.figure(figsize=(8,6),dpi=200)
    qq.plot.scatter(x='expected',y='null',color='k',s=s,alpha=alpha,ax=plt.gca())
    qq.plot.scatter(x='expected',y='data',ax=plt.gca(),s=s,alpha=alpha,color='r',lw = 0);
    plt.ylim([-1, plt.ylim()[1]]);
    xmax = plt.xlim()[1]
    plt.plot([0, xmax], [0, xmax],ls='--', c="k",alpha=0.3)
    plt.xlim([0,xmax])
    plt.xlabel('Expected -log$_{10}$($p$-value)');
    plt.ylabel('Observed -log$_{10}$($p$-value)')
    mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size':26}) ;
    mpl.rc('text', usetex=True)
    if fname is not None: plt.savefig(fname)


def plotSiteReal(site, ax=None, fontsize=8, legend=False, title=None):
    if ax is None:
        dpi = 300
        _, ax = plt.subplots(1, 1, figsize=(3, 2), dpi=dpi, sharex=True, sharey=True)
        sns.set_context("notebook", font_scale=1, rc={"lines.linewidth": 1.2})
    pos = site.name
    site = site.sort_index().groupby(level=[0, 1]).apply(lambda x: (x.iloc[0], x.iloc[1]))
    df = site.apply(lambda x: pd.Series(np.random.binomial(x[1], x[0] / x[1], 10000)) / x[1]).T
    df = df.stack(['REP', 'GEN']).reset_index(['REP', 'GEN'])
    idx = pd.Series(range(site.index.get_level_values('GEN').unique().shape[0]), index=np.sort(site.index.get_level_values('GEN').unique()))
    ax = sns.boxplot(data=df, x='GEN', y=0, hue='REP', width=0.3, ax=ax);
    for i, mybox in enumerate(ax.artists):
        # Change the appearance of that box
        c = mybox.get_facecolor()
        mybox.set_facecolor('None')
        mybox.set_edgecolor(c)
        for j in range(i * 6, i * 6 + 6):
            line = ax.lines[j]
            line.set_color(c)
            line.set_mfc('None')
    for nnn,(_, xxx) in enumerate(site.apply(lambda x: x[0] / x[1]).unstack('REP').iteritems()):
        # print idx.loc[xxx.dropna().index] + (nnn - 1) * (0.1)
        try:
            pd.Series(xxx.dropna().values, index=idx.loc[xxx.dropna().index] + (nnn - 1) * (0.1)).plot(style='-o',
                                                                                                       color=
                                                                                                       sns.color_palette()[
                                                                                                           nnn], ax=ax,
                                                                                                       markersize=3,
                                                                                                       grid=False,
                                                                                                       linewidth=0.5)
        except:
            pass
    handles, labels = ax.get_legend_handles_labels();
    ax.set_xlim([ax.get_xlim()[0] - ax.get_xlim()[1] * 0.03, ax.get_xlim()[1] + ax.get_xlim()[1] * 0.03])
    ax.set_ylim([ax.get_ylim()[0] - ax.get_ylim()[1] * 0.03, ax.get_ylim()[1] + ax.get_ylim()[1] * 0.03])
    if legend:
        ax.legend(handles[3:], map(lambda x: 'Replicate {}'.format(int(x) + 1), labels[3:]), loc='best', title='',
                  fontsize=fontsize - 2)
    else:
        ax.legend_.remove()
    ax.set_ylabel('')
    ax.set_xlabel('Generation')
    setSize(ax, fontsize=fontsize - 2)
    if title is not None:
        ax.set_title('{}:{}'.format(pos[0], pos[1]), fontsize=fontsize)
    ax.xaxis.grid(True, linewidth=6)

def getNameColorMarker(df):
    f = lambda x: x.method.replace('HMM', r'$\mathcal{H}$').replace('MarkovChain', r'$\mathcal{M}')
    # + '$,\pi=$' + str(int(x.q * 100))
    # f = lambda x: x.method.replace('HMM', r'$\mathcal{H}$')
    cols = ['method']
    if 'q' in df.index.names:
        cols = ['q'] + cols
    names = df.unstack('S').reset_index()[cols].drop_duplicates()
    names['name'] = names.apply(f, axis=1)
    names = names.set_index(cols).sort_index(level='method')
    names['marker'] = getMarker(names.shape[0])
    names['color'] = getColorMap(names.shape[0])
    return names


def plotOnePower(df, info, axes, legendSubplot=-1, fontsize=7, markersize=5, ylabel='Hard', panel=list('ABC')):
    for j, (name, dff) in enumerate(df.groupby(level='coverage')):
        dff = dff.unstack('S')
        dff = dff.sortlevel(['method'], ascending=True)
        names = info.loc[dff.reset_index('coverage').index]
        dff.index = names.name

        dff.T.plot(ax=axes[j], legend=False, color=names.color.tolist(), style=names.marker.tolist(),
                   markersize=markersize)
        axes[j].axhline(y=5, color='k');
        setTicks(dff)
        if j == legendSubplot:
            handles, labels = axes[j].get_legend_handles_labels()
            axes[j].legend(handles[::-1], labels[::-1], loc='center left', fontsize=fontsize)
        if name == np.inf:
            name = r'$\infty$'
        else:
            name = '{:.0f}'.format(name)
        if ylabel == 'Hard': axes[j].set_title(r'$\lambda=$' + name, fontsize=fontsize)
        axes[j].set_xlabel(r'$s$')
        axes[j].set_ylabel(r'Power ({} Sweep)'.format(ylabel))
        setSize(axes[j], fontsize=fontsize)


def setSize(ax, fontsize=5):
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(fontsize)
    try:
        for item in ([ax.zaxis.label] + ax.get_zticklabels()):
            item.set_fontsize(fontsize)
    except:
        pass


def setLegendSize(ax, fontsize=5):
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, prop={'size': fontsize}, loc='best')

def setTicks(df):
    plt.xticks(df.columns.values);
    plt.xlim([0.018, 0.105]);
    plt.ylim([-2.0, 105]);
    plt.yticks(np.sort(np.append(np.arange(20, 101, 20), [5])))
    plt.xlabel('')


def savefig(name, dpi,path=utl.paperFiguresPath):
    # plt.tight_layout()
    plt.savefig(path+ name + '.pdf');
    plt.savefig(path+ name + '.tiff', dpi=dpi)


def plotQuantile(df, kde):
    from Utils import Util as utl
    quantiles = np.sort(np.append(np.linspace(0.0, 1, 1000)[:-1], np.linspace(0.999, 1, 10)))
    qq = pd.concat([utl.getQantilePvalues(df.COMALE, kde, quantiles=quantiles),
                    utl.getQantilePvalues(df.COMALENC, kde, quantiles=quantiles)], axis=1);
    qq.columns = ['data', 'null'];
    QQPval(qq, fname=utl.paperFiguresPath + 'qq.pdf')
