"""
Minimal plotting functions implemented through matplotlib.
"""

#############  PACKAGES  #############

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

#matplotlib.rcParams['pdf.fonttype'] = 42

from matplotlib.ticker import FixedLocator
from scipy.stats.kde   import gaussian_kde
from copy import deepcopy

############# PARAMETERS #############

#matplotlib.rcParams['pdf.fonttype'] = 42

blues = {
    3   :   ['#3182bd', '#9ecae1', '#deebf7'],
    4   :   ['#525252', '#969696', '#cccccc', '#f7f7f7'],
    5   :   ['#252525', '#636363', '#969696', '#cccccc', '#f7f7f7'],
    6   :   ['#252525', '#636363', '#969696', '#bdbdbd', '#d9d9d9', '#f7f7f7'],
    7   :   ['#252525', '#525252', '#737373', '#969696', '#bdbdbd', '#d9d9d9', '#f7f7f7'],
    8   :   ['#252525', '#525252', '#737373', '#969696', '#bdbdbd', '#d9d9d9', '#f0f0f0', '#ffffff']
    }
    
greys = {
    3   :   ['#636363', '#bdbdbd', '#f0f0f0'],
    4   :   ['#2171b5', '#6baed6', '#bdd7e7', '#eff3ff'],
    5   :   ['#08519c', '#3182bd', '#6baed6', '#bdd7e7', '#eff3ff'],
    6   :   ['#08519c', '#3182bd', '#6baed6', '#9ecae1', '#c6dbef', '#eff3ff'],
    7   :   ['#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#eff3ff'],
    8   :   ['#084594', '#2171b5', '#4292c6', '#6baed6', '#9ecae1', '#c6dbef', '#deebf7', '#f7fbff']
    }

defcolor     = '#252525'
defplotcolor = '#2171b5'
defpaircolor = ['#2171b5', '#9ecae1']

def cm2inch(x): return float(x)/2.54

goldr = (1.0 + np.sqrt(5)) / 2.0

nbins        = 50
aspect       = 1.0
plotwidth    = cm2inch(8.90)
offset       = 0.02
offsetpt     = 1.0
sublabelx    = -0.10
sublabely    = 0.95
ticklength   = 3
tickpad      = 3
axwidth      = 0.4
sizedot      = 8.
sizeline     = 0.8
textcolor    = defcolor

# paper
fontfamily   = 'Arial'
sizesublabel = 8+1
sizetext     = 6
sizelabel    = 6
sizetick     = 6
smallsizedot = 6.

## grant
#fontfamily   = 'Arial'
#sizesublabel = 12
#sizetext     = 8
#sizelabel    = 8
#sizetick     = 8
#smallsizedot = 6. * 2
#sizeline     = 1

##slides/poster
#fontfamily   = 'Avenir'
#sizesublabel = 18
#sizetext     = 18
#sizelabel    = 18
#sizetick     = 18
#smallsizedot = 6. * 4
#axwidth      = 1.5
#sizeline     = 3.0

def_labelprops = {
    'family' : fontfamily,
    'size'   : sizelabel,
    'color'  : textcolor
    }

def_sublabelprops = {
    'family' : fontfamily,
    'size'   : sizesublabel,
    'weight' : 'bold',
    'ha'     : 'center',
    'va'     : 'center',
    'color'  : 'k'
    }

def_ticklabelprops = {
    'family' : fontfamily,
    'size'   : sizetick,
    'color'  : textcolor
    }

def_axprops = {
    'linewidth' : axwidth,
    'linestyle' : '-',
    'color'     : textcolor
    }

def_tickprops = {
    'length'    : ticklength,
    'width'     : axwidth/2,
    'pad'       : tickpad,
    'axis'      : 'both',
    'direction' : 'out',
    'colors'    : textcolor,
    'bottom'    : True,
    'left'      : True,
    'top'       : False,
    'right'     : False
    }

def_minortickprops = {
    'length'    : ticklength-1.25,
    'width'     : axwidth/2,
    'axis'      : 'both',
    'direction' : 'out',
    'which'     : 'minor',
    'colors'    : textcolor
    }

def_legendprops = {
    'loc'           : 1,
    'frameon'       : False,
    'scatterpoints' : 1,
    'handletextpad' : -0.1,
    'prop'          : {'size' : sizelabel}
    }

def_scatterprops = {
    'lw'        : 0,
    's'         : sizedot,
    'marker'    : 'o'
    }

# Define a new function "callout" which circles points on a scatter plot and labels them
# Need properties for the arrows, the scatter points, and the text
def_callscatterprops = {
    'lw'     : sizeline,
    's'      : sizedot,
    'marker' : 'o'
    }

def_callarrowprops = {
    'lw'     : sizeline,
    's'      : sizedot,
    'marker' : 'o'
    }

def_errorprops = {
    'mew'        : 0,
    'markersize' : smallsizedot/2,
    'fmt'        : 'o',
    'elinewidth' : sizeline/2,
    'capthick'   : 0,
    'capsize'    : 0
    }

def_lineprops = {
    'lw'        : sizeline,
    'ls'        : '-',
    #'drawstyle' : 'steps'
    }

def_fillprops = {
    'lw'          : sizeline,
    'linestyle'   : '-',
    'interpolate' : True
    }

# NEW STYLE SOLID
def_histprops = {
    'histtype'  : 'bar',
    'lw'        : sizeline/2,
    'ls'        : 'solid',
    'edgecolor' : defcolor,
    'align'     : 'left'
    }

# OLD STYLE OUTLINE
#def_histprops = {
#    'histtype'  : 'step',#'bar',
#    'lw'        : sizeline,#/2,
#    'ls'        : 'solid',
#    #'edgecolor' : defcolor
#    }

def_barprops = {
    'lw'          : sizeline/2,
    'width'       : 1,
    'edgecolor'   : defcolor,
    'align'       : 'center', #other option: edge
    'orientation' : 'vertical'
    }

def_violinprops = {
    'vert'        : False,
    'widths'      : 0.5,    # Maximum height for violin, ~ half of vertical space
    'showmeans'   : True,
    'showextrema' : False,
    'showmedians' : False,
    'points'      : 100,    # number of points to evaluate for Gaussian KDE
    'bw_method'   : 'silverman'
    }

def_hexbinprops = {
    'gridsize' : (50,50),
    'mincnt'   : 1,
    'xscale'   : 'linear',
    'yscale'   : 'linear',
    'lw'       : axwidth/4.
    }

def_contourprops = {
    'linewidths' : sizeline/2,
    'linestyles' : 'solid'
    }

def_vplineprops = {
    'lw'    : axwidth,
    'color' : defcolor
    }

def_figprops = {
    'transparent' : True,
    'bbox_inches' : 'tight'
    }

def_tickprops_boxed = {
    'length'    : ticklength/1.75,
    'width'     : axwidth/1.75,
    'pad'       : tickpad,
    'axis'      : 'both',
    'direction' : 'in',
    'colors'    : textcolor,
    'bottom'    : True,
    'left'      : True,
    'top'       : True,
    'right'     : True
    }

def_minortickprops_boxed = {
    'length'    : (ticklength/1.75)-1.25,
    'width'     : axwidth/1.75,
    'axis'      : 'both',
    'direction' : 'in',
    'which'     : 'minor',
    'colors'    : textcolor
    }

def_axprops_ppt = {
    'linewidth' : 2 * axwidth,
    'linestyle' : '-',
    'color'     : textcolor
    }

def_tickprops_ppt = {
    'length'    : ticklength,
    'width'     : axwidth,
    'pad'       : tickpad,
    'axis'      : 'both',
    'direction' : 'out',
    'colors'    : textcolor,
    'bottom'    : True,
    'left'      : True,
    'top'       : False,
    'right'     : False
    }

def_minortickprops_ppt = {
    'length'    : ticklength-1.25,
    'width'     : axwidth,
    'axis'      : 'both',
    'direction' : 'out',
    'which'     : 'minor',
    'colors'    : textcolor
    }

singlevartypes = ['hist', 'kde', 'violin']


############# PLOTTING  FUNCTIONS #############

def plot(**pdata):
    """ Generic plotting routine. Sets basic options then passes parameters to detailed plotting functions. """
    
    # Sanity checks
    
    ndata = 0
    
    assert 'ax' not in pdata or 'save' not in pdata, 'Attempted to save the figure, but the axis has been supplied manually!'
    
    assert ('x' in pdata or 'y' in pdata), 'No data passed to plot!'
    
    if 'x' in pdata and 'y' in pdata:
        ndata = len(np.array(pdata['x']).shape)
        assert ndata==len(np.array(pdata['y']).shape), 'x (dimension %d) and y (dimension %d) have mismatched shapes!' % (ndata, len(np.array(pdata['y']).shape))
        if ndata==1 and not ('type' in pdata and pdata['type']=='circos') and len(np.array(pdata['x'][0]))==1:
            pdata['x'] = [pdata['x']]
            pdata['y'] = [pdata['y']]
    
    elif 'x' in pdata and 'y' not in pdata:
        ndata = len(np.array(pdata['x']).shape)
        assert pdata['type'] in singlevartypes, 'Only one set of values passed, but plot type (%s) requires two sets of values to plot!' % pdata['type']
        if ndata==1 and len(np.array(pdata['x'][0]))==1: pdata['x'] = [pdata['x']]
    
    elif 'y' in pdata and 'x' not in pdata:
        ndata = len(np.array(pdata['y']).shape)
        assert pdata['type'] in singlevartypes, 'Only one set of values passed, but plot type (%s) requires two sets of values to plot!' % pdata['type']
        pdata['x'] = [k for k in pdata['y']]
        pdata['y'] = []
        if ndata==1 and len(np.array(pdata['x'][0]))==1: pdata['x'] = [pdata['x']]

    if 'colors' not in pdata: assert len(pdata['x'])<9, 'No colors provided, but the number of elements to be plotted (%d) is >8, maximum number of default colors!' % len(pdata['x'])

    # If no axis is passed, create the axis
    
    fig = 0

    if 'ax' not in pdata:
        if 'aspect' not in pdata and 'dims' not in pdata:
            fig         = plt.figure(figsize=(aspect * plotwidth, plotwidth))
            pdata['ax'] = plt.subplot(111)
        elif 'aspect' in pdata:
            fig         = plt.figure(figsize=(pdata['aspect'] * plotwidth, plotwidth))
            pdata['ax'] = plt.subplot(111)
        elif 'dims' in pdata:
            fig         = plt.figure(figsize=(pdata['dims'][0], pdata['dims'][1]))
            pdata['ax'] = plt.subplot(111)

    # Fill style parameters if not passed in pdata

    if 'theme' not in pdata: pdata['theme'] = 'open'

    if 'type'           not in pdata: pdata['type']           = 'scatter'
    if 'ticklabelprops' not in pdata: pdata['ticklabelprops'] = def_ticklabelprops

    if 'axprops' not in pdata:
        if 'ppt' in pdata and pdata['ppt']: pdata['axprops'] = def_axprops_ppt
        else:                               pdata['axprops'] = def_axprops
    if 'tickprops' not in pdata:
        if 'ppt' in pdata and pdata['ppt']: pdata['tickprops'] = def_tickprops_ppt
        elif pdata['theme']=='boxed':       pdata['tickprops'] = def_tickprops_boxed
        else:                               pdata['tickprops'] = def_tickprops
    if 'minortickprops' not in pdata:
        if 'ppt' in pdata and pdata['ppt']: pdata['minortickprops'] = def_minortickprops_ppt
        elif pdata['theme']=='boxed':       pdata['minortickprops'] = def_minortickprops_boxed
        else:                               pdata['minortickprops'] = def_minortickprops

    if 'labelprops'     not in pdata:                         pdata['labelprops']     = def_labelprops
    if 'sublabelprops'  not in pdata and 'sublabel' in pdata: pdata['sublabelprops']  = def_sublabelprops
    if 'sublabelcoords' not in pdata and 'sublabel' in pdata: pdata['sublabelcoords'] = [sublabelx, sublabely]
    if 'figprops'       not in pdata and 'save' in pdata:     pdata['figprops']       = def_figprops
    if 'legendprops'    not in pdata and 'legend' in pdata:   pdata['legendprops']    = def_legendprops
    
    if 'colors' not in pdata:
        if   len(pdata['x'])==1: pdata['colors'] = [defplotcolor]
        elif len(pdata['x'])==2: pdata['colors'] = [k for k in defpaircolor]
        else:                    pdata['colors'] = [k for k in blues[len(pdata['x'])]]

    if 'plotprops' not in pdata:
        if pdata['type']=='scatter': pdata['plotprops'] = def_scatterprops
        if pdata['type']=='error':   pdata['plotprops'] = def_errorprops
        if pdata['type']=='line':    pdata['plotprops'] = def_lineprops
        if pdata['type']=='fill':    pdata['plotprops'] = def_fillprops
        if pdata['type']=='hist':    pdata['plotprops'] = def_histprops
        if pdata['type']=='bar':     pdata['plotprops'] = def_barprops
        if pdata['type']=='kde':     pdata['plotprops'] = def_lineprops
        if pdata['type']=='circos':  pdata['plotprops'] = def_lineprops
        if pdata['type']=='violin':  pdata['plotprops'] = def_violinprops
        if pdata['type']=='hexbin':  pdata['plotprops'] = def_hexbinprops
        if pdata['type']=='contour': pdata['plotprops'] = def_contourprops

    # Fill in x axis limits and ticks if not passed in pdata
    
    if pdata['type']=='circos' and 'xticks' not in pdata: pdata['xticks'] = []
    if pdata['type']=='circos' and 'yticks' not in pdata: pdata['yticks'] = []

    if 'xlim' not in pdata: pdata['xlim'] = [np.min([np.min(x) for x in pdata['x']]), np.max([np.max(x) for x in pdata['x']])]

    if 'xticks' not in pdata:
        xtick, xminortick = [], []
        if 'logx' in pdata and pdata['logx']: xtick, xminortick = getlogticks(pdata['xlim'])
        else:                                 xtick, xminortick = getticks(pdata['xlim'])
        pdata['xticks'] = xtick
        if 'xminorticks' not in pdata: pdata['xminorticks'] = xminortick
    if 'xminorticks' not in pdata: pdata['xminorticks'] = []

    if 'bins' not in pdata and pdata['type']=='hist':
        width         = (pdata['xlim'][1]-pdata['xlim'][0])/nbins
        bins          = np.arange(pdata['xlim'][0], pdata['xlim'][1]+width, width)
        pdata['bins'] = bins

    if 'combine' not in pdata and pdata['type']=='hist': pdata['combine'] = False

    # Fill in y axis limits and ticks

    if 'ylim' not in pdata and pdata['type'] not in singlevartypes:
        pdata['ylim'] = [np.min([np.min(y) for y in pdata['y']]), np.max([np.max(y) for y in pdata['y']])]
    elif 'ylim' not in pdata:
        ymin = 1.0
        ymax = 0.0
        if pdata['type']=='hist':
            for i in range(len(pdata['x'])):
                w = np.ones_like(pdata['x'][i])/float(len(pdata['x'][i]))
                h, bin_edges = np.histogram(pdata['x'][i], weights=w, range=pdata['xlim'], bins=pdata['bins'])
                if np.min(h)<ymin: ymin = np.min(h)
                if np.max(h)>ymax: ymax = np.max(h)
        elif pdata['type']=='kde':
           for i in range(len(pdata['x'])):
                ypdf = gaussian_kde(pdata['x'][i], bw_method='silverman')
                x    = np.linspace(np.min(pdata['xticks']), np.max(pdata['xticks']), 300)
                if np.min(ypdf(x))<ymin: ymin = np.min(ypdf(x))
                if np.max(ypdf(x))>ymax: ymax = np.max(ypdf(x))
        elif pdata['type']=='violin':
            ymax = len(pdata['x'])
            if 'positions' in pdata:
                ymin = np.min(pdata['positions'])
                ymax = np.max(pdata['positions'])
        pdata['ylim'] = [ymin, ymax]

    if 'yticks' not in pdata:
        ytick, yminortick = [], []
        if 'logy' in pdata and pdata['logy']: ytick, yminortick = getlogticks(pdata['ylim'])
        else:                                 ytick, yminortick = getticks(pdata['ylim'])
        pdata['yticks'] = ytick
        if 'yminorticks' not in pdata: pdata['yminorticks'] = yminortick
    if 'yminorticks' not in pdata: pdata['yminorticks'] = []

    if pdata['type']=='fill' and 'y2' not in pdata: pdata['y2']=0

    # Sanity check for ticks and log scale

    if 'logx' in pdata and pdata['logx']:
        assert pdata['xlim'][0]>0, 'Attempting to plot x on a log scale, but the lower limit (%lf) is <= 0!' % pdata['xlim'][0]
        if 'xticks' not in pdata:
            xtick, xminortick = getlogticks(pdata['xlim'])
            pdata['xticks'] = xtick
            if 'xminorticks' not in pdata: pdata['xminorticks'] = xminortick
    else: pdata['logx'] = False

    if 'logy' in pdata and pdata['logy'] and pdata['type'] not in singlevartypes:
        assert pdata['ylim'][0]>0, 'Attempting to plot y on a log scale, but the lower limit (%lf) is <= 0!' % pdata['ylim'][0]
        if 'yticks' not in pdata:
            ytick, yminortick = getlogticks(pdata['ylim'])
            pdata['yticks'] = ytick
            if 'yminorticks' not in pdata: pdata['yminorticks'] = yminortick
    else: pdata['logy'] = False

    # Make plot

    if   pdata['type']=='scatter': scatter(**pdata)
    elif pdata['type']=='error':   error(**pdata)
    elif pdata['type']=='line':    line(**pdata)
    elif pdata['type']=='fill':    fill(**pdata)
    elif pdata['type']=='hist':    hist(**pdata)
    elif pdata['type']=='bar':     bar(**pdata)
    elif pdata['type']=='kde':     kde(**pdata)
    elif pdata['type']=='circos':  circos(**pdata)
    elif pdata['type']=='violin':  violin(**pdata)
    elif pdata['type']=='hexbin':  hexbin(**pdata)
    elif pdata['type']=='contour': contour(**pdata)

    # Set plot appearance and plot axes
    
    if ('orientation' in pdata['plotprops']) and (pdata['plotprops']['orientation']=='horizontal'):
        temp                 = [k for k in pdata['xlim']]
        pdata['xlim']        = [k for k in pdata['ylim']]
        pdata['ylim']        = [k for k in temp]
        temp                 = [k for k in pdata['xticks']]
        pdata['xticks']      = [k for k in pdata['yticks']]
        pdata['yticks']      = [k for k in temp]
        temp                 = [k for k in pdata['xminorticks']]
        pdata['xminorticks'] = [k for k in pdata['yminorticks']]
        pdata['yminorticks'] = [k for k in temp]

    if pdata['theme']=='open':
        adjustlim(pdata['xlim'], logscale=pdata['logx'])
        adjustlim(pdata['ylim'], logscale=pdata['logy'])
    if 'nudgex' in pdata: pdata['xlim'][1] *= pdata['nudgex']
    if 'nudgey' in pdata: pdata['ylim'][1] *= pdata['nudgey']
    setappearance(**pdata)

    # Save

    if 'save' in pdata:
        plt.savefig(pdata['save']+'.pdf', **pdata['figprops'])
        plt.close(fig)


def scatter(**pdata):
    """ Generic scatter plot. """

    # Plot data
    
    hollowprops = pdata['plotprops'].copy()
    hollowprops['lw']        = 1
    hollowprops['facecolor'] = 'none'
    
    for i in range(len(pdata['x'])):
        x = pdata['x'][i]
        y = pdata['y'][i]
        
        if 'facecolor' in pdata and pdata['facecolor'][i] and 'edgecolor' in pdata and pdata['edgecolor'][i]:
            cf = pdata['facecolor'][i]
            ce = pdata['edgecolor'][i]
            if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].scatter(x, y,               edgecolor=ce, **hollowprops)
            else:                                        pdata['ax'].scatter(x, y, facecolor=cf, edgecolor=ce, **pdata['plotprops'])
        
        elif 'colors' in pdata and pdata['colors'][i]:
            c = pdata['colors'][i]
            if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].scatter(x, y,              edgecolor=c, **hollowprops)
            else:                                        pdata['ax'].scatter(x, y, facecolor=c, edgecolor=c, **pdata['plotprops'])

        else:
            c = defcolor
            if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].scatter(x, y,              edgecolor=c, **hollowprops)
            else:                                        pdata['ax'].scatter(x, y, facecolor=c, edgecolor=c, **pdata['plotprops'])

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
        
            if 'facecolor' in pdata and pdata['facecolor'][i] and 'edgecolor' in pdata and pdata['edgecolor'][i]:
                cf = pdata['facecolor'][i]
                ce = pdata['edgecolor'][i]
                l  = pdata['legend'][i]
                if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].scatter(x, y,               edgecolor=ce, label=l, zorder=len(pdata['x'])-i+10, **hollowprops)
                else:                                        pdata['ax'].scatter(x, y, facecolor=cf, edgecolor=ce, label=l, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])
            
            else:
                c = defcolor
                if 'colors' in pdata and pdata['colors'][i]: c = pdata['colors'][i]
                l = pdata['legend'][i]
                if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].scatter(x, y,              edgecolor=c, label=l, zorder=len(pdata['x'])-i+10, **hollowprops)
                else:                                        pdata['ax'].scatter(x, y, facecolor=c, edgecolor=c, label=l, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])
                
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)

        #for text in l.get_texts(): text.set_size(sizelabel)


def error(**pdata):
    """ Generic errorbar plot. """

    # Plot data
    
    hollowprops = pdata['plotprops'].copy()
    hollowprops['mew'] = 1
    hollowprops['mfc'] = 'none'

    for i in range(len(pdata['x'])):
        x    = pdata['x'][i]
        y    = pdata['y'][i]
        xerr = None
        yerr = None
        
        if 'xerr' in pdata and np.any(pdata['xerr']) and np.any(pdata['xerr'][i]): xerr = pdata['xerr'][i]
        if 'yerr' in pdata and np.any(pdata['yerr']) and np.any(pdata['yerr'][i]): yerr = pdata['yerr'][i]
        
        if 'facecolor' in pdata and pdata['facecolor'][i] and 'edgecolor' in pdata and pdata['edgecolor'][i]:
            cf = pdata['facecolor'][i]
            ce = pdata['edgecolor'][i]
            if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].errorbar(x, y, xerr=xerr, yerr=yerr,         mec=ce, color=ce, zorder=len(pdata['x'])-i+10, **hollowprops)
            else:                                        pdata['ax'].errorbar(x, y, xerr=xerr, yerr=yerr, mfc=cf, mec=ce, color=ce, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])
        
        elif 'colors' in pdata and pdata['colors'][i]:
            c = pdata['colors'][i]
            if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].errorbar(x, y, xerr=xerr, yerr=yerr,        mec=c, color=c, zorder=len(pdata['x'])-i+10, **hollowprops)
            else:                                        pdata['ax'].errorbar(x, y, xerr=xerr, yerr=yerr, mfc=c, mec=c, color=c, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])

        else:
            c = defcolor
            if 'hollow' in pdata and pdata['hollow'][i]: pdata['ax'].errorbar(x, y, xerr=xerr, yerr=yerr,        mec=c, color=c, zorder=len(pdata['x'])-i+10, **hollowprops)
            else:                                        pdata['ax'].errorbar(x, y, xerr=xerr, yerr=yerr, mfc=c, mec=c, color=c, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])

#    for i in range(len(pdata['x'])):
#        c = pdata['colors'][i]
#        x = pdata['x'][i]
#        y = pdata['y'][i]
#        
#        if 'xerr' in pdata:
#            if 'yerr' in pdata: pdata['ax'].errorbar(x, y, xerr=pdata['xerr'][i], yerr=pdata['yerr'][i], c=c, mfc=c, mec=c, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])
#            else:               pdata['ax'].errorbar(x, y, xerr=pdata['xerr'][i],                        c=c, mfc=c, mec=c, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])
#        elif 'yerr' in pdata:   pdata['ax'].errorbar(x, y,                        yerr=pdata['yerr'][i], c=c, mfc=c, mec=c, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])
#        else:                   pdata['ax'].errorbar(x, y,                                               c=c, mfc=c, mec=c, zorder=len(pdata['x'])-i+10, **pdata['plotprops'])

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].scatter(x, y, facecolor=c, edgecolor=c, label=l, **def_scatterprops)
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)


def line(**pdata):
    """ Generic line plot. """

    # Plot data

    for i in range(len(pdata['x'])):
        c = pdata['colors'][i]
        x = pdata['x'][i]
        y = pdata['y'][i]
        
        if 'zorder' not in pdata['plotprops']: pdata['ax'].plot(x, y, color=c, zorder=len(pdata['x'])-i, **pdata['plotprops'])
        else:                                  pdata['ax'].plot(x, y, color=c,                           **pdata['plotprops'])
        
#        if 'legend' in pdata:
#            pdata['ax'].plot(x, y, color=c, zorder=len(pdata['x'])-i, label=pdata['legend'][i], **pdata['plotprops'])
#        
#        else:
#            pdata['ax'].plot(x, y, color=c, zorder=len(pdata['x'])-i, **pdata['plotprops'])

    # Make legend (optional)

    if 'legend' in pdata:# and 'plotlegend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].plot(x, y, color=c, label=l, **pdata['plotprops'])
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)


def fill(**pdata):
    """ Generic filled line plot. """

    # Plot data

    for i in range(len(pdata['x'])):
        c = pdata['colors'][i]
        x = pdata['x'][i]
        y = pdata['y'][i]
        
        if 'zorder' not in pdata['plotprops']: pdata['ax'].fill_between(x, y, color=c, zorder=len(pdata['x'])-i, **pdata['plotprops'])
        else:                                  pdata['ax'].fill_between(x, y, color=c                            **pdata['plotprops'])
    
    # Make legend (optional)

    if 'legend' in pdata:# and 'plotlegend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].plot(x, y, color=c, label=l, **pdata['plotprops'])
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)


def hist(**pdata):
    """ Generic histogram. """

    # Plot data

    if pdata['combine']: n, bins, patches = pdata['ax'].hist(pdata['x'], weights=[np.ones_like(x)/float(len(x)) for x in pdata['x']], color=pdata['colors'], bins=pdata['bins'], **pdata['plotprops'])

    else:
        for i in range(len(pdata['x'])):
            c = pdata['colors'][i]
            x = pdata['x'][i]
            w = np.ones_like(x)/float(len(x))
            if 'weights' in pdata: w = pdata['weights'][i]
            
            if 'percent' in pdata and pdata['percent']: w *= 100.0
        
            n, bins, patches = pdata['ax'].hist(x, weights=w, facecolor=c, zorder=len(pdata['x'])-i, bins=pdata['bins'], range=(pdata['xlim'][0], pdata['xlim'][1]), **pdata['plotprops'])

            for patch in patches: patch.set_facecolor(c)

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].scatter(x, y, facecolor=c, edgecolor=c, label=l, **def_scatterprops)
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)

    # Annotate (optional)

    if 'annotate' in pdata:
        for i in range(len(pdata['annotate'])):
            pdata['ax'].text(pdata['annotate'][i][0], pdata['annotate'][i][1], pdata['annotate'][i][2], **pdata['ticklabelprops'])


def bar(**pdata):
    """ Generic bar graph. """
    
    # Check alignment

    isHorizontal = False
    if pdata['plotprops']['orientation'] and pdata['plotprops']['orientation']=='horizontal':
        isHorizontal = True
        if 'width' in pdata['plotprops']:
            pdata['plotprops']['height'] = pdata['plotprops']['width']
            del pdata['plotprops']['width']
        del pdata['plotprops']['orientation']
    else: pdata['plotprops']['orientation'] = 'vertical'
    
    # Plot data

    for i in range(len(pdata['x'])):
        c = pdata['colors'][i]
        x = pdata['x'][i]
        y = pdata['y'][i]
        
        xerr   = None
        yerr   = None
        err_kw = {}
        if 'xerr' in pdata and np.any(pdata['xerr']) and np.any(pdata['xerr'][i]): xerr = pdata['xerr'][i]
        if 'yerr' in pdata and np.any(pdata['yerr']) and np.any(pdata['yerr'][i]): yerr = pdata['yerr'][i]
        if xerr or yerr: err_kw = { 'ecolor' : defcolor, 'capsize' : 0 }
        
        if isHorizontal: pdata['ax'].barh(x, y, xerr=xerr, yerr=yerr, color=c, zorder=len(pdata['x'])-1-i, error_kw=err_kw, **pdata['plotprops'])
        else:            pdata['ax'].bar( x, y, xerr=xerr, yerr=yerr, color=c, zorder=len(pdata['x'])-1-i, error_kw=err_kw, **pdata['plotprops'])

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].scatter(x, y, facecolor=c, edgecolor=c, label=l, **def_scatterprops)
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)

    # Annotate (optional)

    if 'annotate' in pdata:
        for i in range(len(pdata['annotate'])):
            pdata['ax'].text(pdata['annotate'][i][0], pdata['annotate'][i][1], pdata['annotate'][i][2], **pdata['ticklabelprops'])

    if isHorizontal: pdata['plotprops']['orientation'] = 'horizontal'


def kde(**pdata):
    """ Generic kde-smoothed plot. """

    # Plot data

    for i in range(len(pdata['x'])):
        c    = pdata['colors'][i]
        y    = pdata['x'][i]
        ypdf = gaussian_kde(y, bw_method='silverman')
        x    = np.linspace(np.min(pdata['xticks']), np.max(pdata['xticks']), 300)
        
        pdata['ax'].plot(x, ypdf(x), color=c, zorder=len(pdata['x'])-i, **pdata['plotprops'])

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].plot(x, y, color=c, label=l, **pdata['plotprops'])
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)


def circos(**pdata):
    """ Generic circos-style plot. """
    
    # Generate Bezier curves using index to angle map
    
    if 't' not in pdata:
        pdata['t'] = np.arange(0., 1.01, 0.01)
    if 'rad' not in pdata:
        pdata['rad'] = [[0.8, 0.8] for i in range(len(pdata['x']))]
    if 'size' not in pdata:
        pdata['size'] = float(np.max([np.max(pdata['x']), np.max(pdata['y'])]))
    if 'bezrad' not in pdata:
        pdata['bezrad'] = 0.
    if not hasattr(pdata['bezrad'], '__len__'):
        temp = pdata['bezrad']
        pdata['bezrad'] = [temp for i in range(len(pdata['rad']))]
    if 'angle' not in pdata:
        pdata['angle'] = [[-2 * np.pi * pdata['x'][i] / pdata['size'], -2 * np.pi * pdata['y'][i] / pdata['size']] for i in range(len(pdata['x']))]
    
    cpolar = bezier(pdata['rad'], pdata['angle'], pdata['t'], pdata['bezrad'])

    # Plot data

    for i in range(len(pdata['x'])):
        c  = pdata['colors'][i]
        x  = cpolar[i][0]
        y  = cpolar[i][1]
        
        if 'arcprops' in pdata: pdata['ax'].plot(x, y, color=c, zorder=len(pdata['x'])-i, **pdata['arcprops'][i])
        else:                   pdata['ax'].plot(x, y, color=c, zorder=len(pdata['x'])-i, **pdata['plotprops'])

    # Add labels (optional)

    if 'ticklabels' in pdata:
        for tl in pdata['ticklabels']:
            t = -2 * np.pi * tl[0] / pdata['size']
            r = 0.85
            l = tl[0]
            if len(tl)>1: r = tl[1]
            if len(tl)>2: l = tl[2]
            if t>-np.pi/2. or t<-3.*np.pi/2.: pdata['ax'].text(r * np.cos(t), r * np.sin(t), l, rotation_mode='anchor', rotation=rad2deg(t),      ha='left',  va='center', **pdata['ticklabelprops'])
            else:                             pdata['ax'].text(r * np.cos(t), r * np.sin(t), l, rotation_mode='anchor', rotation=rad2deg(t)+180., ha='right', va='center', **pdata['ticklabelprops'])

    if 'textlabels' in pdata:
        for tl in pdata['textlabels']:
            t = -2 * np.pi * tl[0] / pdata['size']
            r = tl[1]
            l = tl[2]
            if t>-np.pi/2. or t<-3.*np.pi/2.: pdata['ax'].text(r * np.cos(t), r * np.sin(t), l, rotation_mode='anchor', rotation=rad2deg(t),      ha='left',  va='center', **pdata['labelprops'])
            else:                             pdata['ax'].text(r * np.cos(t), r * np.sin(t), l, rotation_mode='anchor', rotation=rad2deg(t)+180., ha='right', va='center', **pdata['labelprops'])

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].plot(x, y, color=c, label=l, **pdata['plotprops'])
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)

        #for text in l.get_texts(): text.set_size(sizelabel+1)


def violin(**pdata):
    """ Generic violin plot. """

    # Plot data

    for i in range(len(pdata['x'])):
        x   = pdata['x'][i]
        pos = [i+1]
        if 'positions' in pdata: pos = pdata['positions'][i]
        
        #print(pos)
        #print(len(x))
        
        #vp = pdata['ax'].violinplot(x, **pdata['plotprops'])
        vp = pdata['ax'].violinplot(x, positions=pos, **pdata['plotprops'])
    
        # Adjust appearance
    
        c  = defcolor
        fc = defcolor
        ec = defcolor
        
        if 'colors' in pdata:
            c  = pdata['colors'][i]
            fc = pdata['colors'][i]
            ec = pdata['colors'][i]
    
        if 'facecolor' in pdata:
            fc = pdata['facecolor'][i]
        
        if 'edgecolor' in pdata:
            fc = pdata['edgecolor'][i]
    
        for body in vp['bodies']:
            body.set_facecolor(fc)
            body.set_edgecolor(ec)
            body.set_linewidth(axwidth)

        showLines = ['cmins', 'cmaxes', 'cmeans', 'cmedians']
        hideLines = ['cbars']

        for line in [vp[l] for l in showLines if l in vp]:
            line.set_color(c)
            line.set_linewidth(axwidth)
        for line in [vp[l] for l in hideLines if l in vp]:
            line.set_linewidth(0)

    # Make legend (optional)

    if 'legend' in pdata:
        x = [10 * pdata['xlim'][1], 10 * pdata['xlim'][1]]
        y = [10 * pdata['ylim'][1], 10 * pdata['ylim'][1]]
        for i in range(len(pdata['legend'])):
            c = pdata['colors'][i]
            l = pdata['legend'][i]
            pdata['ax'].plot(x, y, color=c, label=l, **pdata['plotprops'])
        l = pdata['ax'].legend(**pdata['legendprops'])
        for text in l.get_texts(): text.set_color(textcolor)


def hexbin(**pdata):
    """ 2D histogram with hexagonal bins. """

    # Plot data

    for i in range(len(pdata['x'])):
        x = pdata['x'][i]
        y = pdata['y'][i]
        
        hb = pdata['ax'].hexbin(x, y, **pdata['plotprops'])

    # Annotate (optional)

    if 'annotate' in pdata:
        for i in range(len(pdata['annotate'])):
            pdata['ax'].text(pdata['annotate'][i][0], pdata['annotate'][i][1], pdata['annotate'][i][2], **pdata['ticklabelprops'])


def contour(**pdata):
    """ Generic contour plot. """

    # Sort out contour shading plot props from line props
    
    cShade = ['cmap', 'hatches']
    cLine  = ['linewidths', 'linestyles']

    # Plot data

    for i in range(len(pdata['x'])):
        x = pdata['x'][i]
        y = pdata['y'][i]
        z = pdata['z'][i]
        
        tempprops = deepcopy(pdata['plotprops'])
        for prop in cLine:
            if prop in tempprops: del tempprops[prop]
        if 'colors' in tempprops and 'cmap' in tempprops: del tempprops['colors']
        
        pdata['ax'].contourf(x, y, z, **tempprops)
        
        tempprops = deepcopy(pdata['plotprops'])
        for prop in cShade:
            if prop in tempprops: del tempprops[prop]
        
        pdata['ax'].contour(x, y, z, **tempprops)

    # Annotate (optional)

    if 'annotate' in pdata:
        for i in range(len(pdata['annotate'])):
            pdata['ax'].text(pdata['annotate'][i][0], pdata['annotate'][i][1], pdata['annotate'][i][2], **pdata['ticklabelprops'])


def setappearance(**pdata):
    """ Set out general plot appearance (axis labels, tick parameters, etc). """

    # Make axis invisible

    if 'noaxes' in pdata and pdata['noaxes']:
        for axis in ['left', 'right', 'bottom', 'top']: pdata['ax'].spines[axis].set_visible(False)
    
    elif pdata['theme']=='boxed':
        for axis in ['left', 'right', 'bottom', 'top']: pdata['ax'].spines[axis].set_visible(True)

    elif pdata['theme']=='open':
        for axis in ['top', 'right']:   pdata['ax'].spines[axis].set_visible(False)
        for axis in ['bottom', 'left']: pdata['ax'].spines[axis].set_visible(True)

    if 'hide' in pdata:
        for axis in pdata['hide']: pdata['ax'].spines[axis].set_visible(False)

    if 'show' in pdata:
        for axis in pdata['show']: pdata['ax'].spines[axis].set_visible(True)
    
    # Add labels
    
    if 'xlabel' in pdata: pdata['ax'].set_xlabel(pdata['xlabel'], **pdata['labelprops'])
    if 'ylabel' in pdata: pdata['ax'].set_ylabel(pdata['ylabel'], **pdata['labelprops'])
    
    if 'sublabel' in pdata:
        pdata['sublabel'] = pdata['sublabel']#.upper()
        pdata['ax'].text(pdata['sublabelcoords'][0], pdata['sublabelcoords'][1], pdata['sublabel'], transform=pdata['ax'].transAxes, **pdata['sublabelprops'])

    # Set axis limits automatically, or from passed parameters
    
    #if 'logx' in pdata and pdata['logx']: pdata['ax'].set_xscale('log')
    #if 'logy' in pdata and pdata['logy']: pdata['ax'].set_yscale('log')
    
    if 'logx' in pdata and pdata['logx']: pdata['ax'].set_xscale('log', nonpositive='clip')  # GitHub
    if 'logy' in pdata and pdata['logy']: pdata['ax'].set_yscale('log', nonpositive='clip')
    
#    if 'logx' in pdata and pdata['logx']: pdata['ax'].set_xscale('log', nonposx='clip')  # Code Ocean
#    if 'logy' in pdata and pdata['logy']: pdata['ax'].set_yscale('log', nonposy='clip')
    
    if 'xlim' in pdata: pdata['ax'].set_xlim(pdata['xlim'][0], pdata['xlim'][1])
    if 'ylim' in pdata: pdata['ax'].set_ylim(pdata['ylim'][0], pdata['ylim'][1])
    
    # Set tick positions
    
    if 'xticks' in pdata: pdata['ax'].set_xticks(pdata['xticks'])
    if 'yticks' in pdata: pdata['ax'].set_yticks(pdata['yticks'])

    if 'xminorticks' in pdata: pdata['ax'].get_xaxis().set_minor_locator(FixedLocator(pdata['xminorticks']))
    if 'yminorticks' in pdata: pdata['ax'].get_yaxis().set_minor_locator(FixedLocator(pdata['yminorticks']))

    # Plot axes

    plotaxes(**pdata)

    # Set tick labels and properties
    
    if 'xticklabels' in pdata: pdata['ax'].set_xticklabels(pdata['xticklabels'])
    if 'yticklabels' in pdata: pdata['ax'].set_yticklabels(pdata['yticklabels'])
    
    if 'xminorticklabels' in pdata: pdata['ax'].set_xticklabels(pdata['xminorticklabels'], minor=True)
    if 'yminorticklabels' in pdata: pdata['ax'].set_yticklabels(pdata['yminorticklabels'], minor=True)

    if 'bottom' in pdata['tickprops'] and pdata['tickprops']['bottom']:
        if 'top' in pdata['tickprops'] and pdata['tickprops']['top']: pdata['ax'].xaxis.set_ticks_position('both')
        else:                                                         pdata['ax'].xaxis.tick_bottom()
    elif 'top' in pdata['tickprops'] and pdata['tickprops']['top']:   pdata['ax'].xaxis.tick_top()

    if 'left' in pdata['tickprops'] and pdata['tickprops']['left']:
        if 'right' in pdata['tickprops'] and pdata['tickprops']['right']: pdata['ax'].yaxis.set_ticks_position('both')
        else:                                                             pdata['ax'].yaxis.tick_left()
    elif 'right' in pdata['tickprops'] and pdata['tickprops']['right']:   pdata['ax'].yaxis.tick_right()

    for label in pdata['ax'].xaxis.get_majorticklabels(): label.set(**pdata['ticklabelprops'])
    for label in pdata['ax'].yaxis.get_majorticklabels(): label.set(**pdata['ticklabelprops'])

    for label in pdata['ax'].xaxis.get_minorticklabels(): label.set(**pdata['ticklabelprops'])
    for label in pdata['ax'].yaxis.get_minorticklabels(): label.set(**pdata['ticklabelprops'])

    # Set tick parameters
    
    pdata['ax'].tick_params(**pdata['tickprops'])
    pdata['ax'].tick_params(**pdata['minortickprops'])

    # Set axis properties

    for axis in ['top','bottom','left','right']:
        if 'color'     in pdata['axprops']: pdata['ax'].spines[axis].set_color(pdata['axprops']['color'])
        if 'linewidth' in pdata['axprops']: pdata['ax'].spines[axis].set_linewidth(pdata['axprops']['linewidth'])


def plotaxes(**pdata):
    """ Add axes to a plot. """
    
    xxmin = 0
    xxmax = 0
    yymin = 0
    yymax = 0
    
    # x axis start location
    
    if 'xaxstart' in pdata: xxmin = pdata['xaxstart']
    else:
        if len(pdata['ax'].get_xticks())>0:
            if len(pdata['ax'].get_xaxis().get_minorticklocs())>0: xxmin = np.min([np.min(pdata['ax'].get_xticks()), np.min(pdata['ax'].get_xaxis().get_minorticklocs())])
            else:                                                  xxmin = np.min(pdata['ax'].get_xticks())
        else:                                                      xxmin = np.min(pdata['ax'].get_xlim())

    # x axis end location
    
    if 'xaxend' in pdata: xxmax = pdata['xaxend']
    else:
        if len(pdata['ax'].get_xticks())>0:
            if len(pdata['ax'].get_xaxis().get_minorticklocs())>0: xxmax = np.max([np.max(pdata['ax'].get_xticks()), np.max(pdata['ax'].get_xaxis().get_minorticklocs())])
            else:                                                  xxmax = np.max(pdata['ax'].get_xticks())
        else:                                                      xxmax = np.max(pdata['ax'].get_xlim())

    # x axis height
    
    if 'xaxy' in pdata: xaxy = pdata['xaxy']
    else:               xaxy = np.min(pdata['ax'].get_ylim())

    # y axis start location

    if 'yaxstart' in pdata: yymin = pdata['yaxstart']
    else:
        if len(pdata['ax'].get_yticks())>0:
            if len(pdata['ax'].get_yaxis().get_minorticklocs())>0: yymin = np.min([np.min(pdata['ax'].get_yticks()), np.min(pdata['ax'].get_yaxis().get_minorticklocs())])
            else:                                                  yymin = np.min(pdata['ax'].get_yticks())
        else:                                                      yymin = np.min(pdata['ax'].get_ylim())

    # y axis end location

    if 'yaxend' in pdata: yymax = pdata['yaxend']
    else:
        if len(pdata['ax'].get_yticks())>0:
            if len(pdata['ax'].get_yaxis().get_minorticklocs())>0: yymax = np.max([np.max(pdata['ax'].get_yticks()), np.max(pdata['ax'].get_yaxis().get_minorticklocs())])
            else:                                                  yymax = np.max(pdata['ax'].get_yticks())
        else:                                                      yymax = np.max(pdata['ax'].get_ylim())

    # y axis position
    
    if 'yaxx' in pdata: yaxx = pdata['yaxx']
    else:               yaxx = np.min(pdata['ax'].get_xlim())

    # Get axes limits

    ax_xx = [xxmin, xxmax]
    ax_xy = [ xaxy,  xaxy]
    ax_yx = [ yaxx,  yaxx]
    ax_yy = [yymin, yymax]
    
    if 'mirrory' in pdata and pdata['mirrory']: ax_yx = [np.max(pdata['ax'].get_xlim()), np.max(pdata['ax'].get_xlim())]

    aspect = pdata['ax'].get_aspect()
    if aspect=='auto':
        pts    = (pdata['ax'].get_position()).get_points()
        aspect = (pts[1][1]-pts[0][1])/(pts[1][0]-pts[0][0])
#    elif aspect=='equal':
#        aspect = 1

    offset  = offsetpt
    if 'axoffset' in pdata: offset = pdata['axoffset']

    offsetl = np.max([offset * aspect, offset])
    offsetb = np.max([offset / aspect, offset])
    
    # Plot axes ('open' style)
    
    if pdata['theme']=='open':
        pdata['ax'].spines['left'].set_bounds(ax_yy[0], ax_yy[1])
        pdata['ax'].spines['left'].set_position(('outward', offsetl))

        pdata['ax'].spines['bottom'].set_bounds(ax_xx[0], ax_xx[1])
        pdata['ax'].spines['bottom'].set_position(('outward', offsetb))

    # Plot axes ('boxed' style)

    elif pdata['theme']=='boxed':
        for axis in ['left', 'right']: pdata['ax'].spines[axis].set_bounds(np.min(pdata['ax'].get_ylim()), np.max(pdata['ax'].get_ylim()))
        for axis in ['bottom', 'top']: pdata['ax'].spines[axis].set_bounds(np.min(pdata['ax'].get_xlim()), np.max(pdata['ax'].get_xlim()))

        #for axis in ['left', 'right']: pdata['ax'].spines[axis].set_position(('outward', offsetl))
        #for axis in ['bottom', 'top']: pdata['ax'].spines[axis].set_position(('outward', offsetb))

    else:
        if 'show' in pdata:
            if 'right' in show:
                pdata['ax'].spines['right'].set_bounds(ax_yy[0], ax_yy[1])
                pdata['ax'].spines['right'].set_position(('outward', offsetl))
            if 'left' in show:
                pdata['ax'].spines['left'].set_bounds(ax_yy[0], ax_yy[1])
                pdata['ax'].spines['left'].set_position(('outward', offsetl))
            if 'top' in show:
                pdata['ax'].spines['top'].set_bounds(ax_xx[0], ax_xx[1])
                pdata['ax'].spines['top'].set_position(('outward', offsetl))
            if 'bottom' in show:
                pdata['ax'].spines['bottom'].set_bounds(ax_xx[0], ax_xx[1])
                pdata['ax'].spines['bottom'].set_position(('outward', offsetl))

    # Fix aspect ratio

#    if aspect==1:
#        pdata['ax'].set_aspect('equal', 'datalim')


############# AUXILIARY FUNCTIONS #############

def bezier(rad, angle, t, bezier_rad):
    """
    Generate quadratic Bezier curves between pairs of points in a circle, along with a third point
    selected according to the input Bezier radius.
    """

    # Convert points to cartesian

    rad    = np.array(rad, float)
    angle  = np.array(angle, float)
    cpolar = []

    if len(rad.T)>0:
        p1x, p1y = polar2cart(rad.T[0], angle.T[0])
        p3x, p3y = polar2cart(rad.T[1], angle.T[1])
        p2x, p2y = polar2cart(np.array([bezier_rad[i] for i in range(len(rad))],float), (angle.T[0]+angle.T[1])/2.)
        
        p1 = [[p1x[i], p1y[i]] for i in range(len(p1x))]
        p2 = [[p2x[i], p2y[i]] for i in range(len(p2x))]
        p3 = [[p3x[i], p3y[i]] for i in range(len(p3x))]

        # Get Bezier curves

        curve  = [bezier_primitive(p1[i], p2[i], p3[i], t) for i in range(len(p1))]

        # Convert to polar and return

        cpolar = [[curve[i][0], curve[i][1]] for i in range(len(curve))]
    
    return cpolar


def bezier_primitive(p1, p2, p3, t):
    """
    Generate a quadratic Bezier curve from three input points, along with list of points to sample.
    """

    p1 = np.array(p1,float)
    p2 = np.array(p2,float)
    p3 = np.array(p3,float)
    t  = np.array(t,float)

    bcurve = np.outer(p1, (1 - t) * (1 - t)) + np.outer(p2, 2 * (1 - t) * t) + np.outer(p3, t * t)
    return bcurve


def rad2deg(x):
    """ Convert radians to degrees. """
    return (180. * x / np.pi)


def deg2rad(x):
    """ Convert degrees to radians. """
    return (np.pi * x / 180.)


def polar2cart(r, t):
    """ Convert polar to cartesian coordinates. """

    x = r * np.cos(t)
    y = r * np.sin(t)

    return x, y


def cart2polar(x, y):
    """ Convert cartesian to polar coordinates. """

    r = np.sqrt((x**2)+(y**2))
    t = np.arccos(x/r)

    for i in range(len(t)):
        if np.isnan(t[i]): t[i]=0

    return r, t


def adjustlim(lim, logscale):
    """ Adjust limits to slightly offset the data from the axes. """

    if logscale:
        lim[0] = np.exp( np.log(lim[0]) - (offset / (1.0 - offset)) * (np.log(lim[1]) - np.log(lim[0])) )

    else:
        tot    = lim[1]-lim[0]
        #delta  = (offset / (1.0 - offset)) * tot
        delta  = (tot / (1.0 - offset)) - tot
        lim[0] = lim[0] - delta


def getticks(lim):
    """ Return a set of linearly spaced ticks, including the upper and lower limits as major ticks. """

    tick      = np.array([lim[0], lim[1]])
    dx        = (tick[1]-tick[0])/5.0
    minortick = np.arange(tick[0]+dx, tick[1], dx)

    return tick, minortick


def getlogticks(lim, insertLow=False, insertHigh=True):
    """ Return a set of logarithmically spaced ticks, including the upper and lower limits as major ticks. """

    tick      = []
    minortick = []

    low  = np.ceil(np.log10(lim[0]))
    high = np.floor(np.log10(lim[1]))
    tick = np.logspace(low, high, num=int(1+high-low))
    
    if lim[0]!=10**low  and insertLow:  tick = np.insert(tick, 0,         lim[0])
    if lim[1]!=10**high and insertHigh: tick = np.insert(tick, len(tick), lim[1])
    
    if tick[0]!=10**low:
        minortick = [k * (tick[0]/10.0) for k in range(1,10) if (k * (tick[0])/10.0)>lim[0]]
        for i in range(1,len(tick)-1): minortick += [k * tick[i]  for k in range(2,10) if (k * tick[i])<lim[1]]
        if tick[-1]!=10**high:         minortick += [k * tick[-1] for k in range(2,10) if (k * tick[-1])<lim[1]]

    else:
        for i in range(len(tick)-1): minortick += [k * tick[i]  for k in range(2,10) if (k * tick[i])<lim[1]]
        if tick[-1]!=10**high:       minortick += [k * tick[-1] for k in range(2,10) if (k * tick[-1])<lim[1]]

    return tick, minortick




