"""
This Python file contains code for generating plots in:
    Resolving genetic linkage reveals patterns of selection in HIV-1 evolution
    bioRxiv
    Muhammad S. Sohail
    Raymond H. Y. Louie
    Matthew R. McKay (m.mckay@ust.hk)
    John P. Barton (john.barton@ucr.edu)

Additional code to pre-process the data and pass it to these plotting
routines is contained in the Jupyter notebook `figures.ipynb`.
"""

#############  PACKAGES  #############

import sys, os
from copy import deepcopy

import numpy as np

import scipy as sp
import scipy.stats as st

import pandas as pd

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plot
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg

import seaborn as sns

from colorsys import hls_to_rgb

import mplot as mp


############# PARAMETERS #############

# GLOBAL VARIABLES

NUC = ['-', 'A', 'C', 'G', 'T']
REF = NUC[0]
WFS_DIR = 'src/wfsim'
MPL_DIR = 'src/MPL'
HIV_DIR = 'data/HIV'
SIM_DIR = 'data/simulation'

# Standard color scheme

BKCOLOR  = '#252525'
LCOLOR   = '#969696'
C_BEN    = '#EB4025' #'#F16913'
C_BEN_LT = '#F08F78' #'#fdd0a2'
C_NEU    =  LCOLOR   #'#E8E8E8' # LCOLOR
C_NEU_LT = '#E8E8E8' #'#F0F0F0' #'#d9d9d9'
C_DEL    = '#3E8DCF' #'#604A7B'
C_DEL_LT = '#78B4E7' #'#dadaeb'

# Plot conventions

def cm2inch(x): return float(x)/2.54
SINGLE_COLUMN   = cm2inch(8.5)
ONE_FIVE_COLUMN = cm2inch(11.4)
DOUBLE_COLUMN   = cm2inch(17.4)

GOLDR        = (1.0 + np.sqrt(5)) / 2.0
FONTFAMILY   = 'Arial'  # paper style
SIZESUBLABEL = 8
SIZELABEL    = 6
SIZETICK     = 6
SMALLSIZEDOT = 6.
SIZELINE     = 0.6
#FONTFAMILY   = 'Arial'  # grant style
#SIZESUBLABEL = 10
#SIZELABEL    = 10
#SIZETICK     = 10
#SMALLSIZEDOT = 6. * 2
#SIZELINE     = 1
#FONTFAMILY   = 'Avenir'  # slides style
#SIZESUBLABEL = 20
#SIZELABEL    = 20
#SIZETICK     = 20
#SMALLSIZEDOT = 6. * 7
#SIZELINE     = 1.5
#SLIDE_WIDTH  = 10.5
TICKLENGTH   = 3
TICKPAD      = 3
AXWIDTH      = 0.4

FIGPROPS = {
    'transparent' : True,
    #'bbox_inches' : 'tight'
}

DEF_BARPROPS = {
    'lw'          : SIZELINE/2,
    'width'       : 0.25,
    'edgecolor'   : BKCOLOR,
    'align'       : 'center', #other option: edge
    'orientation' : 'vertical'
}

DEF_HISTPROPS = {
    'histtype'    : 'bar',
    'lw'          : SIZELINE/2,
    'rwidth'      : 0.8,
    'ls'          : 'solid',
    'edgecolor'   : BKCOLOR,
    'alpha'       : 0.5
}

DEF_ERRORPROPS = {
    'mew'        : AXWIDTH,
    'markersize' : SMALLSIZEDOT/2,
    'fmt'        : 'o',
    'elinewidth' : SIZELINE/2,
    'capthick'   : 0,
    'capsize'    : 0
}

DEF_LINEPROPS = {
    'lw' : SIZELINE,
    'ls' : '-'
}

DEF_LABELPROPS = {
    'family' : FONTFAMILY,
    'size'   : SIZELABEL,
    'color'  : BKCOLOR
}

DEF_SUBLABELPROPS = {
    'family'  : FONTFAMILY,
    'size'    : SIZESUBLABEL,
    'weight'  : 'bold',
    'ha'      : 'center',
    'va'      : 'center',
    'color'   : 'k',
    'clip_on' : False
}

DEF_TICKLABELPROPS = {
    'family' : FONTFAMILY,
    'size'   : SIZETICK,
    'color'  : BKCOLOR
}

DEF_TICKPROPS = {
    'length'    : TICKLENGTH,
    'width'     : AXWIDTH/2,
    'pad'       : TICKPAD,
    'axis'      : 'both',
    'direction' : 'out',
    'colors'    : BKCOLOR,
    'bottom'    : True,
    'left'      : True,
    'top'       : False,
    'right'     : False
}

DEF_MINORTICKPROPS = {
    'length'    : TICKLENGTH-1.25,
    'width'     : AXWIDTH/2,
    'axis'      : 'both',
    'direction' : 'out',
    'which'     : 'minor',
    'color'     : BKCOLOR
}

DEF_AXPROPS = {
    'linewidth' : AXWIDTH,
    'linestyle' : '-',
    'color'     : BKCOLOR
}

PARAMS = {'text.usetex': False, 'mathtext.fontset': 'stixsans', 'mathtext.default' : 'regular'}
plot.rcParams.update(PARAMS)


############# PLOTTING  FUNCTIONS #############

def plot_figure_example_mpl(**pdata):
    """
    Example evolutionary trajectory for a 50-site system and inferred selection coefficients,
    together with aggregate properties across sampling levels.
    """
    
    # unpack data

    n_gen  = pdata['n_gen']
    dg     = pdata['dg']
    N      = pdata['N']
    xfile  = pdata['xfile']
    method = pdata['method']

    n_ben = pdata['n_ben']
    n_neu = pdata['n_neu']
    n_del = pdata['n_del']
    s_ben = pdata['s_ben']
    s_neu = pdata['s_neu']
    s_del = pdata['s_del']
    
    r_seed = 0
    if 'r_seed' in pdata:
        r_seed = pdata['r_seed']
    np.random.seed(r_seed)

    # load and process data files

    data  = np.loadtxt('%s/data/%s.dat' % (WFS_DIR, xfile))
    times = np.unique(data.T[0])
    x     = []
    for i in range(0, n_gen, dg):
        idx    = data.T[0]==times[i]
        t_data = data[idx].T[2:].T
        t_num  = data[idx].T[1].T
        t_freq = np.einsum('i,ij->j', t_num, t_data) / float(np.sum(t_num))
        x.append(t_freq)
    x = np.array(x).T

    s_true = [s_ben for i in range(n_ben)] + [0 for i in range(n_neu)] + [s_del for i in range(n_del)]
    s_inf  = np.loadtxt('%s/out/%s_%s.dat' % (MPL_DIR, xfile.split('wfsim_')[1], method))
    cov    = np.loadtxt('%s/out/covariance-%s.dat' % (MPL_DIR, xfile.split('wfsim_')[1]))
    ds     = np.linalg.inv(cov) / N

    # PLOT FIGURE

    ## set up figure grid

    w     = SINGLE_COLUMN #SLIDE_WIDTH
    goldh = w / 1.4
    fig   = plot.figure(figsize=(w, goldh))

    box_traj = dict(left=0.17, right=0.95, bottom=0.65, top=0.95)
    box_coef = dict(left=0.17, right=0.95, bottom=0.05, top=0.50)

    gs_traj = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_traj)
    gs_coef = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_coef)
    ax_traj = plot.subplot(gs_traj[0, 0])
    ax_coef = plot.subplot(gs_coef[0, 0])
    
    dx = -0.12
    dy =  0.02

    ## a -- all trajectories together

    pprops = { 'xticks':      [0, 100, 200, 300, 400],
               'yticks':      [0, 1],
               'yminorticks': [0.25, 0.5, 0.75],
               'nudgey':      1.1,
               'xlabel':      'Generation',
               'ylabel':      'Allele\nfrequency, ' + r'$x$',
               'plotprops':   {'lw': SIZELINE, 'ls': '-', 'alpha': 0.6 },
               'axoffset':    0.1,
               'theme':       'open' }

    xdat = [range(0, n_gen, dg) for k in range(n_ben)]
    ydat = [k for k in x[:n_ben]]
    pprops['plotprops']['alpha'] = 1
    mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[C_BEN_LT for k in range(len(x))], **pprops)

    xdat = [range(0, n_gen, dg) for k in range(n_neu)]
    ydat = [k for k in x[n_ben:n_ben+n_neu]]
    pprops['plotprops']['alpha'] = 0.4
    mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[C_NEU for k in range(len(x))], **pprops)

    xdat = [range(0, n_gen, dg) for k in range(n_del)]
    ydat = [k for k in x[n_ben+n_neu:]]
    pprops['plotprops']['alpha'] = 1
    mp.plot(type='line', ax=ax_traj, x=xdat, y=ydat, colors=[C_DEL_LT for k in range(len(x))], **pprops)
    
    ax_traj.text(box_traj['left']+dx, box_traj['top']+dy-0.01, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b -- individual beneficial/neutral/deleterious selection coefficients

    sprops = { 'lw' : 0, 's' : 9., 'marker' : 'o' }

    pprops = { 'xlim':        [ -0.3,    5],
               'ylim':        [-0.04, 0.04],
               'yticks':      [-0.04, 0, 0.04],
               'yminorticks': [-0.03, -0.02, -0.01, 0.01, 0.02, 0.03],
               'yticklabels': [-4, 0, 4],
               'xticks':      [],
               'ylabel':      'Inferred selection\ncoefficient, ' + r'$\hat{s}$' + ' (%)',
               'theme':       'open',
               'hide':        ['bottom'] }

    n_coef    = [n_ben, n_neu, n_del]
    c_coef    = [C_BEN, C_NEU, C_DEL]
    c_coef_lt = [C_BEN_LT, C_NEU_LT, C_DEL_LT]
    offset    = [0, n_ben, n_ben+n_neu]

    for k in range(len(n_coef)):
        mp.line(ax=ax_coef, x=[[k-0.35, k+0.35]], y=[[s_true[offset[k]], s_true[offset[k]]]], colors=[BKCOLOR], plotprops=dict(lw=SIZELINE, ls=':'), **pprops)
        plotprops = DEF_ERRORPROPS.copy()
        plotprops['alpha'] = 1
        for i in range(n_coef[k]):
            xdat = [k + np.random.normal(0, 0.08)]
            ydat = [s_inf[offset[k]+i]]
            yerr = np.sqrt(ds[offset[k]+i][offset[k]+i])
            if i==n_coef[k]-1 and k==len(n_coef)-1:
                mp.plot(type='error', ax=ax_coef, x=[xdat], y=[ydat], yerr=[yerr], edgecolor=[c_coef[k]], facecolor=[c_coef_lt[k]], plotprops=plotprops, **pprops)
            else:
                mp.error(ax=ax_coef, x=[xdat], y=[ydat], yerr=[yerr], edgecolor=[c_coef[k]], facecolor=[c_coef_lt[k]], plotprops=plotprops, **pprops)

    coef_legend_x  =  3.8
    coef_legend_d  = -0.15
    coef_legend_dy = -0.00948646125116713 #-0.010540512501296811
    coef_legend_y  = [0.035, 0.035 + coef_legend_dy, 0.035 + (2*coef_legend_dy)]
    coef_legend_t  = ['Beneficial', 'Neutral', 'Deleterious']
    for k in range(len(coef_legend_y)):
        mp.error(ax=ax_coef, x=[[coef_legend_x+coef_legend_d]], y=[[coef_legend_y[k]]], edgecolor=[c_coef[k]], facecolor=[c_coef_lt[k]], plotprops=plotprops, **pprops)
        ax_coef.text(coef_legend_x, coef_legend_y[k], coef_legend_t[k], ha='left', va='center', **DEF_LABELPROPS)

    yy = 0.035 + (3.5 * coef_legend_dy)
    mp.line(ax=ax_coef, x=[[coef_legend_x-0.21, coef_legend_x-0.09]], y=[[yy, yy]], colors=[BKCOLOR], plotprops=dict(lw=SIZELINE, ls=':'), **pprops)
    ax_coef.text(coef_legend_x, yy, 'True selection\ncoefficient', ha='left', va='center', **DEF_LABELPROPS)

    ax_coef.text(box_coef['left']+dx, box_coef['top']+dy, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)


#    ## FORMATTING CHECK
#    # Get the x and y data and transform it into pixel coordinates
#    xy_pixels = []
#    xy_pixels.append(ax_coef.transData.transform((coef_legend_x,      coef_legend_y[0])))
#    xy_pixels.append(ax_coef.transData.transform((coef_legend_x-0.21, coef_legend_y[1])))
#    xy_pixels.append(ax_coef.transData.transform((coef_legend_x-0.09, coef_legend_y[1])))
#    print(xy_pixels)
#
#    dx1 = xy_pixels[0][0]-xy_pixels[1][0]
#    dx2 = xy_pixels[0][0]-xy_pixels[2][0]
#    dy  = xy_pixels[0][1]-xy_pixels[1][1]
#
#    print(dx1, dx2, dy)
#
#    invt = ax_coef.transData.inverted()
#    xy1  = invt.transform((0,0))
#    xy2  = invt.transform((0,9))
#
#    print(xy1[1]-xy2[1])
#    ## FORMATTING CHECK

    # SAVE FIGURE

    plot.savefig('figures/fig1-example-mpl.pdf', dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    #plot.savefig('figures/fig1.png', dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figure 1 done.')


def plot_figure_hiv_summary(**pdata):
    """
    Patterns of selection in HIV-1 across patients, including a) the classification of top beneficial mutations,
    b) enrichment in CD8+ T cell escapes, and c) enrichment in reversions.
    """

    # unpack data
    
    n_poly    = pdata['n_poly']
    x_enr     = pdata['x_enr']
    y_CD8_MPL = pdata['y_CD8_MPL']
    y_CD8_SL  = pdata['y_CD8_SL']
    y_rev_MPL = pdata['y_rev_MPL']
    y_rev_SL  = pdata['y_rev_SL']

    # PLOT FIGURE

    ## set up figure grid

    w       = ONE_FIVE_COLUMN # 1.5*SLIDE_WIDTH
    hshrink = 0.55 #0.85
    goldh   = w * hshrink
    fig     = plot.figure(figsize=(w, goldh))
    
    box_top  = 0.88 #0.93
    box_circ = dict(left=0.13, right=0.36, bottom=box_top-(0.23/hshrink)-(0.08/hshrink), top=box_top-(0.08/hshrink))
    box_epit = dict(left=0.63, right=0.93, bottom=box_top-(0.15/hshrink),                top=box_top)
    box_reve = dict(left=0.63, right=0.93, bottom=box_top-(0.15/hshrink)-(0.22/hshrink), top=box_top-(0.22/hshrink))
    
    gs_circ = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ)
    gs_epit = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_epit)
    gs_reve = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_reve)
    ax_circ = plot.subplot(gs_circ[0, 0])
    ax_epit = plot.subplot(gs_epit[0, 0])
    ax_reve = plot.subplot(gs_reve[0, 0])

    ## a -- circle plot
    
    arc_colors = ['#FF6E1B', '#FF8F28', '#FFD200', '#FFE552', '#969696', '#cccccc', '#00A3AD', '#1295D8']
    arc_r      = 1
    arc_props  = dict(center=[0,0], r=arc_r, width=arc_r, lw=AXWIDTH)

    arc_list = [  ]
    arc_degs = [90]
    curr_deg =  90
    for i in range(len(n_poly)):
        curr_deg += 360 * n_poly[i]
        arc_degs.append(curr_deg)
    
    for i in range(len(n_poly)):
        arc_list.append(matplotlib.patches.Wedge(theta1=arc_degs[i], theta2=arc_degs[i+1], fc=arc_colors[i], ec='w', **arc_props))
    
    arc_list.append(matplotlib.patches.Wedge(theta1=0, theta2=360, fc='none', ec=BKCOLOR, **arc_props))
    
    for arc in arc_list:
        ax_circ.add_artist(arc)

    dr = 0.15
    label = ['Env exposed (%.1f%%)'                     % (100*n_poly[0]),
             '±N-linked\nglycosylation\nmotif (%.1f%%)' % (100*n_poly[1]),
             'CD8+ T cell\nescape\n(%.1f%%)'            % (100*n_poly[2]),
             'Flanking\nCD8+ T cell\nepitope (%.1f%%)'  % (100*n_poly[3]),
             'Other\nsynonymous\n(%.1f%%)'              % (100*n_poly[4]),
             'Synonymous\nreversion (%.1f%%)'           % (100*n_poly[5]),
             'Other\nnonsynonymous\nreversion (%.1f%%)' % (100*n_poly[6]),
             'Other\nnonsynonymous\n(%.1f%%)'           % (100*n_poly[7])]
    for i in range(len(n_poly)):
        if i in [-1]:
            continue
        else:
            center  = np.pi * (arc_degs[i] + arc_degs[i+1]) / 360
            label_x = [(arc_r - dr) * np.cos(center), (arc_r + dr) * np.cos(center)]
            label_y = [(arc_r - dr) * np.sin(center), (arc_r + dr) * np.sin(center)]
            off_x   = 0.05 * np.cos(center)
            off_y   = 0.05 * np.sin(center)

            ddx       = 0.03
            ddy       = 0.10
            txtprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
            plotprops = dict(lw=AXWIDTH/2, ls='-', clip_on=False)
            
            if label_x[0]<0:
                txtprops['ha'] = 'right'
            else:
                txtprops['ha'] = 'left'

            if label_y[0]<-0.7:
                txtprops['ha'] = 'center'
                txtprops['va'] = 'top'
            if label_y[0]>0.80:
                txtprops['ha'] = 'center'
                txtprops['va'] = 'bottom'

            # 1% positions
            if 'Env' in label[i]:
                off_x += 0.80
            if 'glycosylation' in label[i]:
                txtprops['ha'] = 'right'
                off_x -= 0.05
                off_y -= 0.22
                label_x.append(label_x[-1]-dr)
                label_y.append(label_y[-1])
            if 'Flanking' in label[i]:
                off_x -= 0.64
                off_y += 0.45
                label_x.append(label_x[-1]-dr)
                label_y.append(label_y[-1])
#            if 'Other\nsynonymous' in label[i]:
#                txtprops['ha'] = 'right'
#                off_x += 0.08
#                off_y -= 0.05
            if 'Other\nsynonymous' in label[i]:
                txtprops['ha'] = 'right'
                off_x -= 0.06
                off_y += 0.22
                label_x.append(label_x[-1]-dr)
                label_y.append(label_y[-1])
            if 'Synonymous\nreversion' in label[i]:
                txtprops['ha'] = 'left'
                off_x -= 0.13
                off_y -= 0.05
            if 'Other\nnonsynonymous\nreversion' in label[i]:
                txtprops['ha'] = 'center'
                off_x += 0.78
                off_y += 0.45
                label_x.append(label_x[-1]+dr)
                label_y.append(label_y[-1])

#            if i==0:
#                label_x[-1] += 0.85 * dr * np.cos(center)
#                label_y[-1] += 0.85 * dr * np.sin(center)

            if 'Flanking' in label[i]:
                continue

            ax_circ.text(label_x[-1]+off_x, label_y[-1]+off_y, label[i], **txtprops)

            if i==len(n_poly)-1:
                pprops = dict(xlim=[-1.1, 1.1], ylim=[-1.1, 1.1], xticks=[], yticks=[], noaxes=True)
                mp.plot(type='line', ax=ax_circ, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops, **pprops)
            else:
                mp.line(ax=ax_circ, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)

    dx = -0.08
    dy =  0.025/hshrink
    ax_circ.text(box_circ['left']+dx, box_epit['top']+dy, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b -- enrichment in CD8+ T cell escape mutations
    
    lineprops = { 'lw': SIZELINE*1.2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }
    fillprops = { 'lw': 0, 'alpha': 0.2, 'interpolate': True, 'step': 'mid' }
    
    pprops = { 'xlim':        [0.01,  0.199],
               'ylim':        [0, 20],
               'yticks':      [0, 20],
               'yminorticks': [5, 10, 15],
               'xticklabels': [1, 10],
               'ylabel':      'Fold enrichment\nin CD8+ T cell\nescape mutations',
               'logx':        True,
               'theme':       'open' }

    C_MPL = '#FFB511'
    C_SL  = C_NEU
    
    pprops['colors'] = [C_MPL]
    mp.line(ax=ax_epit, x=[x_enr], y=[y_CD8_MPL], plotprops=lineprops, **pprops)
    mp.fill(ax=ax_epit, x=[x_enr], y=[y_CD8_MPL], plotprops=fillprops, **pprops)

    pprops['colors'] = [C_SL]
    mp.line(             ax=ax_epit, x=[x_enr], y=[y_CD8_SL], plotprops=lineprops, **pprops)
    mp.plot(type='fill', ax=ax_epit, x=[x_enr], y=[y_CD8_SL], plotprops=fillprops, **pprops)

    dx = -0.10
    dy =  0.025/hshrink
    ax_epit.text(box_epit['left']+dx, box_epit['top']+dy, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## c -- enrichment in reversions

    pprops = { 'xlim':        [0.01,  0.199],
               'ylim':        [0, 40],
               'yticks':      [0, 40],
               'yminorticks': [10, 20, 30],
               'xticklabels': [1, 10],
               'xlabel':      'Fraction of most positively\nselected variants (%)',
               'ylabel':      'Fold enrichment\nin reversions\noutside epitopes',
               'logx':        True,
               'theme':       'open' }

    C_MPL = '#FFB511'
    C_SL  = C_NEU

    pprops['colors'] = [C_MPL]
    mp.line(ax=ax_reve, x=[x_enr], y=[y_rev_MPL], plotprops=lineprops, **pprops)
    mp.fill(ax=ax_reve, x=[x_enr], y=[y_rev_MPL], plotprops=fillprops, **pprops)

    pprops['colors'] = [C_SL]
    mp.line(             ax=ax_reve, x=[x_enr], y=[y_rev_SL], plotprops=lineprops, **pprops)
    mp.plot(type='fill', ax=ax_reve, x=[x_enr], y=[y_rev_SL], plotprops=fillprops, **pprops)

    ax_reve.text(box_reve['left']+dx, box_reve['top']+dy, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # legend

    invt = ax_circ.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = xy1[0]-xy2[0]
    legend_dx2 = xy1[0]-xy3[0]
    legend_dy  = xy1[1]-xy2[1]

    enr_legend_x = 5.60 #1.75
    enr_legend_y = 1.90 #1.65
    enr_legend_t = ['MPL', 'Independent\nmodel']
    enr_legend_c = [C_MPL, C_SL]
    for k in range(len(enr_legend_t)):
        mp.line(ax=ax_circ, x=[[enr_legend_x + legend_dx1, enr_legend_x + legend_dx2]],
                y=[[enr_legend_y + (1.5 * k * legend_dy), enr_legend_y + (1.5 * k * legend_dy)]],
                colors=[enr_legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_circ.text(enr_legend_x, enr_legend_y + (1.5 * k * legend_dy), enr_legend_t[k], ha='left', va='center', **DEF_LABELPROPS)

    # SAVE FIGURE
    
    plot.savefig('figures/fig2-hiv-summary.pdf', dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figure 2 done.')


def plot_figure_ch77_kf9(**pdata):
    """
    CH77 KF9 epitope escape mutation frequencies, inferred selection coefficients, and linkage.
    """

    # unpack data
    
    patient       = pdata['patient']
    region        = pdata['region']
    inf_idxs      = pdata['inf_idxs']
    epitope       = pdata['epitope']
    epitope_range = pdata['epitope_range']
    epitope_label = pdata['epitope_label']
    cov_label     = pdata['cov_label']
    label2ddr     = pdata['label2ddr']
    legend_loc    = pdata['legend_loc']
    traj_ticks    = pdata['traj_ticks']
    sel_ticks     = pdata['sel_ticks']
    sel_minors    = pdata['sel_minors']
    sel_space     = pdata['sel_space']
    fig_title     = pdata['fig_title']
    tag           = patient+'-'+region
    
    # process stored data
    
    df_poly = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_esc  = df_poly[(df_poly.epitope==epitope) & (df_poly.nucleotide!=df_poly.TF)]

    times = [int(i.split('_')[-1]) for i in df_esc.columns if 'f_at_' in i]
    times.sort()
    
    var_tag  = []
    var_smpl = []
    var_sind = []
    var_traj = []
    var_idxs = []
    curr_HXB2 = 8988
    curr_aln  = 4015
    curr_char = 'a'
    for df_iter, df_entry in df_esc.iterrows():
        if df_entry.nucleotide=='-':
            continue
        if pd.notnull(df_entry.HXB2_index):
            var_tag.append(str(int(df_entry.HXB2_index))+df_entry.nucleotide)
            var_idxs.append([int(df_entry.polymorphic_index), str(df_entry.nucleotide)])
            curr_HXB2 = int(df_entry.HXB2_index)
            curr_aln  = int(df_entry.alignment_index)
        else:
            if int(df_entry.alignment_index)!=curr_aln:
                curr_aln  = int(df_entry.alignment_index)
                curr_char = chr(ord(curr_char) + 1)
            var_tag.append(str(int(curr_HXB2))+curr_char+df_entry.nucleotide)
            var_idxs.append([int(df_entry.polymorphic_index), str(df_entry.nucleotide)])
        var_traj.append([df_entry['f_at_%d' % t] for t in times])
        var_smpl.append(df_entry.s_MPL)
        var_sind.append(df_entry.s_SL)

    var_c = sns.husl_palette(len(var_traj))

    df_ds = pd.read_csv('%s/analysis/%s-delta-s.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_ds = df_ds[~((df_ds.mask_polymorphic_index==df_ds.target_polymorphic_index)
                     & (df_ds.mask_nucleotide==df_ds.target_nucleotide))]

    ds_matrix = []
    for i in range(len(var_idxs)):
        ds_vec = []
        sub_df_ds = df_ds[(df_ds.target_polymorphic_index==var_idxs[i][0]) & (df_ds.target_nucleotide==var_idxs[i][1])]
        for j in range(len(var_idxs)):
            if j==i:
                ds_vec.append(0)
            else:
                ds_vec.append(sub_df_ds[(sub_df_ds.mask_polymorphic_index==var_idxs[j][0]) & (sub_df_ds.mask_nucleotide==var_idxs[j][1])].iloc[0].effect)
        for j in range(len(inf_idxs)):
            if var_idxs[i][0]==inf_idxs[j][0] and var_idxs[i][1]==inf_idxs[j][1]:
                ds_vec.append(0)
            else:
                ds_vec.append(sub_df_ds[(sub_df_ds.mask_polymorphic_index==inf_idxs[j][0]) & (sub_df_ds.mask_nucleotide==inf_idxs[j][1])].iloc[0].effect)
        ds_matrix.append(ds_vec)

    inf_tag = []
    for i in range(len(inf_idxs)):
        temp_df = df_poly[(df_poly.polymorphic_index==inf_idxs[i][0]) & (df_poly.nucleotide==inf_idxs[i][1])]
        if pd.notnull(temp_df.iloc[0].HXB2_index):
            inf_tag.append(str(int(temp_df.iloc[0].HXB2_index))+str(temp_df.iloc[0].nucleotide))
        else:
            inf_tag.append('8865yA') # hard-coded insertion in DG9 epitope

    for i in range(len(ds_matrix)):
        print('%s' % ('\t'.join(['%.4f' % ds for ds in ds_matrix[i]])))
    print(np.max(ds_matrix))

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN # SLIDE_WIDTH
    hshrink = 0.5 #0.77
    goldh   = 1.0 * w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top  = 0.91
    dy       = 0.14/hshrink
    y0       = 0.08/hshrink
    y1       = (0.10/hshrink) * (0.205 / 0.225)
    y2       = 0.42/hshrink
    box_traj = dict(left=0.16/2, right=0.40/2, bottom=box_top-y0, top=box_top)
    box_smpl = dict(left=0.61/2, right=(0.61+(0.032*len(var_c)))/2, bottom=box_top-y0, top=box_top)
    box_sind = dict(left=(0.61+sel_space+(0.032*len(var_c)))/2, right=(0.61+sel_space+(2*0.032*len(var_c)))/2, bottom=box_top-y0, top=box_top)
    box_ds   = dict(left=0.16/2, right=0.57/2, bottom=box_top-y0-y1-dy, top=box_top-y0-dy)
    box_circ = dict(left=0.54, right=0.96, bottom=box_top-y2, top=box_top)

    gs_traj = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_traj)
    gs_smpl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_smpl)
    gs_sind = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_sind)
    gs_ds   = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_ds)
    gs_circ = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ)
    ax_traj = plot.subplot(gs_traj[0, 0])
    ax_smpl = plot.subplot(gs_smpl[0, 0])
    ax_sind = plot.subplot(gs_sind[0, 0])
    ax_ds   = plot.subplot(  gs_ds[0, 0])
    ax_circ = plot.subplot(gs_circ[0, 0])

    dx = -0.06
    dy =  0.025/hshrink

    ## a -- trajectory plot

    pprops = { 'xticks':      traj_ticks,
               'yticks':      [0, 1],
               'yminorticks': [0.25, 0.5, 0.75],
               'nudgey':      1.1,
               'xlabel':      'Time (days)',
               'ylabel':      'Variant frequency\nin %s epitope\n' % cov_label,
               'plotprops':   {'lw': SIZELINE*1.5, 'ls': '-', 'alpha': 1.0 },
               'axoffset':    0.1,
               'theme':       'open' }

    for i in range(len(var_tag)-1):
        xdat = [times]
        ydat = [var_traj[i]]
        mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[var_c[i]], **pprops)

    xdat = [times]
    ydat = [var_traj[len(var_tag)-1]]
    mp.plot(type='line', ax=ax_traj, x=xdat, y=ydat, colors=[var_c[len(var_tag)-1]], **pprops)

    ax_traj.text(box_traj['left']+dx, box_traj['top']+dy, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b1 -- selection coefficients inferred by MPL

    xmin, xmax =  len(var_tag), 0
    ymin, ymax =  0, 0.08

    hist_props = dict(lw=SIZELINE/2, width=0.6, align='center', orientation='vertical',
                      edgecolor=[BKCOLOR for i in range(len(var_tag))])

    bar_x  = [i+0.5 for i in range(len(var_tag))]
    pprops = { 'colors':      [var_c],
               'xlim':        [xmax, xmin],
               'ylim':        [np.min(sel_ticks), np.max(sel_ticks)],
               'xticks':      bar_x,
               'xticklabels': var_tag,
               'yticks':      sel_ticks,
               'yticklabels': [int(100*l) for l in sel_ticks],
               'yminorticks': sel_minors,
               'ylabel':      'Inferred selection\ncoefficient, '+r'$\hat{s}$'+' (%)',
               'theme':       'open',
               'hide' :       [] }

    mp.plot(type='bar', ax=ax_smpl, x=[bar_x], y=[var_smpl], plotprops=hist_props, **pprops)
    plot.setp(ax_smpl.xaxis.get_majorticklabels(), rotation=90)

    transFigureInv = fig.transFigure.inverted()
    labelprops     = dict(color=BKCOLOR, ha='center', va='top', family=FONTFAMILY, size=SIZELABEL,
                          clip_on=False, transform=fig.transFigure)
    ax_smpl.text(box_smpl['left']+(box_smpl['right']-box_smpl['left'])/2, box_smpl['top'], 'MPL', **labelprops)

    ax_smpl.text(box_smpl['left']+dx, box_smpl['top']+dy, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b2 -- selection coefficients inferred by SL

    pprops = { 'colors':      [var_c],
               'xlim':        [xmax, xmin],
               'ylim':        [np.min(sel_ticks), np.max(sel_ticks)],
               'xticks':      bar_x,
               'xticklabels': var_tag,
               'yticks':      [],
               'theme':       'open',
               'hide' :       ['left','right'] }

    mp.plot(type='bar', ax=ax_sind, x=[bar_x], y=[var_sind], plotprops=hist_props, **pprops)
    plot.setp(ax_sind.xaxis.get_majorticklabels(), rotation=90)

    ax_sind.text(box_sind['left']+(box_sind['right']-box_sind['left'])/2, box_sind['top'], 'Independent\nmodel', **labelprops)

    # add background

    cBG = '#F5F5F5'
    bg  = ax_sind.axis()
    ddx = 0.01
    ddy = 0.01
    rec = matplotlib.patches.Rectangle(xy=(box_sind['left']-(0.1*ddx),
                                           box_sind['bottom']-(0.2*ddy)), #box_sind['bottom']-(0.7*ddy)),
                                       width=box_sind['right']-box_sind['left']+(0.2*ddx),
                                       height=box_sind['top']-box_sind['bottom']+(1.7*ddy), transform=fig.transFigure, ec=None, fc=cBG, clip_on=False, zorder=-100)
    rec = ax_sind.add_patch(rec)

    ## c -- effects of linkage on selection for sites within epitope

    site_rec_props = dict(height=1, width=1, ec=None, lw=AXWIDTH/2, clip_on=False)
    rec_patches    = []
    rec_space      = 0.0

    for i in range(len(ds_matrix)):
        for j in range(len(ds_matrix[i])):
            temp_ds = ds_matrix[i][j]
            t       = temp_ds / 0.02
            if np.fabs(t)>1:
                t /= np.fabs(t)
            if t>0:
                c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
            else:
                c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
            rec = matplotlib.patches.Rectangle(xy=((1 + rec_space) * j + (j>=len(ds_matrix)) * 0.3, len(ds_matrix) - i), fc=c, **site_rec_props)
            rec_patches.append(rec)

    epitope_rec_props = dict(height=len(ds_matrix), width=(1 + rec_space)*len(ds_matrix)-rec_space, ec=BKCOLOR,
                             fc='none', lw=AXWIDTH/2, clip_on=False)
    outside_rec_props = dict(height=len(ds_matrix), width=(1 + rec_space)*(len(ds_matrix[0]) - len(ds_matrix))-rec_space,
                             ec=BKCOLOR, fc='none', lw=AXWIDTH/2, clip_on=False)
    rec_patches.append(matplotlib.patches.Rectangle(xy=(0, 1), **epitope_rec_props))
    rec_patches.append(matplotlib.patches.Rectangle(xy=((1 + rec_space)*len(ds_matrix)+0.3, 1), **outside_rec_props))

    txtprops = dict(ha='right', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    for i in range(len(ds_matrix)):
        ax_ds.text(-0.2, len(ds_matrix) - i + 0.5, '%s' % var_tag[i], **txtprops)
    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=90)
    for i in range(len(var_tag)):
        ax_ds.text((1 + rec_space) * (i + 0.5), 0.5, '%s' % var_tag[i], **txtprops)
    for i in range(len(inf_tag)):
        ax_ds.text((1 + rec_space) * (len(var_tag) + i + 0.5) + 0.3, 0.5, '%s' % inf_tag[i], **txtprops)

    txtprops = dict(ha='center', va='bottom', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    ax_ds.text((1 + rec_space) * (2.0 + 0.5) + 0.00,  6.3, 'KF9 escape\nmutations', **txtprops)
    ax_ds.text((1 + rec_space) * (6.5 + 0.5) + 0.30,  6.3, 'Most influential\nvariants outside\nKF9 epitope', **txtprops)
    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    ax_ds.text((1 + rec_space) * (12.5 + 0.5) + 0.9, 4.0, ('Effect of variant i on inferred\nselection coefficient '
                                                            + r'$\hat{s}_{j}$' + ', ' + r'$\Delta\hat{s}_{ij}$' + ' (%)'), **txtprops)
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    ax_ds.text((1 + rec_space) * (4.0 + 0.5) + 0.15, -2, 'Variant i', **txtprops)
    ax_ds.text((1 + rec_space) * ( 9.5 + 0.5) + 0.9, 3.0+1.5, -2, **txtprops)
    ax_ds.text((1 + rec_space) * (12.5 + 0.5) + 0.9, 3.0+1.5,  0, **txtprops)
    ax_ds.text((1 + rec_space) * (15.5 + 0.5) + 0.9, 3.0+1.5,  2, **txtprops)
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=90)
    ax_ds.text(-2.70, 3, 'Target variant j', **txtprops)

    for i in range(-3, 3+1, 1):
        c = ''
        t = i/3
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=((1 + rec_space) * (12.5 + i) + 0.9, 3.5+1.5), fc=c, **site_rec_props)
        rec_patches.append(rec)

    for patch in rec_patches:
        ax_ds.add_artist(patch)

    pprops = { 'colors': [BKCOLOR],
               'xlim': [0, len(ds_matrix[0]) + 1],
               'ylim': [1, len(ds_matrix)],
               'xticks': [],
               'yticks': [],
               'plotprops': dict(lw=0, s=0.2*SMALLSIZEDOT, marker='o', clip_on=False),
               'ylabel': '',
               'theme': 'open',
               'hide' : ['top', 'bottom', 'left', 'right'] }

    mp.plot(type='scatter', ax=ax_ds, x=[[-5]], y=[[-5]], **pprops)

    ax_ds.text(box_traj['left']+dx, box_ds['top']+0.03/hshrink, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## d -- circle plot

    sig_s, sig_site_real, sig_nuc_idx, epitope_start, epitope_end = plot_circle(ax_circ, tag, epitope_range, epitope_label, cov_label, label2ddr)
    
    ax_circ.text(box_circ['left']-0.02, box_circ['top']+dy, 'd'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # MAKE LEGEND

    invt = ax_circ.transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((0,9))
    
    coef_legend_x  = -1.40
    coef_legend_dx = 0.05
    coef_legend_y  = -0.89
    coef_legend_dy = (xy1[1]-xy2[1])
    
    if legend_loc=='top':
        coef_legend_x = 0.65 # 0.55
        coef_legend_y = 1.45 # 1.60

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    s_mult   = 20*SMALLSIZEDOT
    ex_s     = [-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]
    show_s   = [   1,     0,     0,     0,     0, 1,    0,    0,    0,    0,    1]
    c_s      = [C_DEL, C_BEN]
    for i in range(len(ex_s)):
        plotprops = dict(lw=0, marker='o', s=np.fabs(ex_s[i])*s_mult, clip_on=False)
        mp.scatter(ax=ax_circ, x=[[coef_legend_x + i*coef_legend_dx]], y=[[coef_legend_y]], colors=[c_s[ex_s[i]>0]], plotprops=plotprops)
        if show_s[i]:
            ax_circ.text(coef_legend_x + i*coef_legend_dx, coef_legend_y + 0.75*coef_legend_dy, '%d' % (100*ex_s[i]), **txtprops)
    ax_circ.text(coef_legend_x + 5*coef_legend_dx, coef_legend_y + (2.25*coef_legend_dy), 'Inferred selection\ncoefficient, $\hat{s}$ (%)', **txtprops)

    if legend_loc=='top':
        coef_legend_x = 0.80 # 0.65
    coef_legend_y = coef_legend_y - (1.5*coef_legend_dy)
    txtprops = dict(ha='left', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    arc_props = dict(lw=AXWIDTH/2, ec=BKCOLOR, fc='#f2f2f2', alpha=0.8, clip_on=False)
    epatch    = matplotlib.patches.Rectangle(xy=(coef_legend_x-0.025, coef_legend_y+(0.0*coef_legend_dy)), width=0.05, height=-1.*coef_legend_dy, **arc_props)
    ax_circ.add_artist(epatch)
    ax_circ.text(coef_legend_x+(1*coef_legend_dx), coef_legend_y-(coef_legend_dy/2), 'CD8+ T cell\nepitope', **txtprops)

    # SAVE FIGURE

    plot.savefig('figures/%s.pdf' % fig_title, dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figure 3 done.')


def plot_figure_cap256_vrc26(**pdata):
    """
    CAP256 mutation frequencies, inferred selection coefficients, and linkage.
    """

    # unpack data
    
    patient = pdata['patient']
    region  = pdata['region']
    tag     = patient+'-'+region

    # process stored data
    
    df_index = pd.read_csv('%s/processed/%s-SU-%s-index.csv' % (HIV_DIR, patient, region), comment='#', memory_map=True)
    df_poly  = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_sub   = df_poly[(df_poly.nucleotide!=df_poly.TF) & (df_poly.HXB2_index>=6702) & (df_poly.HXB2_index<=6737)]

    times = [int(i.split('_')[-1]) for i in df_sub.columns if 'f_at_' in i]
    times.sort()
    
    c_vals   = ['#FFB511', C_NEU]
    var_c    = []
    var_tag  = []
    var_smpl = []
    var_sind = []
    var_traj = []
    for df_iter, df_entry in df_sub.iterrows():
        var_tag.append(str(int(df_entry.polymorphic_index))+df_entry.nucleotide)
        var_traj.append([df_entry['f_at_%d' % t] for t in times])
        var_smpl.append(df_entry.s_MPL)
        var_sind.append(df_entry.s_SL)
        df_temp = df_index[df_index.alignment==df_entry.alignment_index].iloc[0]
        if df_temp.SU==df_entry.nucleotide:
            var_c.append(c_vals[0])
        else:
            var_c.append(c_vals[-1])

    # print fraction synonymous
    df = df_poly[df_poly.nucleotide!=df_poly.TF]
    top1         = int(np.round(0.01*len(df)))
    s_MPL_sorted = np.argsort(df.s_MPL)[::-1]
    s_SL_sorted  = np.argsort(df.s_SL)[::-1]
    print('synonymous:\t%d (MPL)\t%d (SL)\t of %d total' %
          (np.sum(df.iloc[s_MPL_sorted[:top1]].nonsynonymous==0), np.sum(df.iloc[s_SL_sorted[:top1]].nonsynonymous==0), top1))
          
    # print special selection coefficients
    print('variant\ts_MPL')
    print('6709C\t%.3f' % (df_poly[(df_poly.HXB2_index==6709) & (df_poly.nucleotide=='C')].s_MPL))
    print('6717T\t%.3f' % (df_poly[(df_poly.HXB2_index==6717) & (df_poly.nucleotide=='T')].s_MPL))
    print('6730C\t%.3f' % (df_poly[(df_poly.HXB2_index==6730) & (df_poly.nucleotide=='C')].s_MPL))
    print('6730G\t%.3f' % (df_poly[(df_poly.HXB2_index==6730) & (df_poly.nucleotide=='G')].s_MPL))
    print('6730T\t%.3f' % (df_poly[(df_poly.HXB2_index==6730) & (df_poly.nucleotide=='T')].s_MPL))
    print('6731T\t%.3f' % (df_poly[(df_poly.HXB2_index==6731) & (df_poly.nucleotide=='T')].s_MPL))
    print('')

    # PLOT FIGURE

    ## set up figure grid

    w       = SINGLE_COLUMN # SLIDE_WIDTH
    hshrink = 0.91
    goldh   = 1.9 * w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top  = 0.95
    dy       = 0.11 / 2
    y0       = 0.10 / hshrink
    y1       = 0.4444444444444445 / hshrink
    y2       = 0.88 * (2.2 / (36 + 11 * 0.3)) / hshrink

    box_traj = dict(left=0.16, right=0.94, bottom=box_top-y0, top=box_top)
    box_circ = dict(left=0.08, right=0.92, bottom=box_top-y0-y1-(1.2*dy), top=box_top-y0-(1.2*dy))
    box_smpl = dict(left=0.06, right=0.94, bottom=box_top-y0-y1-y2-(3.05*dy), top=box_top-y0-y1-(3.05*dy))

    gs_traj = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_traj)
    gs_circ = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ)
    gs_smpl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_smpl)

    ax_traj = plot.subplot(gs_traj[0, 0])
    ax_circ = plot.subplot(gs_circ[0, 0])
    ax_smpl = plot.subplot(gs_smpl[0, 0])

    #dx = -0.08
    dy =  0.02/hshrink

    ## a -- trajectory plot

    pprops = { 'xticks':      [30, 190, 350, 510, 670],
               'yticks':      [0, 1],
               'yminorticks': [0.25, 0.5, 0.75],
               'nudgey':      1.1,
               'xlabel':      'Time (days)',
               'ylabel':      'Variant frequency\nin VRC26 epitope',
               'plotprops':   {'lw': 1.5*SIZELINE, 'ls': '-', 'alpha': 0.75 },
               'axoffset':    0.1,
               'theme':       'open' }

    for i in range(len(var_tag)-1):
        xdat = [times]
        ydat = [var_traj[i]]
        mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[var_c[i]], **pprops)

    xdat = [times]
    ydat = [var_traj[len(var_tag)-1]]
    mp.plot(type='line', ax=ax_traj, x=xdat, y=ydat, colors=[var_c[len(var_tag)-1]], **pprops)

    ax_traj.text(0.06, box_traj['top']+dy, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    invt = ax_traj.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1*(xy1[1]-xy2[1])

    traj_legend_x =  80
    traj_legend_y =  0.95
    traj_legend_d = -0.05
    traj_legend_t = ['Superinfecting\nvariant', 'Other']
    traj_legend_c = ['#FFB511', C_NEU]
    for k in range(len(traj_legend_t)):
        mp.line(ax=ax_traj, x=[[traj_legend_x + legend_dx1, traj_legend_x + legend_dx2]],
                y=[[traj_legend_y + (1.5 * k * legend_dy), traj_legend_y + (1.5 * k * legend_dy)]],
                colors=[traj_legend_c[k]], plotprops=dict(lw=1.5*SIZELINE, ls='-', clip_on=False))
        ax_traj.text(traj_legend_x, traj_legend_y + (1.5 * k * legend_dy), traj_legend_t[k], ha='left', va='center', **DEF_LABELPROPS)

    ## b -- circle plot

    epitope_range = [[6702, 6737]]
    epitope_label = ['VRC26 epitope']
    cov_label     = 'VRC26 epitope'
    label2ddr     = { 'tat exon 2': 0.13,
                      'rev exon 2': 0.10,
                      'DG9':        0.06,
                      'KF9':        0.06,
                      'DI9':        0.13  }

#    # _ SLIDES
#    w       = 8
#    fig     = plot.figure(figsize=(w, w))
#    ax_circ = plot.subplot(111)
#    # ^ SLIDES

    sig_s, sig_site_real, sig_nuc_idx, epitope_start, epitope_end = plot_circle(ax_circ, tag, epitope_range, epitope_label, cov_label, label2ddr)

#    # _ SLIDES
#    plot.savefig('figures/new-slides-cap256-circle.pdf', dpi=1000, **FIGPROPS)
#    plot.savefig('figures/new-slides-cap256-circle.png', dpi=1000, **FIGPROPS)
#    plot.close(fig)
#    w       = 8
#    fig     = plot.figure(figsize=(w, 2.85))
#    ax_smpl = plot.subplot(311)
#    # ^ SLIDES

    dx =  0.04
    dy = -0.02
    ax_circ.text(0.06, box_circ['top']+dy, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## c -- selection in the VRC26 epitope

    site_rec_props  = dict(height=1, width=1, ec=None, lw=AXWIDTH/2, clip_on=False)
    codon_rec_props = dict(height=4, width=3, ec=BKCOLOR, fc='none', lw=AXWIDTH/2, clip_on=False)
    cBG             = '#F5F5F5'
    rec_patches     = []
    TF_dots_x       = [ 0.5]
    TF_dots_y       = [-3.5]

    sig_s         = np.array(sig_s)
    sig_site_real = np.array(sig_site_real)
    sig_nuc_idx   = np.array(sig_nuc_idx)
    eidx          = 0
    sub_box       = 0
    sub_dx        = 3 + 0.3
    for i in range(epitope_start[eidx], epitope_end[eidx]+1, 3):
        for sub_i in range(3):
            TF_dots_x.append(sub_box*sub_dx + sub_i + 0.5)
            TF_dots_y.append(4-NUC.index(df_index.iloc[i+sub_i].TF)+0.5)
            idxs     = sig_site_real==i+sub_i
            temp_s   = sig_s[idxs]
            temp_nuc = sig_nuc_idx[idxs]
            for j in range(len(NUC)-1):
                if df_index.iloc[i+sub_i].TF==NUC[j+1]:
                    continue
                c = cBG
                if j in temp_nuc:
                    t = temp_s[list(temp_nuc).index(j)] / 0.05
                    if np.fabs(t)>1:
                        t /= np.fabs(t)
                    if t>0:
                        c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
                    else:
                        c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
                rec = matplotlib.patches.Rectangle(xy=(sub_box*sub_dx + sub_i, 3-j), fc=c, **site_rec_props)
                rec_patches.append(rec)
        rec = matplotlib.patches.Rectangle(xy=(sub_box*sub_dx, 0), **codon_rec_props)
        rec_patches.append(rec)
        sub_box += 1

    for i in range(-5, 5+1, 1):
        c = cBG
        t = i/5
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=(sub_box*sub_dx - 6.5 + i, -4), fc=c, **site_rec_props)
        rec_patches.append(rec)

    invt = ax_smpl.transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((0,9))
    coef_legend_dy = (xy1[1]-xy2[1]) # multiply by 3 for slides/poster
    c   = cBG
    rec = matplotlib.patches.Rectangle(xy=(0, -4 + 4*coef_legend_dy), fc=c, **site_rec_props) # paper
#    rec = matplotlib.patches.Rectangle(xy=(0, -4 + 6*coef_legend_dy), fc=c, **site_rec_props) # slides
    rec_patches.append(rec)

    for patch in rec_patches:
        ax_smpl.add_artist(patch)

    pprops = { 'colors': [BKCOLOR],
               'xlim': [0, epitope_end[eidx]-epitope_start[eidx] + (sub_box * (sub_dx - 3)) + 1],
               'ylim': [0, 4.1],
               'xticks': [],
               'yticks': [],
               'plotprops': dict(lw=0, s=0.2*SMALLSIZEDOT, marker='o', clip_on=False),
               'ylabel': '',
               'theme': 'open',
               'hide' : ['top', 'bottom', 'left', 'right'] }

    mp.plot(type='scatter', ax=ax_smpl, x=[TF_dots_x], y=[TF_dots_y], **pprops)

    ax_smpl.text(0.06, box_smpl['top']+0.02, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    for i in range(len(NUC)-1):
        ax_smpl.text(-0.85, 3-i+0.5, NUC[i+1], clip_on=False, **txtprops)

    txtprops['ha'] = 'left'
    ax_smpl.text(1.3, -3.5, 'TF nucleotide', clip_on=False, **txtprops)
    ax_smpl.text(1.3, -3.5 + 4*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # paper
#    ax_smpl.text(1.3, -3.5 + 6*coef_legend_dy, 'Not observed', clip_on=False, **txtprops) # slides

    txtprops['ha'] = 'center'
    txtprops['va'] = 'top'
    for i in range((epitope_end[eidx]-epitope_start[eidx]+1)//3):
        ax_smpl.text(1.5+i*sub_dx, -0.5, 160+i, clip_on=False, **txtprops)

    ax_smpl.text(sub_box*sub_dx - 11, -4.5, -5, clip_on=False, **txtprops)
    ax_smpl.text(sub_box*sub_dx -  6, -4.5,  0, clip_on=False, **txtprops)
    ax_smpl.text(sub_box*sub_dx -  1, -4.5,  5, clip_on=False, **txtprops)
    ax_smpl.text(sub_box*sub_dx -  5.5, -6.0, 'Inferred selection\ncoefficient, $\hat{s}$ (%)', clip_on=False, **txtprops)

#    # _ SLIDES
#    plot.savefig('figures/new-slides-cap256-selection.pdf', dpi=1000, **FIGPROPS)
#    plot.savefig('figures/new-slides-cap256-selection.png', dpi=1000, **FIGPROPS)
#    plot.close(fig)
#    # ^ SLIDES

    # MAKE LEGEND

    invt = ax_circ.transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((0,9))

    coef_legend_x  = -1.00
    coef_legend_dx = 0.05
    coef_legend_y  = -1.05
    coef_legend_dy = xy1[1]-xy2[1]

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    s_mult   = 20*SMALLSIZEDOT
    ex_s     = [-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]
    show_s   = [   1,     0,     0,     0,     0, 1,    0,    0,    0,    0,    1]
    c_s      = [C_DEL, C_BEN]
    for i in range(len(ex_s)):
        plotprops = dict(lw=0, marker='o', s=np.fabs(ex_s[i])*s_mult, clip_on=False)
        mp.scatter(ax=ax_circ, x=[[coef_legend_x + i*coef_legend_dx]], y=[[coef_legend_y]], colors=[c_s[ex_s[i]>0]], plotprops=plotprops)
        if show_s[i]:
            ax_circ.text(coef_legend_x + i*coef_legend_dx, coef_legend_y + 0.75*coef_legend_dy, '%d' % (100*ex_s[i]), **txtprops)
    ax_circ.text(coef_legend_x + 5*coef_legend_dx, coef_legend_y + (2.25*coef_legend_dy), 'Inferred selection\ncoefficient, $\hat{s}$ (%)', **txtprops)

    # SAVE FIGURE

    plot.savefig('figures/fig4-cap256-vrc26.pdf', dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figure 4 done.')


def plot_circle(ax, tag, epitope_range, epitope_label, cov_label, label2ddr):
    """
    Create a circle plot showing selection, large elements of the inverse covariance, and HIV sequence features.
    """

    df_info  = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_index = pd.read_csv('%s/processed/%s-index.csv' % (HIV_DIR, tag), comment='#', memory_map=True)

    cov  = np.loadtxt('%s/HIV/covariance-%s-poly-seq2state.dat' % (MPL_DIR, tag))
    num  = np.loadtxt('%s/HIV/numerator-%s-poly-seq2state.dat' % (MPL_DIR, tag))
    cinv = np.linalg.inv(cov)
    ds   = cinv / 1e4

    # Check for selection coefficients that are significant: |s| > mult * sigma_s

    mult  = 0.0
    n_sig = 0

    sig_s         = []
    sig_site_real = []
    sig_nuc_idx   = []

    print('%s\nvariant\ts\t(sigma)' % (tag))
    for df_iter, df_entry in df_info.iterrows():
        if df_entry.nucleotide!='-':
            s_val   = df_entry.s_MPL
            idx     = (df_entry.polymorphic_index * 4) + NUC.index(df_entry.nucleotide) - 1
            sigma_s = np.sqrt(ds[idx][idx])
            if np.fabs(s_val)>mult*sigma_s:
                #print('%d\t%.1e (%.1e)' % (idx, s_val, sigma_s))
                n_sig += 1
                sig_s.append(s_val)
                sig_site_real.append(df_entry.alignment_index)
                sig_nuc_idx.append(NUC.index(df_entry.nucleotide) - 1)

    print('')
    print('%d/%d (%d%%) significant at %.2f sigma' % (n_sig, len(df_info), 100*n_sig/len(df_info), mult))
    print('')

    # Sequence landmarks

    seq_range = [[ 790, 2292], [2085, 5096], [5041, 5619], [5559, 5850], [5831, 6045], [5970, 6045], [6045, 6310]]
    seq_range = seq_range +   [[6225, 8795], [8379, 8469], [8379, 8653], [8797, 9417], [9086, 9719]]
    seq_label = [       'gag',        'pol',        'vif',        'vpr', 'tat exon 1', 'rev exon 1',        'vpu']
    seq_label = seq_label +  [        'env', 'tat exon 2', 'rev exon 2',        'nef',     "3' LTR"]

    landmark_start = []
    landmark_end   = []
    landmark_label = []

    idx = df_index.iloc[0].HXB2
    for i in range(len(df_index)):
        if pd.notnull(df_index.iloc[i].HXB2):
            idx = df_index.iloc[i].HXB2
        for j in range(len(seq_range)):
            if idx==seq_range[j][0] or (i==0 and idx>seq_range[j][0] and idx<seq_range[j][1]):
                landmark_start.append(i)
                landmark_end.append(i)
                landmark_label.append(seq_label[j])
            if idx==seq_range[j][1] or (i==len(df_index)-1 and idx>seq_range[j][0] and df_index.iloc[0].HXB2<=seq_range[j][0]
                                        and idx<seq_range[j][1]):
                landmark_end[landmark_label.index(seq_label[j])] = i

    # Epitope labels

    seq_range = epitope_range.copy()
    seq_label = epitope_label.copy()

    epitope_start = []
    epitope_end   = []
    epitope_label = []
    epitope_sites = []
    site2epitope  = {}

    idx = df_index.iloc[0].HXB2
    for i in range(len(df_index)):
        if pd.notnull(df_index.iloc[i].HXB2):
            idx = df_index.iloc[i].HXB2
        for j in range(len(seq_range)):
            if idx==seq_range[j][0] or (i==0 and idx>seq_range[j][0] and idx<seq_range[j][1]):
                epitope_start.append(i)
                epitope_end.append(i)
                epitope_label.append(seq_label[j])
                if seq_label[j]=='DG9':
                    epitope_start[-1] -= 9 # account for DEP insertion
            if idx==seq_range[j][1] or (i==len(df_index)-1 and idx>seq_range[j][0] and df_index.iloc[0].HXB2<=seq_range[j][0]
                                        and idx<seq_range[j][1]):
                epitope_end[epitope_label.index(seq_label[j])] = i
                iix = epitope_label.index(seq_label[j])
                for k in range(epitope_start[iix],epitope_end[iix]):
                    epitope_sites.append(k)
                    if k in site2epitope:
                        print('Unexpected! Overlapping epitopes.')
                    else:
                        site2epitope[k] = seq_label[j]
                print(''.join(list(df_index.iloc[epitope_start[iix]:epitope_end[iix]].TF)))

    # Populate links

    inv_cov = []
    idx_1   = []
    idx_2   = []

    eidx = -1
    if cov_label:
        eidx = epitope_label.index(cov_label)
    cov_cutoff = 0.007
    
    for i in range(len(cinv)):
        for j in range(i+1,len(cinv)):
            if np.fabs(cinv[i][j])>cov_cutoff and (i//len(NUC) != j//len(NUC)):
                ii = df_info[df_info.polymorphic_index==i//len(NUC)].iloc[0].alignment_index
                jj = df_info[df_info.polymorphic_index==j//len(NUC)].iloc[0].alignment_index
                if eidx==-1 or ((ii>=epitope_start[eidx] and ii<=epitope_end[eidx]) or (jj>=epitope_start[eidx] and jj<=epitope_end[eidx])):
                    inv_cov.append(cinv[i][j] / cov_cutoff)
                    idx_1.append(ii)
                    idx_2.append(jj)

    # Dot plot for significant selection coefficients

    c_dot  = { True : C_BEN, False : C_DEL }
    s_mult = 20*SMALLSIZEDOT

    start_idx = 0
    end_idx   = len(df_index)
    level_rad = [0.9, 0.85, 0.80, 0.75]

    def_buff = 100
    def site2angle(site, deg=False, buffer=def_buff):
        if deg: return    -360.*(site+(buffer/2))/(end_idx-start_idx+buffer)
        else:   return -2*np.pi*(site+(buffer/2))/(end_idx-start_idx+buffer)

    dot_colors = [c_dot[s>0]          for s in sig_s]
    dot_sizes  = [s_mult * np.fabs(s) for s in sig_s]
    dot_x = [level_rad[sig_nuc_idx[i]] * np.cos(site2angle(sig_site_real[i])) for i in range(len(sig_s))]
    dot_y = [level_rad[sig_nuc_idx[i]] * np.sin(site2angle(sig_site_real[i])) for i in range(len(sig_s))]

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    for i in range(1, len(NUC)):
        ax.text(level_rad[i-1], 0, NUC[i], **txtprops)

    # Lines for polymorphic sites

    line_r = [0.67, 0.71]
    line_x = [[line_r[0] * np.cos(site2angle(i)), line_r[1] * np.cos(site2angle(i))] for i in np.unique(df_info.alignment_index)]
    line_y = [[line_r[0] * np.sin(site2angle(i)), line_r[1] * np.sin(site2angle(i))] for i in np.unique(df_info.alignment_index)]
    line_c = [LCOLOR                                                                 for i in np.unique(df_info.alignment_index)]

    # Arcs for  sequence landmarks

    arc_r            = [0.96, 0.99, 1.02]
    arc_dr           = 0.01
    track            = 0
    landmark_patches = []
    landmark_text_d  = { 'pol'    : -40,
                         'vif'    :  10,
                         "3' LTR" :  40 }

    for i in range(len(landmark_label)):
        if i>0 and landmark_start[i]<landmark_end[i-1]:
            track = (track + 1)%len(arc_r)
        
        arc_props = dict(center=[0,0], r=arc_r[track], width=arc_dr, lw=AXWIDTH/2, ec=BKCOLOR, fc='none',
                         theta1=site2angle(landmark_end[i], deg=True), theta2=site2angle(landmark_start[i], deg=True))
        landmark_patches.append(matplotlib.patches.Wedge(**arc_props))

        # label with line
        if landmark_label[i] in ['tat exon 1', 'tat exon 2', 'rev exon 1', 'rev exon 2']:
            txtprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
            plotprops = dict(lw=AXWIDTH/2, ls='-', clip_on=False)

            label_x = [arc_r[track]*np.cos(site2angle(landmark_start[i]))]
            label_y = [arc_r[track]*np.sin(site2angle(landmark_start[i]))]
            ddr     = 0.11
            if landmark_label[i] in label2ddr:
                ddr = label2ddr[landmark_label[i]]

            ddx1 =  ddr * np.cos(site2angle(landmark_start[i]))
            ddy  =  ddr * np.sin(site2angle(landmark_start[i]))
            ddx2 = ddx1 + np.sign(ddx1)*0.03

            label_x = label_x + [label_x[0] + ddx1, label_x[0] + ddx2]
            label_y = label_y + [label_y[0] +  ddy, label_y[0] +  ddy]
            if label_x[0]<0:
                txtprops['ha'] = 'right'
            else:
                txtprops['ha'] = 'left'

            ax.text(label_x[-1], label_y[-1], landmark_label[i], **txtprops)
            mp.line(ax=ax, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)

        # plot normally
        else:
            delta_site = 50
            if landmark_label[i] in ['vif']:
                delta_site += landmark_text_d[landmark_label[i]] + 25
            elif landmark_label[i]=='vpu' and 'cap256' in tag:
                delta_site = 25
            elif landmark_label[i]=='env' and 'cap256' in tag:
                delta_site = 75
            elif landmark_label[i] in ["3' LTR"] or (landmark_label[i]=='pol' and '700010077' in tag):
                delta_site += landmark_text_d[landmark_label[i]]
            txt_dr   = 0.04
            txt_ang  = site2angle(landmark_start[i]+delta_site)
            txt_rot  = site2angle(landmark_start[i]+delta_site, deg=True)-90
            txt_x    = (arc_r[track] + (arc_dr/2) + txt_dr)*np.cos(txt_ang)
            txt_y    = (arc_r[track] + (arc_dr/2) + txt_dr)*np.sin(txt_ang)
            txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY,
                            size=SIZELABEL, rotation=txt_rot)
            ax.text(txt_x, txt_y, landmark_label[i], **txtprops)

    # Arcs for TGCA selection tracks
    
    for i in range(len(level_rad)):
        arc_props = dict(center=[0,0], r=level_rad[i], width=0, lw=AXWIDTH/2, ec=C_NEU_LT, fc='none', theta1=0, theta2=360)
        landmark_patches.append(matplotlib.patches.Wedge(**arc_props))

    # Arcs for epitopes

    arc_r           = 1.04
    arc_dr          = 0.32
    epitope_patches = []

    for i in range(len(epitope_label)):
        
        arc_props = dict(center=[0,0], r=arc_r, width=arc_dr, lw=AXWIDTH/2, ec=BKCOLOR, fc='#f2f2f2', alpha=0.8,
                         theta1=site2angle(epitope_end[i], deg=True), theta2=site2angle(epitope_start[i], deg=True))
        epitope_patches.append(matplotlib.patches.Wedge(**arc_props))

        # label epitopes
        if True:
            mid       = (site2angle(epitope_end[i], deg=True)+site2angle(epitope_start[i], deg=True))/2.
            label_x   = [arc_r * np.cos(mid * 2 * np.pi / 360.)]
            label_y   = [arc_r * np.sin(mid * 2 * np.pi / 360.)]
            ddr       = 0.06
            if epitope_label[i] in label2ddr:
                ddr = label2ddr[epitope_label[i]]
            ddx1 =  ddr * np.cos(mid * 2 * np.pi / 360.)
            ddy  =  ddr * np.sin(mid * 2 * np.pi / 360.)
            ddx2 = ddx1 + np.sign(ddx1)*0.03
            txtprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
            plotprops = dict(lw=AXWIDTH/2, ls='-', clip_on=False)

            label_x = label_x + [label_x[0] + ddx1, label_x[0] + ddx2]
            label_y = label_y + [label_y[0] +  ddy, label_y[0] +  ddy]
            if label_x[0]<0:
                txtprops['ha'] = 'right'
            else:
                txtprops['ha'] = 'left'

            ax.text(label_x[-1], label_y[-1], epitope_label[i], **txtprops)
            mp.line(ax=ax, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)

    # Arc plot for large values of integrated covariance

    c_pos  = LCOLOR
    c_neg  = LCOLOR
    c_circ = { True : c_neg, False : c_pos }

    arc_rad   = 0.65
    arc_mult  = SIZELINE/2
    arc_alpha = 0.1

    circ_color = [c_circ[ic>0]                                     for ic in inv_cov]
    #circ_color = [c_circ[idx_1[i] in epitope_sites or idx_2[i] in epitope_sites] for i in range(len(inv_cov))]
    circ_rad   = [[arc_rad, arc_rad]                               for ic in inv_cov]
    circ_arc   = [dict(lw=arc_mult * np.fabs(ic), alpha=arc_alpha) for ic in inv_cov]
    circ_x     = [(i+(def_buff/2)) for i in idx_1]
    circ_y     = [(i+(def_buff/2)) for i in idx_2]

    # Make plot

    pprops = { 'colors':   circ_color,
               'xlim':     [-1.05, 1.05],
               'ylim':     [-1.05, 1.05],
               'size':     end_idx-start_idx+def_buff,
               'bezrad':   0.0,
               'rad':      circ_rad,
               'arcprops': circ_arc,
               'noaxes':   True }

    for i in range(len(dot_x)):
        plotprops = dict(lw=0, marker='o', s=dot_sizes[i], zorder=9999)
        mp.scatter(ax=ax, x=[[dot_x[i]]], y=[[dot_y[i]]], colors=[dot_colors[i]], plotprops=plotprops)

    for i in range(len(np.unique(df_info.alignment_index))):
        plotprops = dict(lw=AXWIDTH, ls='-')
        mp.line(ax=ax, x=[line_x[i]], y=[line_y[i]], colors=[line_c[i]], plotprops=plotprops)

    mp.plot(type='circos', ax=ax, x=circ_x, y=circ_y, **pprops)

    for patch in landmark_patches:
        ax.add_artist(patch)

    for patch in epitope_patches:
        ax.add_artist(patch)

    return sig_s, sig_site_real, sig_nuc_idx, epitope_start, epitope_end


#########################
# SUPPLEMENTARY FIGURES #
#########################

def plot_supplementary_figure_example_mpl(**pdata):
    """
    Example evolutionary trajectory for a 50-site system and inferred selection coefficients,
    together with aggregate properties across sampling levels.
    """
    
    # unpack passed data

    n_gen  = pdata['n_gen']
    dg     = pdata['dg']
    N      = pdata['N']
    xfile  = pdata['xfile']
    method = pdata['method']

    n_ben = pdata['n_ben']
    n_neu = pdata['n_neu']
    n_del = pdata['n_del']
    s_ben = pdata['s_ben']
    s_neu = pdata['s_neu']
    s_del = pdata['s_del']

    # load and process data files

    data  = np.loadtxt('%s/data/%s.dat' % (WFS_DIR, xfile))
    times = np.unique(data.T[0])
    x     = []
    for i in range(0, n_gen, dg):
        idx    = data.T[0]==times[i]
        t_data = data[idx].T[2:].T
        t_num  = data[idx].T[1].T
        t_freq = np.einsum('i,ij->j', t_num, t_data) / float(np.sum(t_num))
        x.append(t_freq)
    x = np.array(x).T

    s_true = [s_ben for i in range(n_ben)] + [0 for i in range(n_neu)] + [s_del for i in range(n_del)]
    s_inf  = np.loadtxt('%s/out/%s_%s.dat' % (MPL_DIR, xfile.split('wfsim_')[1], method))
    cov    = np.loadtxt('%s/out/covariance-%s.dat' % (MPL_DIR, xfile.split('wfsim_')[1]))
    ds     = np.linalg.inv(cov) / N

    # PLOT FIGURE

    ## set up figure grid

    w     = DOUBLE_COLUMN
    goldh = w / 1.08
    fig   = plot.figure(figsize=(w, goldh))

    n_rows = [   n_ben,    n_neu,       n_del]
    offset = [       0,    n_ben, n_ben+n_neu]
    tag    = [   'ben',    'neu',       'del']
    colors = [   C_BEN,    C_NEU,       C_DEL]
    fc     = [C_BEN_LT, C_NEU_LT,    C_DEL_LT]

    htot = 0.80 + 0.01
    dh   = (htot - 0.08) / float(n_ben + n_neu + n_del + 2)

    boxl = [0.06, 0.06, 0.06, 0.06]
    boxr = [0.34, 0.34, 0.34, 0.34]
    boxb = [htot - (dh * n_ben), htot - (dh * (n_ben + n_neu + 1)), htot - (dh * (n_ben + n_neu + n_del + 2))]
    boxt = [               htot,         htot - (dh * (n_ben + 1)), htot - (dh * (n_ben + n_neu + 2))]
    gs   = [0 for k in range(len(tag))]

    ## a -- all trajectories together

    dx = 0.03
    dy = 0.06
    tempgs = gridspec.GridSpec(1, 1)
    tempgs.update(left=boxl[0], right=boxr[0], bottom=htot+dy+0.02, top=0.97)
    ax     = plot.subplot(tempgs[0, 0])

    lineprops = { 'lw' : SIZELINE*1.5, 'ls' : '-', 'alpha' : 0.6 }

    pprops = { 'xticks'      : [0, 50, 100, 150, 200, 250, 300, 350, 400],
               'yticks'      : [0, 1],
               'yminorticks' : [0.25, 0.5, 0.75],
               'nudgey'      : 1.1,
               'xlabel'      : 'Generation',
               'ylabel'      : 'Allele\nfrequency, '+r'$x$',
               'plotprops'   : lineprops,
               'axoffset'    : 0.1,
               'theme'       : 'open' }

    xdat = [range(0, n_gen, dg) for k in range(len(x))]
    ydat = [k for k in x]
    mp.plot(type='line', ax=ax, x=xdat, y=ydat, colors=[LCOLOR for k in range(len(x))], **pprops)

    ## b -- individual beneficial/neutral/deleterious trajectories and selection coefficients

    idx = 0
    for k in range(len(tag)):
        
        ### trajectories
        
        gs[k] = gridspec.GridSpec(n_rows[k], 2)
        gs[k].update(left=boxl[k], right=boxr[k], bottom=boxb[k], top=boxt[k], wspace=0.05)
        ax = [[plot.subplot(gs[k][i, 0]), plot.subplot(gs[k][i, 1])] for i in range(n_rows[k])]
        
        legendprops = { 'loc' : 4, 'frameon' : False, 'scatterpoints' : 1, 'handletextpad' : 0.1,
                        'prop' : {'size' : SIZELABEL}, 'ncol' : 1 }
        lineprops   = { 'lw' : SIZELINE*1.5, 'linestyle' : '-', 'alpha' : 1.0 }
        fillprops   = { 'lw' : 0, 'alpha' : 0.3, 'interpolate' : True }
        
        pprops = { 'xticks' : [],
                   'yticks' : [],
                   'hide'   : ['top','bottom','left','right'],
                   'theme'  : 'open' }
        
        for i in range(n_rows[k]):
            pprops['colors'] = [colors[k]]
            pprops['xlim']   = [    0,  400]
            pprops['ylim']   = [-0.08, 1.08]
            ydat             = x[offset[k]+i]
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']   = [0, 100, 200, 300, 400]
                pprops['xlabel']   = 'Generation'
                pprops['hide']     = ['top','left','right']
                pprops['axoffset'] = 0.3
            mp.line(             ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=lineprops, **pprops)
            mp.plot(type='fill', ax=ax[i][0], x=[range(0, n_gen, dg)], y=[ydat], plotprops=fillprops, **pprops)
        
        ### selection coefficient estimates
        
        sprops = { 'lw' : 0, 's' : 9., 'marker' : 'o' }
        
        pprops = { 'yticks'  : [],
            'xticks'  : [],
            'hide'    : ['top','bottom','left','right'],
            'theme'   : 'open' }
        
        for i in range(n_rows[k]):
            pprops['xlim'] = [-0.04, 0.04]
            pprops['ylim'] = [  0.5,  1.5]
            ydat           = [1]
            xdat           = [s_inf[offset[k]+i]]
            xerr           = np.sqrt(ds[offset[k]+i][offset[k]+i])
            if (i==n_rows[k]-1) and (k==len(tag)-1):
                pprops['xticks']      = [-0.04,   -0.02,      0,   0.02,   0.04]
                pprops['xticklabels'] = [  ' ', r'$-2$', r'$0$', r'$2$', r'$4$']
                pprops['xlabel']      = 'Inferred selection\ncoefficient, ' + r'$\hat{s}$' + ' (%)'
                pprops['hide']        = ['top','left','right']
                pprops['axoffset']    = 0.3
            mp.plot(type='error', ax=ax[i][1], x=[xdat], y=[ydat], xerr=[xerr], colors=[colors[k]], **pprops)
            ax[i][1].axvline(x=s_true[idx], ls=':', lw=SIZELINE, color=BKCOLOR)
            idx += 1

    ### bounding boxes

    ax[0][0].text(boxl[0]-0.03,    0.98, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax[0][0].text(boxl[0]-0.03, boxt[0], 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    lineprops = { 'lw' : AXWIDTH/2., 'ls' : '-', 'alpha' : 1.0 }
    pprops    = { 'xlim' : [0, 1], 'xticks' : [], 'ylim' : [0, 1], 'yticks' : [],
        'hide' : ['top','bottom','left','right'], 'plotprops' : lineprops }
    txtprops = {'ha' : 'right', 'va' : 'center', 'color' : BKCOLOR, 'family' : FONTFAMILY,
        'size' : SIZELABEL, 'rotation' : 90, 'transform' : fig.transFigure}
    ax[0][0].text(boxl[0]-0.01, (boxb[0]+boxt[0])/2.,  'Beneficial', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[1]+boxt[1])/2.,     'Neutral', **txtprops)
    ax[0][0].text(boxl[0]-0.01, (boxb[2]+boxt[2])/2., 'Deleterious', **txtprops)

    boxprops = {'ec' : BKCOLOR, 'lw' : SIZELINE/2., 'fc' : 'none', 'clip_on' : False, 'zorder' : -100}

    dx  = 0.005
    dxl = 0.000
    dxr = 0.001
    dy  = 0.003
    for k in range(len(tag)):
        ll = boxl[k] + dxl                     # left box left
        lb = boxb[k] - dy                      # left box bottom
        rl = (boxl[k] + boxr[k])/2. + dx + dxl # right box left
        wd = (boxr[k] - boxl[k])/2. - dxr      # box width
        ht = (boxt[k] - boxb[k]) + (2. * dy)   # box height
        
        recL = matplotlib.patches.Rectangle(xy=(ll, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recL = ax[0][0].add_patch(recL)
        recR = matplotlib.patches.Rectangle(xy=(rl, lb), width=wd, height=ht, transform=fig.transFigure, **boxprops)
        recR = ax[0][0].add_patch(recR)

    ## c -- histogram of selection coefficients with perfect sampling

    ### set up grid

    box_l  = dict(left=0.44, right=0.58, bottom=0.05, top=0.96)
    box_ru = dict(left=0.69, right=0.97, bottom=0.63, top=0.96)
    box_rl = dict(left=0.69, right=0.97, bottom=0.14, top=0.47)

    gs_l  = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_l)
    gs_ru = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_ru)
    gs_rl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_rl)
    ax_l  = plot.subplot(gs_l[ 0, 0])
    ax_ru = plot.subplot(gs_ru[0, 0])
    ax_rl = plot.subplot(gs_rl[0, 0])

    ### plot histogram

    df   = pd.read_csv('%s/MPL_example_collected_extended.csv.gz' % (SIM_DIR), memory_map=True)
    df   = df[df.method==method]
    df_s = df[(df.deltat==1) & (df.ns==1000)]

    ben_cols = ['s%d' % i for i in range(n_ben)]
    neu_cols = ['s%d' % i for i in range(n_ben, n_ben+n_neu)]
    del_cols = ['s%d' % i for i in range(n_ben+n_neu, n_ben+n_neu+n_del)]

    colors     = [C_BEN, C_NEU, C_DEL]
    tags       = ['beneficial', 'neutral', 'deleterious']
    cols       = [ben_cols, neu_cols, del_cols]
    s_true_loc = [s_ben, s_neu, s_del]

    dashlineprops = { 'lw' : SIZELINE * 2.0, 'ls' : ':', 'alpha' : 0.5, 'color' : BKCOLOR }
    histprops = dict(histtype='bar', lw=SIZELINE/2, rwidth=0.8, ls='solid', alpha=0.7, edgecolor='none',
                     orientation='horizontal')
    pprops = { 'xlim'        : [-0.04, 0.04],
               'xticks'      : [  -0.04,   -0.03,   -0.02,   -0.01,     0.,   0.01,   0.02,   0.03,   0.04],
               'yticklabels' : [r'$-4$', r'$-3$', r'$-2$', r'$-1$', r'$0$', r'$1$', r'$2$', r'$3$', r'$4$'],
               'ylim'        : [0., 0.10],
               'yticks'      : [0., 0.05, 0.10],
               'xlabel'      : 'Frequency',
               'ylabel'      : 'Inferred selection coefficient, ' + r'$\hat{s}$' + ' (%)',
               'bins'        : np.arange(-0.04, 0.04+0.001, 0.001),
               'combine'     : True,
               'plotprops'   : histprops,
               'axoffset'    : 0.1,
               'theme'       : 'boxed' }

    for i in range(len(tags)):
        x = [np.array(df_s[cols[i]]).flatten()]
        tprops = dict(ha='left', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=270, clip_on=False)
        ax_l.text(0.102, s_true_loc[i], r'$s_{%s}$' % (tags[i]), color=colors[i], **tprops)
        ax_l.axhline(y=s_true_loc[i], **dashlineprops)
        if i<len(tags)-1: mp.hist(             ax=ax_l, x=x, colors=[colors[i]], **pprops)
        else:             mp.plot(type='hist', ax=ax_l, x=x, colors=[colors[i]], **pprops)

    ## d, e -- AUCs for inferring beneficial/deleterious mutations

    ns_vals = [10, 20, 30, 40, 50, 80, 100]
    dt_vals = [1, 5, 10, 20, 50]

    AUC_matrix_bn = np.zeros((len(dt_vals), len(ns_vals)))
    AUC_matrix_nd = np.zeros((len(dt_vals), len(ns_vals)))

    for i in range(len(dt_vals)):
        for j in range(len(ns_vals)):
            df_AUC = df[(df.deltat==dt_vals[i]) & (df.ns==ns_vals[j])]
            AUC_matrix_bn[i, j] = np.mean(df_AUC.AUROC_ben)
            AUC_matrix_nd[i, j] = np.mean(df_AUC.AUROC_del)

    pprops = { 'xlim'        : [0, len(dt_vals)],
               'xticks'      : np.arange(len(dt_vals))+0.5,
               'xticklabels' : [int(k) for k in dt_vals],
               'ylim'        : [0, len(ns_vals)],
               'yticks'      : np.arange(len(ns_vals))+0.5,
               'yticklabels' : [int(k) for k in ns_vals],
               'xlabel'      : 'Time between samples, '+r'$\Delta t$' + ' (generations)',
               'ylabel'      : 'Number of samples, '+r'$n_s$',
               'theme'       : 'boxed' }
    tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, clip_on=False)

    ax_ru.pcolor(AUC_matrix_bn.T, vmin=0.75, vmax=1.0, cmap='GnBu', alpha=0.75)
    for i in range(len(AUC_matrix_bn)):
        for j in range(len(AUC_matrix_bn[0])):
            tc = 'k'
            if AUC_matrix_bn[i,j]>0.96: tc = 'white'
            ax_ru.text(i+0.5, j+0.5, '%.2f' % (AUC_matrix_bn[i,j]), color=tc, **tprops)
    mp.plot(type='scatter', ax=ax_ru, x=[[-1]], y=[[-1]], colors=[BKCOLOR], **pprops)

    ax_rl.pcolor(AUC_matrix_nd.T, vmin=0.75, vmax=1.0, cmap='GnBu', alpha=0.75)
    for i in range(len(AUC_matrix_nd)):
        for j in range(len(AUC_matrix_nd[0])):
            tc = 'k'
            if AUC_matrix_nd[i,j]>0.96: tc = 'white'
            ax_rl.text(i+0.5, j+0.5, '%.2f' % (AUC_matrix_nd[i,j]), color=tc, **tprops)
    mp.plot(type='scatter', ax=ax_rl, x=[[-1]], y=[[-1]], colors=[BKCOLOR], **pprops)

    ## outside text labels

    tprops = dict(color=BKCOLOR, ha='center', va='center', family=FONTFAMILY, size=SIZELABEL,
                  clip_on=False, transform=fig.transFigure)
    dx = -0.03
    dy =  0.02

    ax_ru.text((box_ru['right']-box_ru['left'])/2+box_ru['left'],  box_l['top']+dy, 'Mean AUROC (beneficial)',  **tprops)
    ax_ru.text((box_rl['right']-box_rl['left'])/2+box_rl['left'], box_rl['top']+dy, 'Mean AUROC (deleterious)', **tprops)

    ax_l.text(  box_l['left']+dx,  box_l['top']+dy, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_ru.text(box_ru['left']+dx, box_ru['top']+dy, 'd'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_rl.text(box_rl['left']+dx, box_rl['top']+dy, 'e'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # SAVE FIGURE

    plot.savefig('figures/figs1-example-mpl.pdf', dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figure S1 done.')


def plot_supplementary_figure_performance(**pdata):
    """
    Comparisons versus existing methods of selection inference.
    """
    
    # unpack data

    test_sets = pdata['test_sets']
    traj_file = pdata['traj_file']
    t_ticks   = pdata['t_ticks']
    n_ben     = pdata['n_ben']
    n_neu     = pdata['n_neu']
    n_del     = pdata['n_del']
    x_ben     = pdata['x_ben']
    y_ben     = pdata['y_ben']
    x_del     = pdata['x_del']
    y_del     = pdata['y_del']
    x_err     = pdata['x_err']
    y_err     = pdata['y_err']
    x_t       = pdata['x_t']
    y_t       = pdata['y_t']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.90
    goldh   = 0.90 * w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top   = 0.93
    box_left  = 0.25
    box_right = 0.75
    ddy      = 0.06 / hshrink
    dy       = 0.10 / hshrink

    box_traj = dict(left=box_left, right=box_right, bottom=box_top-(1*dy)-(0*ddy), top=box_top-(0*dy)-(0*ddy))
    box_cben = dict(left=box_left, right=box_right, bottom=box_top-(2*dy)-(1*ddy), top=box_top-(1*dy)-(1*ddy))
    box_cdel = dict(left=box_left, right=box_right, bottom=box_top-(3*dy)-(2*ddy), top=box_top-(2*dy)-(2*ddy))
    box_rmse = dict(left=box_left, right=box_right, bottom=box_top-(4*dy)-(3*ddy), top=box_top-(3*dy)-(3*ddy))
    box_time = dict(left=box_left, right=box_right, bottom=box_top-(5*dy)-(4*ddy), top=box_top-(4*dy)-(4*ddy))

    gs_traj = gridspec.GridSpec(1, len(test_sets), wspace=0.05, **box_traj)
    gs_cben = gridspec.GridSpec(1, len(test_sets), wspace=0.05, **box_cben)
    gs_cdel = gridspec.GridSpec(1, len(test_sets), wspace=0.05, **box_cdel)
    gs_rmse = gridspec.GridSpec(1, len(test_sets), wspace=0.05, **box_rmse)
    gs_time = gridspec.GridSpec(1, len(test_sets), wspace=0.05, **box_time)

    ax_traj = [plot.subplot(gs_traj[0, i]) for i in range(len(test_sets))]
    ax_cben = [plot.subplot(gs_cben[0, i]) for i in range(len(test_sets))]
    ax_cdel = [plot.subplot(gs_cdel[0, i]) for i in range(len(test_sets))]
    ax_rmse = [plot.subplot(gs_rmse[0, i]) for i in range(len(test_sets))]
    ax_time = [plot.subplot(gs_time[0, i]) for i in range(len(test_sets))]

    ## set colors and methods list

    hc        = '#FFB511'
    nc        = C_NEU_LT
    methods   = ['MPL', 'FIT', 'LLS', 'CLEAR', 'EandR', 'ApproxWF', 'WFABC',  'IM']
    labels    = ['MPL',   '1',   '2',     '3',     '4',        '5',     '6',   '7']
    colorlist = [   hc,    nc,    nc,      nc,      nc,         nc,      nc,    nc]

    hist_props = dict(lw=SIZELINE/2, width=0.75, align='center', orientation='vertical',
                      edgecolor=[BKCOLOR for i in range(len(methods))])

    ## a -- example trajectory

    for k in range(len(test_sets)):
        data  = np.loadtxt(traj_file[k])
        times = np.unique(data.T[0])
        x     = []
        for i in range(len(times)):
            idx    = data.T[0]==times[i]
            t_data = data[idx].T[2:].T
            t_num  = data[idx].T[1].T
            t_freq = np.einsum('i,ij->j', t_num, t_data) / float(np.sum(t_num))
            x.append(t_freq)
        x = np.array(x).T

        pprops = { 'xlim':      [np.min(times), np.max(times)],
                   'ylim':      [0, 1.05],
                   'xticks':    t_ticks[k],
                   'yticks':    [],
                   'xlabel':    'Generation',
                   'plotprops': {'lw': SIZELINE, 'ls': '-', 'alpha': 0.6 },
                   'axoffset':  0.1,
                   'theme':     'open',
                   'hide':      ['left','right'] }

        if k==0:
            pprops['yticks']      = [0, 1]
            pprops['yminorticks'] = [0.25, 0.5, 0.75]
            pprops['ylabel']      = 'Allele\nfrequency, ' + r'$x$'
            pprops['hide']        = []

        xdat = [times for i in range(n_ben[k])]
        ydat = [i for i in x[:n_ben[k]]]
        pprops['plotprops']['alpha'] = 1
        mp.line(ax=ax_traj[k], x=xdat, y=ydat, colors=[C_BEN_LT for i in range(len(x))], **pprops)

        xdat = [times for i in range(n_neu[k])]
        ydat = [i for i in x[n_ben[k]:n_ben[k]+n_neu[k]]]
        pprops['plotprops']['alpha'] = 0.4
        mp.line(ax=ax_traj[k], x=xdat, y=ydat, colors=[C_NEU for i in range(len(x))], **pprops)

        xdat = [times for i in range(n_del[k])]
        ydat = [i for i in x[n_ben[k]+n_neu[k]:]]
        pprops['plotprops']['alpha'] = 1
        mp.plot(type='line', ax=ax_traj[k], x=xdat, y=ydat, colors=[C_DEL_LT for i in range(len(x))], **pprops)

    ## b -- classification of beneficial mutants

    for k in range(len(test_sets)):
        xmin, xmax = -0.6, len(x_ben[k])-0.6
        ymin, ymax =  0.5, 1.0

        pprops = { 'colors':      [colorlist],
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      x_ben[k],
                   'xticklabels': labels,
                   'yticks':      [],
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks'] = [0.5, 0.75, 1.0]
            pprops['ylabel'] = 'Classification of\nbeneficial alleles\n(AUROC)'
            pprops['hide']   = []
        
        mp.plot(type='bar', ax=ax_cben[k], x=[x_ben[k]], y=[y_ben[k]], plotprops=hist_props, **pprops)

    ## c -- classification of deleterious mutants

    for k in range(len(test_sets)):
        xmin, xmax = -0.6, len(x_err[k])-0.6
        ymin, ymax =  0.5, 1.0

        pprops = { 'colors':      [colorlist],
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      x_del[k],
                   'xticklabels': labels,
                   'yticks':      [],
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks'] = [0.5, 0.75, 1.0]
            pprops['ylabel'] = 'Classification of\ndeleterious alleles\n(AUROC)'
            pprops['hide']   = []
        
        mp.plot(type='bar', ax=ax_cdel[k], x=[x_del[k]], y=[y_del[k]], plotprops=hist_props, **pprops)

    ## d - NRMSE

    for k in range(len(test_sets)):
        xmin, xmax = -0.6, len(x_err[k])-0.6
        ymin, ymax =  0., 3.
        
        print('%s\tNRMSE between %.2f and %.2f' % (test_sets[k], np.min(y_err[k]), np.max(y_err[k])))

        pprops = { 'colors':      [colorlist],
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      x_err[k],
                   'xticklabels': labels,
                   'yticks':      [],
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks']      = [0, 1, 2, 3]
            pprops['ylabel']      = 'Error on inferred\nselection coefficients\n(NRMSE)'
            pprops['hide']        = []
        
        mp.plot(type='bar', ax=ax_rmse[k], x=[x_err[k]], y=[y_err[k]], plotprops=hist_props, **pprops)

    ## e - run time

    for k in range(len(test_sets)):
        xmin, xmax = -0.6, len(x_t[k])-0.6
        ymin, ymax = -3, 6
        
        print('%s\ttime between %.2f and %.2f' % (test_sets[k], np.min(y_t[k]), np.max(y_t[k])))

        pprops = { 'colors':      [colorlist],
                   'xlim':        [xmin, xmax],
                   'ylim':        [ymin, ymax],
                   'xticks':      x_t[k],
                   'xticklabels': labels,
                   'yticks':      [],
                   'theme':       'open',
                   'hide':        ['left','right'] }

        if k==0:
            pprops['yticks']      = [-3, 0, 3, 6]
            pprops['ylabel']      = 'Run time\n' + '(log' +r'$_{10}$' '(seconds))'
            pprops['hide']        = []
        
        mp.plot(type='bar', ax=ax_time[k], x=[x_t[k]], y=[y_t[k]], plotprops=hist_props, **pprops)

        ax_time[k].text(box_left + (1.04*k+0.5)*(box_right - box_left)/len(test_sets), box_traj['top']+0.02,
                        test_sets[k].split('_')[1].capitalize(),
                        ha='center', va='center', transform=fig.transFigure, clip_on=False, **DEF_LABELPROPS)

    # labels and legend

    labelx = 0.18
    ax_traj[0].text(labelx, box_traj['top'] + 0.01, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_cben[0].text(labelx, box_cben['top'] + 0.02, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_cdel[0].text(labelx, box_cdel['top'] + 0.02, 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_rmse[0].text(labelx, box_rmse['top'] + 0.02, 'd'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)
    ax_time[0].text(labelx, box_time['top'] + 0.02, 'e'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    invt = ax_traj[0].transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = xy1[0]-xy2[0]
    legend_dx2 = xy1[0]-xy3[0]
    legend_dy  = xy1[1]-xy2[1]

    legend_x  =  700
    legend_d  = -10
    legend_y  =  1.0
    legend_t  = ['Beneficial', 'Neutral', 'Deleterious']
    legend_c  = [   C_BEN,    C_NEU,    C_DEL]
    legend_cl = [C_BEN_LT, C_NEU_LT, C_DEL_LT]
    plotprops = DEF_ERRORPROPS.copy()
    plotprops['clip_on'] = False
    for k in range(len(legend_t)):
        mp.error(ax=ax_traj[0], x=[[legend_x + legend_d]], y=[[legend_y + (k * legend_dy)]],
                 edgecolor=[legend_c[k]], facecolor=[legend_cl[k]], plotprops=plotprops, **pprops)
        ax_traj[0].text(legend_x, legend_y + (k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)

    for k in range(1, len(methods)):
        ax_traj[0].text(legend_x + legend_d, legend_y + ((k + 3) * legend_dy), labels[k], ha='center', va='center', **DEF_LABELPROPS)
        ax_traj[0].text(legend_x, legend_y + ((k + 3) * legend_dy), methods[k], ha='left', va='center', **DEF_LABELPROPS)

    # Save figure

    plot.savefig('figures/figs2-performance.pdf', **FIGPROPS)
    plot.close(fig)

    print('Figure S2 done.')


def plot_supplementary_figure_absolute_delta_s(**pdata):
    """
    Histogram of absolute values of \Delta s for masked variants in each patient/genomic region.
    """
    
    # unpack data
    
    patient_list = pdata['patient_list']
    region_list  = pdata['region_list']
    ds_values    = pdata['ds_values']
    fig_title    = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.94
    goldh   = w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top    = 0.98
    box_bottom = 0.11
    box_left   = 0.10
    box_right  = 0.90
    ddy      = 0.06 / hshrink
    dy       = 0.10 / hshrink

    box_hist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_hist  = gridspec.GridSpec(6, 5, wspace=0.15, hspace=0.20, **box_hist)
    i_idx    = 0
    j_idx    = 0

    # iterate through patients/regions and plot |\Delta s| histograms

    for k in range(len(patient_list)):
    
        ds_max = 1.281
        #ds_max = 0.006
        n_bins = 50
        bin_ds = ds_max/n_bins
        bins   = np.arange(0, ds_max+bin_ds, bin_ds)
        hist_y = [np.sum((bins[i]<=ds_values[k]) & (ds_values[k]<bins[i+1])) for i in range(len(bins)-1)] #/len(ds_values[k])
        hist_y.append(np.sum(ds_values[k]>=bins[-1]))

        for i in range(len(hist_y)):
            if hist_y[i]>0:
                hist_y[i] = 1 + np.log10(hist_y[i])
            else:
                hist_y[i] = 0

        hist_props = dict(lw=AXWIDTH/2, width=0.9*ds_max/n_bins, align='center', orientation='vertical', edgecolor=[BKCOLOR])

        pprops = { 'xlim':        [0, 1.3],
                   'xticks':      [ 0, 0.4, 0.8, 1.2, 1.3],
                   'xticklabels': ['',  '',  '',  '',  ''],
                   'ylim':        [0, 4],
                   'yticks':      [],
                   'colors':      ['#FFB511'],
                   'plotprops':   hist_props,
                   'axoffset':    0.1,
                   'theme':       'open',
                   'hide':        ['left'] }

        if j_idx==0:
            pprops['yticks']      = [0, 1, 2, 3, 4]
            pprops['yticklabels'] = [0, 1, 10, 100, 1000]
            pprops['ylabel']      = 'Counts'
            pprops['hide'].remove('left')
        
        if i_idx==5 or (j_idx>1 and i_idx==4):
            pprops['xticklabels'] = [0, 40, 80, 120, '']
            pprops['xlabel']      = 'Sum of absolute values\nof effects on inferred\nselection coefficients,\n' + r'$\sum_j\|\Delta \hat{s}_{ij}\|$' + ' (%)'

        ax = plot.subplot(gs_hist[i_idx, j_idx])
        mp.plot(type='bar', ax=ax, x=[bins+(bin_ds/2)], y=[hist_y], **pprops)
        
        tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
        ax.text(0.65, 3.25, patient_list[k].upper() + ' ' + (r'$%d\prime$' % int(region_list[k])), **tprops)
        
        if j_idx<4:
            j_idx += 1
        else:
            i_idx += 1
            j_idx  = 0

    # SAVE FIGURE

    plot.savefig('figures/%s.pdf' % (fig_title), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Delta s histogram done.')


def plot_supplementary_figure_delta_s_distance(**pdata):
    """
    Distribution of \Delta s values by genomic distance between the masked and target variants.
    """
    
    # unpack data
    
    ds_values   = np.array(pdata['ds_values'])
    ds_distance = np.array(pdata['ds_distance'])
    fig_title   = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top    = 0.95
    box_bottom = 0.12
    box_left   = 0.25
    box_right  = 0.75

    box_dist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_dist  = gridspec.GridSpec(1, 1, **box_dist)
    ax_dist  = plot.subplot(gs_dist[0, 0])

    # get \Delta s distribution by distance group
    
    edges = [0, 1, 10, 30, 100, 10000]
    group_ds = []
    
    for i in range(len(edges)-1):
        idxs = (ds_distance>=edges[i]) & (ds_distance<edges[i+1])
        group_ds.append(ds_values[idxs])
    
    ds_max = 0.01
    n_bins = 30
    bin_ds = ds_max / n_bins
    bins = np.arange(0, ds_max+bin_ds, bin_ds)
    bin_x = [bins for i in range(len(group_ds))]
    bin_y = []
    y_min = -7

    for i in range(len(group_ds)):
        temp_y = []
        for j in range(len(bins)-1):
            temp_y.append(np.sum((group_ds[i]>=bins[j]) & (group_ds[i]<bins[j+1])))
        temp_y.append(np.sum(group_ds[i]>=bins[-1]))
        temp_y = np.array(temp_y) / np.sum(temp_y)
        for j in range(len(temp_y)):
            if temp_y[j]>0:
                temp_y[j] = np.log10(temp_y[j])
            else:
                temp_y[j] = y_min
        bin_y.append(temp_y)

    # MAKE STEP STYLE HISTOGRAMS FOR EACH DISTANCE GROUP

    lineprops = { 'lw': SIZELINE*2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }
    #fillprops = { 'lw': 0, 'alpha': 0.2, 'interpolate': True, 'step': 'mid' }

    pprops = { 'xlim':        [-bin_ds, np.max(bins)+bin_ds],
               'xticks':      [0, 0.0025, 0.005, 0.0075, 0.01],
               #'xticklabels': [r'$0$', r'$2.5\times10^{-3}$', r'$5\times10^{-3}$', r'$7.5\times10^{-3}$', r'$\geq 10^{-2}$'],
               'xticklabels': [0, 0.25, 0.5, 0.75, r'$\geq 1$'],
               'ylim':        [-5, 0],
               'yticks':      [-5, -4, -3, -2, -1, 0],
               'yticklabels': ['$10^{-5}$', '$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'],
               'xlabel':      'Absolute value of effect on inferred selection coefficient, ' + r'$\|\Delta \hat{s}_{ij}\|$' + ' (%)',
               'ylabel':      'Frequency',
               'bins':        bins,
               'combine':     True,
               'plotprops':   lineprops,
               'axoffset':    0.1,
               'theme':       'open' }

    h0 =  0.11
    dh =  0 - 0.11
    l0 =  0.53
    dl =  0.73 - 0.53
    s0 =  1.00
    ds = -1.00

    for k in range(len(group_ds)-1):
        x = k/(len(group_ds) - 1)
        c = hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x))
        # Yellow 0.11 0.53 1.00
        # Gray   0.00 0.59 0.00
        mp.line(ax=ax_dist, x=[bin_x[k]], y=[bin_y[k]], colors=[c], zorder=-k, **pprops)

    c = hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)
    mp.plot(type='line', ax=ax_dist, x=[bin_x[-1]], y=[bin_y[-1]], colors=[c], zorder=-len(group_ds), **pprops)

    # Plot legend

    invt = ax_dist.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1*(xy1[1]-xy2[1])

    legend_x =  0.012
    legend_y = -0.5
    legend_d = -0.0005
    legend_t = ['$0$', '$1-10$', '$11-30$', '$31-100$', '$>100$']
    legend_c = [hls_to_rgb(h0 + (dh * x), l0 + (dl * x), s0 + (ds * x)) for x in np.arange(0, 1+1/len(group_ds), 1/len(group_ds))]
    for k in range(len(legend_t)):
        mp.line(ax=ax_dist, x=[[legend_x + legend_dx1, legend_x + legend_dx2]],
                y=[[legend_y + (k * legend_dy), legend_y + (k * legend_dy)]],
                colors=[legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_dist.text(legend_x, legend_y + (k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)
    ax_dist.text(legend_x+0.0005, legend_y - (1.75*legend_dy), 'Distance between variant i\nand target variant j (bp)', ha='center', va='center', **DEF_LABELPROPS)

    # SAVE FIGURE

    plot.savefig('figures/%s.pdf' % (fig_title), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Delta s distance distribution done.')


def plot_supplementary_figure_delta_s_hive(**pdata):
    """
    Hive plot of \Delta s effects.
    """

    # unpack data
    
    patient_list = pdata['patient_list']
    region_list  = pdata['region_list']
    fig_title    = pdata['fig_title']

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 1.2
    goldh   = w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top    = 1.00
    box_bottom = 0.05
    box_left   = 0.025
    box_right  = 0.975

    box_hive = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_hive  = gridspec.GridSpec(6, 5, wspace=0, hspace=0, **box_hive)
    i_idx    = 0
    j_idx    = 0

    # iterate through patients/regions and plot hive plots

    for k in range(len(patient_list)):
    
        ax = plot.subplot(gs_hive[i_idx, j_idx])
        plot_hive(ax, patient_list[k]+'-'+region_list[k])
        
        tprops = dict(ha='center', va='center', family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
        ax.text(0, -1.15, patient_list[k].upper() + ' ' + (r'$%d\prime$' % int(region_list[k])), **tprops)

        if i_idx<4 and j_idx<4:
            j_idx += 1
        elif i_idx>=4 and j_idx<2:
            j_idx += 1
        else:
            i_idx += 1
            j_idx  = 0

    # plot legend

    ax_hive = plot.subplot(gs_hive[-2:, -2:])
    
    idx_mask_pos  = [90]
    idx_delta_pos = [90]
    ds_pos        = [ 1]
    idx_mask_neg  = [70]
    idx_delta_neg = [70]
    ds_neg        = [ 1]

    # arc plot for large values of Delta s

    L         = 100
    ds_cutoff = 0.004
    ds_max    = 0.01
    r_min     = 0.025
    r_max     = 0.50
    r_norm    = float(L) / (r_max - r_min)
    arc_mult  = SIZELINE * 3 * ds_max
    arc_alpha = 1

    mask_angle_p =  -np.pi/2
    mask_angle_n = 3*np.pi/2
    ds_pos_angle =   np.pi/6
    ds_neg_angle = 5*np.pi/6
    da           = 0.005

    circ_color  = [hls_to_rgb(0.02, 0.53 * ds + 1. * (1 - ds), 0.83) for ds in ds_pos]
    circ_color += [hls_to_rgb(0.58, 0.53 * ds + 1. * (1 - ds), 0.60) for ds in ds_neg]
    circ_rad    = [[r_min + idx_mask_pos[i]/r_norm, r_min + idx_delta_pos[i]/r_norm] for i in range(len(ds_pos))]
    circ_rad   += [[r_min + idx_mask_neg[i]/r_norm, r_min + idx_delta_neg[i]/r_norm] for i in range(len(ds_neg))]
    circ_arc    = [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_pos]
    circ_arc   += [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_neg]
    circ_angle  = [[mask_angle_p+da/r[0], ds_pos_angle-da/r[1]] for r in circ_rad[:len(ds_pos)]]
    circ_angle += [[mask_angle_n-da/r[0], ds_neg_angle+da/r[1]] for r in circ_rad[len(ds_pos):]]
    circ_bezrad = [(r[0]+r[1])/2 for r in circ_rad]
    circ_x      = [i for i in ds_pos] + [i for i in ds_neg]
    circ_y      = [i for i in ds_pos] + [i for i in ds_neg]
    
    # CD8+ T cell epitope legend

    r_mid = r_min + (20 / r_norm)
    scatter_x = [r_mid * np.cos(mask_angle_p),
                 r_mid * np.cos(ds_pos_angle),
                 r_mid * np.cos(ds_neg_angle) ]
    scatter_y = [r_mid * np.sin(mask_angle_p),
                 r_mid * np.sin(ds_pos_angle),
                 r_mid * np.sin(ds_neg_angle) ]
    smallprops = dict(lw=0, marker='o', s=1.0*SMALLSIZEDOT, zorder=9999, clip_on=False)
    bigprops   = dict(lw=0, marker='o', s=1.7*SMALLSIZEDOT, zorder=9998, clip_on=False)
    mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=['#FFFFFF'], plotprops=smallprops)
    mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=[BKCOLOR],   plotprops=  bigprops)
    
    # add legend labels

    ## epitope label
    label_r  = r_mid
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(ds_neg_angle)]
    label_y  = [label_r * np.sin(ds_neg_angle) + 0.5*ddy]
    label_x  = label_x + [label_x[0]]
    label_y  = label_y + [label_y[0] + 2*ddy]

    txtprops  = dict(ha='center', va='bottom', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1], label_y[-1] + 0.125*ddy, 'CD8+ T cell\nepitope', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## 5' -- 3' label
    arrow_x  = np.array([r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)]) - 0.03/2
    arrow_y  = np.array([r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)]) + 0.10/2
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=180*ds_pos_angle/np.pi)
    
    ddr = 0.05/2
    ax_hive.text(arrow_x[0] + 1.0 * ddr * np.cos(ds_pos_angle), arrow_y[0] + 0.1/2 + 1.0 * ddr * np.sin(ds_pos_angle), r'$5\prime$', **txtprops)
    ax_hive.text(arrow_x[1] - 2.5 * ddr * np.cos(ds_pos_angle), arrow_y[1] + 0.1/2 - 2.5 * ddr * np.sin(ds_pos_angle), r'$3\prime$', **txtprops)

    ax_hive.annotate('', xy=(arrow_x[1], arrow_y[1]), xytext=(arrow_x[0], arrow_y[0]),
                     arrowprops=dict(color=BKCOLOR, lw=0, shrink=0.0, width=AXWIDTH, headwidth=AXWIDTH*5, headlength=AXWIDTH*5),)
    
    ## focus/target labels
    label_r = 1.15/2
    tprops  = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(label_r * np.cos(mask_angle_n), label_r * np.sin(mask_angle_n), 'Variant i', **tprops)
    label_r = 1.1/2
    tprops['ha'] = 'left'
    ax_hive.text(label_r * np.cos(ds_pos_angle), label_r * np.sin(ds_pos_angle), 'Target\nvariant j', **tprops)
    tprops['ha'] = 'right'
    ax_hive.text(label_r * np.cos(ds_neg_angle), label_r * np.sin(ds_neg_angle), 'Target\nvariant j', **tprops)
    
    ## negative link label
    label_r  = r_min + idx_mask_neg[0]/r_norm - 0.11/2
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(np.pi)]
    label_y  = [label_r * np.sin(np.pi) - 0.5*ddy]
    label_x  = label_x + [label_x[0] - ddx, label_x[0] - 1.5*ddx]
    label_y  = label_y + [label_y[0] - ddy, label_y[0] -     ddy]
    
    txtprops  = dict(ha='right', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1] - 0.5*ddx, label_y[-1], 'Variant i decreases\ntarget ' + r'$\hat{s}_j$' + ', ' + r'$\Delta \hat{s}_{ij}<0$', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## positive link label
    label_r  = r_min + idx_mask_pos[0]/r_norm - 0.15/2
    ddx, ddy = 0.09/2, 0.15/2
    label_x  = [label_r * np.cos(0)]
    label_y  = [label_r * np.sin(0) - 0.5*ddy]
    label_x  = label_x + [label_x[0] + ddx, label_x[0] + 1.5*ddx]
    label_y  = label_y + [label_y[0] - ddy, label_y[0] -     ddy]
    
    txtprops  = dict(ha='left', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0)
    plotprops = dict(lw=AXWIDTH, ls='-', clip_on=False)
    
    ax_hive.text(label_x[-1] + 0.5*ddx, label_y[-1], 'Variant i increases\ntarget ' + r'$\hat{s}_j$' + ', ' + r'$\Delta \hat{s}_{ij}>0$', **txtprops)
    mp.line(ax=ax_hive, x=[label_x], y=[label_y], colors=[BKCOLOR], plotprops=plotprops)
    
    ## color bar
    ddx      =  0.425
    ddy      = -0.10
    label_x  =  0
    label_y  = -0.96
    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(label_x+ddx, label_y+ddy, 'Effect on inferred selection\ncoefficient, ' + r'$\Delta \hat{s}_{ij}$' + ' (%)', **txtprops)
    
    dbox      = 0.05
    rec_props = dict(height=dbox, width=dbox, ec=None, lw=AXWIDTH/2, clip_on=False)
    for i in range(-5, 5+1, 1):
        c = BKCOLOR
        t = i/5
        if t>0:
            c = hls_to_rgb(0.02, 0.53 * t + 1. * (1 - t), 0.83)
        else:
            c = hls_to_rgb(0.58, 0.53 * np.fabs(t) + 1. * (1 - np.fabs(t)), 0.60)
        rec = matplotlib.patches.Rectangle(xy=((i-0.5)*dbox + ddx, -0.75 + ddy), fc=c, **rec_props)
        ax_hive.add_artist(rec)

    txtprops = dict(ha='center', va='top', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL, rotation=0, clip_on=False)
    ax_hive.text(-5.25*dbox + ddx, -0.78 + ddy, -1, **txtprops)
    ax_hive.text(             ddx, -0.78 + ddy,  0, **txtprops)
    ax_hive.text( 5.00*dbox + ddx, -0.78 + ddy,  1, **txtprops)
    
    # Make plot

    pprops = { 'colors':   circ_color,
               'xlim':     [-1.0, 1.0],
               'ylim':     [-1.0, 1.0],
               'size':     L,
               'rad':      circ_rad,
               'arcprops': circ_arc,
               'angle':    circ_angle,
               'bezrad':   circ_bezrad,
               'noaxes':   True }

    plotprops = dict(lw=AXWIDTH, ls='-', zorder=999)
    line_x = [[r_min * np.cos(mask_angle_p), r_max * np.cos(mask_angle_p)],
              [r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)],
              [r_min * np.cos(ds_neg_angle), r_max * np.cos(ds_neg_angle)] ]
    line_y = [[r_min * np.sin(mask_angle_p), r_max * np.sin(mask_angle_p)],
              [r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)],
              [r_min * np.sin(ds_neg_angle), r_max * np.sin(ds_neg_angle)] ]
    mp.line(ax=ax_hive, x=line_x, y=line_y, colors=[BKCOLOR for i in range(len(line_x))], plotprops=plotprops)

    mp.plot(type='circos', ax=ax_hive, x=circ_x, y=circ_y, **pprops)

    # SAVE FIGURE

    plot.savefig('figures/%s.pdf' % (fig_title), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Hive plots done.')


def plot_supplementary_figure_max_dx(**pdata):
    """
    Distribution of maximum observed change in variant frequency per generation.
    """
    
    # unpack data
    
    dx_sweep  = np.array(pdata['dx_sweep'])
    dx_other  = np.array(pdata['dx_other'])
    fig_title = pdata['fig_title']

    # PLOT FIGURE
    
    # SHOW FREQUENCY INSTEAD

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.55
    goldh   = w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top    = 0.95
    box_bottom = 0.12
    box_left   = 0.25
    box_right  = 0.75

    box_dist = dict(left=box_left, right=box_right, bottom=box_bottom, top=box_top)
    gs_dist  = gridspec.GridSpec(1, 1, **box_dist)
    ax_dist  = plot.subplot(gs_dist[0, 0])

    # get dx bins
    
    dx_max  = 0.143
    n_bins  = 20
    bin_dx  = dx_max/n_bins
    bins    = np.arange(0, dx_max+bin_dx, bin_dx)
    hist_ys = []
    for dx_values in [dx_sweep, dx_other]:
        dx_values = np.array(dx_values)
        hist_y = [np.sum((bins[i]<=dx_values) & (dx_values<bins[i+1])) for i in range(len(bins)-1)]
        hist_y.append(np.sum(dx_values>=bins[-1]))
        hist_ys.append(np.array(hist_y)/len(dx_values))

    for k in range(len(hist_ys)):
        for i in range(len(hist_ys[k])):
            if hist_ys[k][i]>0:
                hist_ys[k][i] = np.log10(hist_ys[k][i])
            else:
                hist_ys[k][i] = -10

    # make figure

    lineprops = { 'lw': SIZELINE*2, 'linestyle': '-', 'alpha': 1.0, 'drawstyle': 'steps-mid' }
    hist_props = dict(lw=AXWIDTH/2, width=0.9*dx_max/n_bins, align='center',
                      orientation='vertical', edgecolor=[BKCOLOR])

    pprops = { 'xlim':        [-bin_dx, np.max(bins)+bin_dx],
               'xticks':      [ 0, 0.03, 0.06, 0.09, 0.12, 0.15],
               'xticklabels': [ 0,    3,    6,    9,   12,   15],
               'ylim':        [-4,  0],
               'yticks':      [-4, -3, -2, -1, 0],
               'yticklabels': ['$10^{-4}$', '$10^{-3}$', '$10^{-2}$', '$10^{-1}$', '1'],
               'xlabel':      'Maximum observed change in variant frequency per day (%)',
               'ylabel':      'Frequency',
               #'colors':      ['#FFB511'],
               'plotprops':   lineprops,
               'axoffset':    0.1,
               'theme':       'open' }

#    mp.plot(type='bar', ax=ax, x=[bins+(bin_ds/2), bins+(bin_ds/2)], y=hist_ys, **pprops)

    h0 =  0.11
    dh =  0 - 0.11
    l0 =  0.53
    dl =  0.73 - 0.53
    s0 =  1.00
    ds = -1.00

    c = hls_to_rgb(h0, l0, s0)
    mp.line(ax=ax_dist, x=[bins], y=[hist_ys[0]], colors=[c], **pprops)

    c = hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)
    mp.plot(type='line', ax=ax_dist, x=[bins], y=[hist_ys[-1]], colors=[c], **pprops)

    # Plot legend

    invt = ax_dist.transData.inverted()
    xy1  = invt.transform((  0, 0))
    xy2  = invt.transform((7.5, 9))
    xy3  = invt.transform((3.0, 9))

    legend_dx1 = 1*(xy1[0]-xy2[0])
    legend_dx2 = 1*(xy1[0]-xy3[0])
    legend_dy  = 1.2*(xy1[1]-xy2[1])

    legend_x =  0.17
    legend_y = -0.5
    legend_d = -0.0005
    legend_t = ['Highly influential variants', 'Other variants']
    legend_c = [hls_to_rgb(h0, l0, s0), hls_to_rgb(h0 + dh, l0 + dl, s0 + ds)]
    for k in range(len(legend_t)):
        mp.line(ax=ax_dist, x=[[legend_x + legend_dx1, legend_x + legend_dx2]],
                y=[[legend_y + (1.0 * k * legend_dy), legend_y + (1.0 * k * legend_dy)]],
                colors=[legend_c[k]], plotprops=dict(lw=2*SIZELINE, ls='-', clip_on=False))
        ax_dist.text(legend_x, legend_y + (1.0 * k * legend_dy), legend_t[k], ha='left', va='center', **DEF_LABELPROPS)
    #ax_dist.text(legend_x+0.0005, legend_y - (1.75*legend_dy), 'Distance between variant i\nand target variant j (bp)', ha='center', va='center', **DEF_LABELPROPS)

    # SAVE FIGURE

    plot.savefig('figures/%s.pdf' % (fig_title), dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Delta x distribution done.')


def plot_supplementary_figure_epitope(**pdata):
    """
    Alternate epitope escape mutation frequencies, inferred selection coefficients, and linkage.
    """

    # unpack data
    
    patient       = pdata['patient']
    region        = pdata['region']
    epitope       = pdata['epitope']
    epitope_range = pdata['epitope_range']
    epitope_label = pdata['epitope_label']
    cov_label     = pdata['cov_label']
    label2ddr     = pdata['label2ddr']
    legend_loc    = pdata['legend_loc']
    traj_ticks    = pdata['traj_ticks']
    sel_ticks     = pdata['sel_ticks']
    sel_minors    = pdata['sel_minors']
    sel_space     = pdata['sel_space']
    fig_title     = pdata['fig_title']
    tag           = patient+'-'+region
    
    # process stored data
    
    df_poly = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_esc  = df_poly[(df_poly.epitope==epitope) & (df_poly.nucleotide!=df_poly.TF)]

    times = [int(i.split('_')[-1]) for i in df_esc.columns if 'f_at_' in i]
    times.sort()
    
    var_tag  = []
    var_smpl = []
    var_sind = []
    var_traj = []
    curr_HXB2 = 8988
    curr_aln  = 4015
    curr_char = 'a'
    for df_iter, df_entry in df_esc.iterrows():
        if df_entry.nucleotide=='-':
            continue
        if pd.notnull(df_entry.HXB2_index):
            var_tag.append(str(int(df_entry.HXB2_index))+df_entry.nucleotide)
            curr_HXB2 = int(df_entry.HXB2_index)
            curr_aln  = int(df_entry.alignment_index)
        else:
            if int(df_entry.alignment_index)!=curr_aln:
                curr_aln  = int(df_entry.alignment_index)
                curr_char = chr(ord(curr_char) + 1)
            var_tag.append(str(int(curr_HXB2))+curr_char+df_entry.nucleotide)
        var_traj.append([df_entry['f_at_%d' % t] for t in times])
        var_smpl.append(df_entry.s_MPL)
        var_sind.append(df_entry.s_SL)

    var_c = sns.husl_palette(len(var_traj))

    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN
    hshrink = 0.93
    goldh   = 1 * w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_top  = 0.95
    dy       = 0.10
    y0       = 0.09/hshrink
    y1       = 0.08/hshrink
    y2       = 0.50/hshrink
    s_spc    = np.max([5, len(var_c)]) * 0.02
    box_traj = dict(left=0.30,                 right=0.70,                     bottom=box_top-y0,                top=box_top)
    box_smpl = dict(left=0.30,                 right=0.30+s_spc,               bottom=box_top-y0-y1-(0.8*dy),    top=box_top-y0-(0.8*dy))
    box_sind = dict(left=0.30+s_spc+sel_space, right=0.30+(2*s_spc)+sel_space, bottom=box_top-y0-y1-(0.8*dy),    top=box_top-y0-(0.8*dy))
    box_circ = dict(left=0.25,                 right=0.75,                     bottom=box_top-y0-y1-y2-(1.9*dy), top=box_top-y0-y1-(1.9*dy))

    gs_traj = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_traj)
    gs_smpl = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_smpl)
    gs_sind = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_sind)
    gs_circ = gridspec.GridSpec(1, 1, width_ratios=[1.0], height_ratios=[1.0], **box_circ)
    ax_traj = plot.subplot(gs_traj[0, 0])
    ax_smpl = plot.subplot(gs_smpl[0, 0])
    ax_sind = plot.subplot(gs_sind[0, 0])
    ax_circ = plot.subplot(gs_circ[0, 0])

    dy = 0.02/hshrink

    ## a -- trajectory plot

    pprops = { 'xticks':      traj_ticks,
               'yticks':      [0, 1],
               'yminorticks': [0.25, 0.5, 0.75],
               'nudgey':      1.1,
               'xlabel':      'Time (days)',
               'ylabel':      'Variant frequency\nin %s epitope\n' % cov_label,
               'plotprops':   {'lw': SIZELINE*1.5, 'ls': '-', 'alpha': 1.0 },
               'axoffset':    0.1,
               'theme':       'open' }

    for i in range(len(var_tag)-1):
        xdat = [times]
        ydat = [var_traj[i]]
        mp.line(ax=ax_traj, x=xdat, y=ydat, colors=[var_c[i]], **pprops)

    xdat = [times]
    ydat = [var_traj[len(var_tag)-1]]
    mp.plot(type='line', ax=ax_traj, x=xdat, y=ydat, colors=[var_c[len(var_tag)-1]], **pprops)

    ax_traj.text(0.25, box_traj['top']+dy, 'a'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b1 -- selection coefficients inferred by MPL

    xmin, xmax = 0, len(var_tag)
    diff = 5. - len(var_c)
    if diff>0:
        xmin -= diff/2
        xmax += diff/2

    hist_props = dict(lw=SIZELINE/2, width=0.6, align='center', orientation='vertical',
                      edgecolor=[BKCOLOR for i in range(len(var_tag))])

    bar_x  = [i+0.5 for i in range(len(var_tag))]
    pprops = { 'colors':      [var_c],
               'xlim':        [xmin, xmax],
               'ylim':        [np.min(sel_ticks), np.max(sel_ticks)],
               'xticks':      bar_x,
               'xticklabels': var_tag,
               'yticks':      sel_ticks,
               'yticklabels': [int(100*l) for l in sel_ticks],
               'yminorticks': sel_minors,
               'ylabel':      'Inferred selection\ncoefficient, '+r'$\hat{s}$'+' (%)',
               'theme':       'open',
               'hide' :       [] }

    mp.plot(type='bar', ax=ax_smpl, x=[bar_x], y=[var_smpl], plotprops=hist_props, **pprops)
    plot.setp(ax_smpl.xaxis.get_majorticklabels(), rotation=90)

    transFigureInv = fig.transFigure.inverted()
    labelprops     = dict(color=BKCOLOR, ha='center', va='top', family=FONTFAMILY, size=SIZELABEL,
                          clip_on=False, transform=fig.transFigure)
    ax_smpl.text(box_smpl['left']+(box_smpl['right']-box_smpl['left'])/2, box_smpl['top'], 'MPL', **labelprops)

    ax_smpl.text(0.25, box_smpl['top']+dy, 'b'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    ## b2 -- selection coefficients inferred by SL

    pprops = { 'colors':      [var_c],
               'xlim':        [xmin, xmax],
               'ylim':        [np.min(sel_ticks), np.max(sel_ticks)],
               'xticks':      bar_x,
               'xticklabels': var_tag,
               'yticks':      [],
               'theme':       'open',
               'hide' :       ['left','right'] }

    mp.plot(type='bar', ax=ax_sind, x=[bar_x], y=[var_sind], plotprops=hist_props, **pprops)
    plot.setp(ax_sind.xaxis.get_majorticklabels(), rotation=90)

    ax_sind.text(box_sind['left']+(box_sind['right']-box_sind['left'])/2, box_sind['top'], 'Independent\nmodel', **labelprops)

    # add background

    cBG = '#F5F5F5'
    bg  = ax_sind.axis()
    ddx = 0.01
    ddy = 0.01
    rec = matplotlib.patches.Rectangle(xy=(box_sind['left']-(0.1*ddx), box_sind['bottom']-(0.7*ddy)), width=box_sind['right']-box_sind['left']+(0.2*ddx),
                                       height=box_sind['top']-box_sind['bottom']+(1.7*ddy), transform=fig.transFigure, ec=None, fc=cBG, clip_on=False, zorder=-100)
    rec = ax_sind.add_patch(rec)

    ## c -- circle plot

    sig_s, sig_site_real, sig_nuc_idx, epitope_start, epitope_end = plot_circle(ax_circ, tag, epitope_range, epitope_label, cov_label, label2ddr)
    
    ax_circ.text(0.25, box_circ['top'], 'c'.upper(), transform=fig.transFigure, **DEF_SUBLABELPROPS)

    # MAKE LEGEND

    invt = ax_circ.transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((0,9))
    
    coef_legend_x  = 1.30
    coef_legend_dx = 0.05
    coef_legend_y  = 1.00
    coef_legend_dy = xy1[1]-xy2[1]

    txtprops = dict(ha='center', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    s_mult   = 20*SMALLSIZEDOT
    ex_s     = [-0.1, -0.08, -0.06, -0.04, -0.02, 0, 0.02, 0.04, 0.06, 0.08, 0.10]
    show_s   = [   1,     0,     0,     0,     0, 1,    0,    0,    0,    0,    1]
    c_s      = [C_DEL, C_BEN]
    for i in range(len(ex_s)):
        plotprops = dict(lw=0, marker='o', s=np.fabs(ex_s[i])*s_mult, clip_on=False)
        mp.scatter(ax=ax_circ, x=[[coef_legend_x + i*coef_legend_dx]], y=[[coef_legend_y]], colors=[c_s[ex_s[i]>0]], plotprops=plotprops)
        if show_s[i]:
            ax_circ.text(coef_legend_x + i*coef_legend_dx, coef_legend_y + 0.75*coef_legend_dy, '%d' % (100*ex_s[i]), **txtprops)
    ax_circ.text(coef_legend_x + 5*coef_legend_dx, coef_legend_y + (2.25*coef_legend_dy), 'Inferred selection\ncoefficient, $\hat{s}$ (%)', **txtprops)

    coef_legend_x = 1.30
    coef_legend_y = coef_legend_y + 5.5 * coef_legend_dy
    txtprops = dict(ha='left', va='center', color=BKCOLOR, family=FONTFAMILY, size=SIZELABEL)
    arc_props = dict(lw=AXWIDTH/2, ec=BKCOLOR, fc='#f2f2f2', alpha=0.8, clip_on=False)
    epatch    = matplotlib.patches.Rectangle(xy=(coef_legend_x-0.025, coef_legend_y+(0.0*coef_legend_dy)), width=0.05, height=-1.*coef_legend_dy, **arc_props)
    ax_circ.add_artist(epatch)
    ax_circ.text(coef_legend_x+(1*coef_legend_dx), coef_legend_y-(coef_legend_dy/2), 'CD8+ T cell\nepitope', **txtprops)

    # SAVE FIGURE

    plot.savefig('figures/%s.pdf' % fig_title, dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figures S7-S8 done.')


def plot_supplementary_figure_cap256_recombination(**pdata):
    """
    Plot patterns of recombination in CAP256.
    """

    # unpack data
    
    patient = pdata['patient']
    region  = pdata['region']
    tag     = patient+'-'+region

    # process stored data
    
    df_index = pd.read_csv('%s/processed/%s-SU-%s-index.csv' % (HIV_DIR, patient, region), comment='#', memory_map=True)
    df_poly  = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    seqs     = np.loadtxt('%s/HIV/%s-poly-seq2state.dat' % (MPL_DIR, tag), dtype='int')
    times    = np.unique(seqs.T[0])
    
    c_vals   = [C_DEL, C_NEU]
    m_loc    = [[] for s in seqs]
    m_c      = [[] for s in seqs]
    for i in range(2, len(seqs[0])):
        df_temp = df_index[df_index.polymorphic==(i-2)]
        SU  = str(df_temp.iloc[0].SU)
        TF  = str(df_temp.iloc[0].TF)
        loc = int(df_temp.iloc[0].alignment)
        for j in range(len(seqs)):
            n = NUC[seqs[j][i]]
            # SU, and SU does not match TF
            if n==SU and SU!=TF:
                m_loc[j].append(loc)
                m_c[j].append(c_vals[0])
            # Other difference from TF
            elif n!=TF:
                m_loc[j].append(loc)
                m_c[j].append(c_vals[1])

    m_loc = np.array(m_loc)
    m_c   = np.array(m_c)
    
    # PLOT FIGURE

    ## set up figure grid

    w       = DOUBLE_COLUMN #SLIDE_WIDTH
    hshrink = 1
    goldh   = 0.66 * w * hshrink
    fig     = plot.figure(figsize=(w, goldh))

    box_hist = dict(left=0.15, right=0.85, bottom=0.05, top=0.95)
    gs_hist  = gridspec.GridSpec(len(times), 1, hspace=0.15, wspace=0.15, height_ratios=[np.sum(seqs.T[0]==t) for t in times], **box_hist)
    ax_hist  = [plot.subplot(gs_hist[i, 0]) for i in range(len(times))]

    for k in range(len(times)):
        
        # Subplot - recombination
        
        sp         = 0.2
        rect_props = dict(width=1, height=1, lw=0, alpha=1, ec='none')
        patches    = []
        
        same_time = seqs.T[0]==times[k]
        sub_m_loc = m_loc[same_time]
        sub_m_c   = m_c[same_time]
        
        for i in range(len(sub_m_loc)):
            rect_props['width'] = 1
            for j in range(len(sub_m_loc[i])):
                patches.append(matplotlib.patches.Rectangle(xy=(sub_m_loc[i][j], (1+sp)*i), fc=sub_m_c[i][j], **rect_props))

            start = 0
            for j in range(len(sub_m_loc[i])):
                rect_props['width'] = sub_m_loc[i][j]-start
                patches.append(matplotlib.patches.Rectangle(xy=(start, (1+sp)*i), fc='#f7f7f7', **rect_props))
                start = sub_m_loc[i][j]+1
            rect_props['width'] = 2645-start
            patches.append(matplotlib.patches.Rectangle(xy=(start, (1+sp)*i), fc='#f7f7f7', **rect_props))

        pprops = { 'colors': [BKCOLOR],
                   'xlim': [-1, 2645],
                   'ylim': [ 0, (1+sp)*len(sub_m_loc)],
                   'xticks': [],
                   'yticks': [],
                   'hide': ['top', 'bottom', 'left', 'right']}

        mp.plot(type='scatter', ax=ax_hist[k], x=[[-3]], y=[[-3]], **pprops)

        for patch in patches:
            ax_hist[k].add_artist(patch)

        rec = matplotlib.patches.Rectangle(xy=(0, 0), height=(1+sp)*len(sub_m_loc), width=2645, ec=BKCOLOR, fc='none', lw=AXWIDTH/2, clip_on=False)
        ax_hist[k].add_artist(rec)

    # outside text

    invt = ax_hist[0].transData.inverted()
    xy1  = invt.transform((0,0))
    xy2  = invt.transform((9,9))
    legend_dx = 0.5 * (xy1[0] - xy2[0])
    legend_dy = xy1[1] - xy2[1]
    legend_x  = 2645 + 120
    legend_y  = (1 + 0.2) * np.sum(seqs.T[0]==times[0]) + 0.75 * legend_dy

    rec_props = dict(height=0.5*np.fabs(legend_dy), width=0.5*np.fabs(legend_dx), clip_on=False)
    rec = matplotlib.patches.Rectangle(xy=(legend_x, legend_y + (-0.25*legend_dy)), fc=C_DEL, **rec_props)
    ax_hist[0].add_artist(rec)
    rec = matplotlib.patches.Rectangle(xy=(legend_x, legend_y + ( 1.25*legend_dy)), fc=C_NEU, **rec_props)
    ax_hist[0].add_artist(rec)

    tprops = dict(ha='left', va='center', family=FONTFAMILY, size=SIZELABEL, color=BKCOLOR, clip_on=False)

    ax_hist[0].text(legend_x - legend_dx,  legend_y + (-0.5*legend_dy), 'Superinfecting\nvariants', **tprops)
    ax_hist[0].text(legend_x - legend_dx,  legend_y + ( 1.0*legend_dy), 'Other variants', **tprops)

    tprops['ha'] = 'center'
    for i in range(len(times)):
        ax_hist[i].text(-30, 0.5*(1+sp)*np.sum(seqs.T[0]==times[i]), times[i], rotation=90, **tprops)

    ax_hist[-1].text(0.13,  0.50, 'Time (days)', rotation=90, transform=fig.transFigure, **tprops)
    ax_hist[-1].text(0.50,  0.025, 'HXB2 index', transform=fig.transFigure, **tprops)
    ax_hist[-1].text(   0, -4, '6225', **tprops)
    ax_hist[-1].text(2645, -4, '8794', **tprops)

    plot.savefig('figures/figs9-cap256-recombination.pdf', dpi = 1000, facecolor = fig.get_facecolor(), edgecolor=None, **FIGPROPS)
    plot.close(fig)

    print('Figure S9 done.')


def plot_hive(ax_hive, tag, **pdata):
    """
    Hive plot of \Delta s effects.
    """

    # import stored data
    
    df_ds    = pd.read_csv('%s/analysis/%s-delta-s.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_poly  = pd.read_csv('%s/analysis/%s-poly.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    df_index = pd.read_csv('%s/processed/%s-index.csv' % (HIV_DIR, tag), comment='#', memory_map=True)
    L        = len(df_index)
    
    # set epitope labels

    epitope_start = []
    epitope_end   = []
    epitope_label = list(set(list(df_index[pd.notnull(df_index.epitope)].epitope)))
    print(tag, epitope_label)
    
    for e in epitope_label:
        epitope_start.append(np.min(df_index[df_index.epitope==e].alignment))
        epitope_end.append(np.max(df_index[df_index.epitope==e].alignment))

    # PLOT FIGURE

    # populate links

    idx_mask_pos  = []
    idx_delta_pos = []
    ds_pos        = []
    idx_mask_neg  = []
    idx_delta_neg = []
    ds_neg        = []

    c_cutoff  = 0.003
    ds_cutoff = 0.004
    ds_max    = 0.010

    df_ds = df_ds[np.fabs(df_ds.effect)>ds_cutoff]
    for it, entry in df_ds.iterrows():
        if (entry.mask_polymorphic_index==entry.target_polymorphic_index#):
#            continue
            and entry.mask_nucleotide==entry.target_nucleotide):
            continue
        if entry.effect>0:
            mask_alignment_index   = df_index[df_index.polymorphic==  entry.mask_polymorphic_index].iloc[0].alignment
            target_alignment_index = df_index[df_index.polymorphic==entry.target_polymorphic_index].iloc[0].alignment
            idx_mask_pos.append(mask_alignment_index)
            idx_delta_pos.append(target_alignment_index)
            ds_pos.append(entry.effect)
        elif entry.effect<0:
            mask_alignment_index   = df_index[df_index.polymorphic==  entry.mask_polymorphic_index].iloc[0].alignment
            target_alignment_index = df_index[df_index.polymorphic==entry.target_polymorphic_index].iloc[0].alignment
            idx_mask_neg.append(mask_alignment_index)
            idx_delta_neg.append(target_alignment_index)
            ds_neg.append(entry.effect)
    
    ds_pos           = (np.array(ds_pos) - c_cutoff) / ds_max
    ds_pos[ds_pos>1] = 1
    pos_sort         = np.argsort(ds_pos)[::-1]
    idx_mask_pos     = np.array(idx_mask_pos)[pos_sort]
    idx_delta_pos    = np.array(idx_delta_pos)[pos_sort]
    ds_pos           = ds_pos[pos_sort]

    ds_neg           = (np.fabs(np.array(ds_neg)) - c_cutoff) / ds_max
    ds_neg[ds_neg>1] = 1
    neg_sort         = np.argsort(ds_neg)[::-1]
    idx_mask_neg     = np.array(idx_mask_neg)[neg_sort]
    idx_delta_neg    = np.array(idx_delta_neg)[neg_sort]
    ds_neg           = ds_neg[neg_sort]

    # arc plot for large values of Delta s

    c_pos  = LCOLOR
    c_neg  = LCOLOR
    c_circ = { True : c_neg, False : c_pos }

    r_min     = 0.05
    r_max     = 1.00
    r_norm    = float(L) / (r_max - r_min)
    arc_mult  = SIZELINE * ds_max / 1.5
    arc_alpha = 1

    mask_angle_p =  -np.pi/2
    mask_angle_n = 3*np.pi/2
    ds_pos_angle =   np.pi/6
    ds_neg_angle = 5*np.pi/6
    da           = 0.005

    circ_color  = [hls_to_rgb(0.02, 0.53 * ds + 1. * (1 - ds), 0.83) for ds in ds_pos]
    circ_color += [hls_to_rgb(0.58, 0.53 * ds + 1. * (1 - ds), 0.60) for ds in ds_neg]
    circ_rad    = [[r_min + idx_mask_pos[i]/r_norm, r_min + idx_delta_pos[i]/r_norm] for i in range(len(ds_pos))]
    circ_rad   += [[r_min + idx_mask_neg[i]/r_norm, r_min + idx_delta_neg[i]/r_norm] for i in range(len(ds_neg))]
    circ_arc    = [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_pos]
    circ_arc   += [dict(lw=arc_mult * np.fabs(0.1) / ds_cutoff, alpha=arc_alpha) for ds in ds_neg]
    circ_angle  = [[mask_angle_p+da/r[0], ds_pos_angle-da/r[1]] for r in circ_rad[:len(ds_pos)]]
    circ_angle += [[mask_angle_n-da/r[0], ds_neg_angle+da/r[1]] for r in circ_rad[len(ds_pos):]]
    circ_bezrad = [(r[0]+r[1])/2 for r in circ_rad]
    circ_x      = [i for i in ds_pos] + [i for i in ds_neg]
    circ_y      = [i for i in ds_pos] + [i for i in ds_neg]
    
    # mark epitopes

    for i in range(len(epitope_label)):
        r_mid = r_min + (epitope_start[i] + epitope_end[i]) / (2 * r_norm)
        scatter_x = [r_mid * np.cos(mask_angle_p),
                     r_mid * np.cos(ds_pos_angle),
                     r_mid * np.cos(ds_neg_angle) ]
        scatter_y = [r_mid * np.sin(mask_angle_p),
                     r_mid * np.sin(ds_pos_angle),
                     r_mid * np.sin(ds_neg_angle) ]
        smallprops = dict(lw=0, marker='o', s=1.0*SMALLSIZEDOT, zorder=9999, clip_on=False)
        bigprops   = dict(lw=0, marker='o', s=1.7*SMALLSIZEDOT, zorder=9998, clip_on=False)
        mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=['#FFFFFF'], plotprops=smallprops)
        mp.scatter(ax=ax_hive, x=[scatter_x], y=[scatter_y], colors=[BKCOLOR],   plotprops=  bigprops)

    # Make plot

    pprops = { 'colors':   circ_color,
               'xlim':     [-1.0, 1.0],
               'ylim':     [-1.0, 1.0],
               'size':     L,
               'rad':      circ_rad,
               'arcprops': circ_arc,
               'angle':    circ_angle,
               'bezrad':   circ_bezrad,
               'noaxes':   True }

    plotprops = dict(lw=AXWIDTH, ls='-', zorder=999)
    line_x = [[r_min * np.cos(mask_angle_p), r_max * np.cos(mask_angle_p)],
              [r_min * np.cos(ds_pos_angle), r_max * np.cos(ds_pos_angle)],
              [r_min * np.cos(ds_neg_angle), r_max * np.cos(ds_neg_angle)] ]
    line_y = [[r_min * np.sin(mask_angle_p), r_max * np.sin(mask_angle_p)],
              [r_min * np.sin(ds_pos_angle), r_max * np.sin(ds_pos_angle)],
              [r_min * np.sin(ds_neg_angle), r_max * np.sin(ds_neg_angle)] ]
    mp.line(ax=ax_hive, x=line_x, y=line_y, colors=[BKCOLOR for i in range(len(line_x))], plotprops=plotprops)

    mp.plot(type='circos', ax=ax_hive, x=circ_x, y=circ_y, **pprops)
