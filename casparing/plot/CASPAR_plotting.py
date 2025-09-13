'''
Module with functionalities to assist with plotting different features of CASPAR

WIP!
'''

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as COLORS

def plot_cluster(dist=False):
    '''
    function to provide associations and distances in CASPAR and give them a unique marker for plotting.
    
    Parameters
    -------
    dist: bool, default=False
        boolean to include distances of each association
    
    Returns
    -------
    cluster: list
        list of association names
    marker: list
        list of markers for scatter plotting
    distance: list, (only if dist=True)
        list of distances for each assocation
    '''
    
    cluster = ['25 Orionis','118 Tau', 'Argus','B Pictoris', 'Chamaeleon I','Corona-Australis', 'IC 348', 'Lagoon Nebula', 
               'Lupus', 'NGC 2024','n Chamaeleontis','p Oph','Sh 2-284', 'Taurus', 'Tucana-Horologium','TW Hya',
            'Upper Centaurus Lupus', 'Upper Scorpius', 'σ Ori', 
            'Field']
    marker = ['o', '1', 'v', 'H','^', 's','2', 'p', 'P', '*', 'h','3', 'x', 'D', '4',
             'X', '<', 'd',  '8', '>']
    if dist:
        distance = [354, 101, 80, 40, 190, 130, 321, 1326, 158, 420, 94, 150, 135, 49, 51, 130, 141, 406, np.nan]
        return cluster, marker, distance
    else:
        return cluster, marker

    
def plot_age():
    '''
    function to provide markers, colors, and darker versions of each color for age bins [0-1, 1-3, 3-8, 8-45, 45+]
    
    Parameters
    -------
    NoneType
        None
    
    Returns
    -------
    age: list
        list of ages in Myr
    marker: list
        list of matplotlib.plot acceptable markers for each age
    c: list
        list of colors for each age
    darkc: list
        list of slightly darker colors for each age.
    '''
    
    age = [0,1, 3,8,45]
    marker = ['o', '*', 's', '^', 'D', '^']
    c = ['#FFB000', '#FE6100', '#DC267F', '#785EF0', '#648FFF', 'gray']
    darkc = ['#CC8D00', '#CC4E00', '#AA2761', '#5C48B7', '#336CFF', 'k']
    return age, marker, c, darkc

def plot_objtype_legend(ax, cS, cBD, cP, fs=12,loc='upper right', labels_color=True):
    '''
    function to create legend specifically when stars, brown dwarfs, and planetary mass companions are separate colors and marker types

    Parameters
    -------
    ax: matplotlib.axis.Axes
        single matplotlib ``Axis`` object 
    cS: str
        color of star datapoints when plotting
    cBD: str
        color of brown dwarf datapoints when plotting
    cP: str
        color of planetary mass companion datapoints when plotting
    fs: int, optional (default=12)
        fontsize of axis labels in main plot.  The legend labels will be fs/2.
    loc: str, optional (default='upper right')
        location of legend in main plot.  Use matplotlib.pyplot.legend loc locations.
    labels_color: bool, optinal (default=True)
        make the text the same color as the marker colors
    
    Returns
    -------
    NoneType
        None
    '''
    
    s1= ax.scatter([],[], color=cS, marker='o')
    s2= ax.scatter([],[], color=cBD,  marker='D')
    s3= ax.scatter([],[], color=cP,  marker='s')
    if labels_color:
        labelcolor=[cS, cBD, cP]
    else:
        labelcolor='k'
    leg = ax.legend([s1, s2, s3], ['Star', 'Brown Dwarf', 'Planetary Mass\nCompanion'], 
    ncol=1, fontsize=fs/2, loc=loc,title_fontsize=fs/2, markerscale=0.5, 
    labelcolor=labelcolor, frameon=True, columnspacing=0.3,
    handletextpad=0.1) 
    leg.set_zorder(1000)
      
def plot_createSinglePlot(figsize=(4,3)):
    '''   
    function to create a single matplotlib.pyplot figure and axis.
    
    Parameters
    -------
    figsize: tuple(float, float), optional (default=(4,3))
        figure size.  Should follow matplotlib.pyplot figsize convention.
    
    Returns
    -------
    matplotlib.figure.Figure
        The ``Figure`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The ``Axis`` object that can be used for further
        customization of the plot.
    '''
    
    fig = plt.figure(figsize=figsize, dpi=300)
    ax = fig.add_subplot(111)
    return fig,ax

def plot_create6multiplots(cols=2):
    '''   
    function to create a six panel matplotlib.pyplot figure and axis.
    
    Parameters
    -------
    cols: int, optional (default=2)
        number of columns in figure.  rows will be 6/cols.  Must be either 2 or 3.
    
    Returns
    -------
    matplotlib.figure.Figure
        The ``Figure`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The first ``Axis`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The second ``Axis`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The third ``Axis`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The forth ``Axis`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The fifth ``Axis`` object that can be used for further
        customization of the plot.
    matplotlib.axis.Axes
        The sixth ``Axis`` object that can be used for further
        customization of the plot.
    '''
    if cols==2:
        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(figsize=(4.5,7), nrows=3,ncols=2,dpi=150, sharey=True)
    elif cols==3:
        fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(figsize=(6.4,3.7), nrows=2,ncols=3,dpi=150, sharey=True)

    else:
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(figsize=(11,2.6), nrows=1,ncols=6,dpi=300, sharey=True)
    return fig, ax1, ax2, ax3, ax4, ax5, ax6
    
def plot_tickparams(ax, fs=12, minor='on'):   
    ax.tick_params(which='major',length=7, direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
    ax.tick_params(which='minor',length=4,  direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
    if minor=='on':
        ax.minorticks_on()
    else:
        ax.minorticks_off()

def plot_MMdotFeatures(fig, ax, fs=12, burning_lines=True, burning_lines_labels=True, tickparams=True, **kwargs):
    if burning_lines_labels:
        fs_BL=kwargs.get('fs_BL', 5) 
        ax.text(-1.124938736, 0.95, 'Hydrogen\nBurning Limit', color='dimgray', rotation=90, fontsize=fs_BL, ha='center', transform=ax.get_xaxis_transform(), va='top')

        ax.text(-1.906371, 0.95, 'Deuterium\nBurning Limit', color='dimgray', rotation=90, fontsize=fs_BL, ha='center', transform=ax.get_xaxis_transform(), va='top')
        
    if burning_lines:
        ax.axvline(x = -1.124938736, color = 'dimgray', linestyle='-', alpha=0.8, linewidth=0.75)
        ax.axvline(x = -1.906371, color = 'dimgray', linestyle='-', alpha=0.8, linewidth=0.75)
    
    xlabel = kwargs.get('xlabel', 'log Mass ($M_{\odot}$)')
    ax.set_xlabel(xlabel, fontsize = fs)
    
    ylabel = kwargs.get('ylabel', 'log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)')
    ax.set_ylabel(ylabel, fontsize = fs)
    
    if tickparams:
        fs_ticks = kwargs.get('fs_ticks', fs)
        ax.tick_params(which='major',length=7, direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs_ticks)
        ax.tick_params(which='minor',length=4,  direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs_ticks)
        ax.minorticks_on()
    
    xlim = kwargs.get('xlim', (-2.5, 0.3))
    ylim = kwargs.get('ylim', (-13.95, -5))
    ax.set_xlim(xlim)    
    ax.set_ylim(ylim)
    if fig is not None:
        fig.tight_layout()
        fig.subplots_adjust(hspace=0)


def plot_MMdot(ax, XXreg, YYreg, XXupp, YYupp, color_cmap=False, **kwargs):
    if 'marker' in kwargs:
        marker = kwargs['marker']
    else:
        marker= 'o'
    if 's' in kwargs:
        s = kwargs['s']
    else:
        s=15
    if 'edgecolor' in kwargs:
        edgecolor = kwargs['edgecolor']
    else:
        edgecolor=None

    if 'linewidth' in kwargs:
        linewidth=kwargs['linewidth']
    else:
        linewidth=1
    if 'alpha' in kwargs:
        alpha = kwargs['alpha']
    else:
        alpha=1
    if 'color' in kwargs:
        color=kwargs['color']
    else:
        color='k'
    if 'zorder' in kwargs:
        zorder=kwargs['zorder']
    else:
        zorder=5
    if 'uppcolor' in kwargs:
        uppcolor=kwargs['uppcolor']
    else:
        uppcolor=color
    u = np.zeros_like(XXupp.values)
    v = -np.ones_like(XXupp.values)
    if color_cmap==False:
        pp = ax.scatter(XXreg, YYreg,  c=color, marker=marker, s=s, facecolor=color, edgecolor=edgecolor, linewidth=linewidth, alpha=alpha, zorder=zorder)

        ax.quiver(XXupp.values, YYupp.values, u,v, headwidth=4, headaxislength=3,headlength=3, width=.005, scale=30, alpha=alpha, color=uppcolor, zorder=zorder, edgecolor=edgecolor, linewidth=linewidth)
    else:
        if 'cmap' in kwargs:
            cmap=kwargs['cmap']
        else:
            cmap=cm.jet
        if 'vmin' in kwargs:
            vmin = kwargs['vmin']
        else:
            vmin=35
        if 'vmax' in kwargs:
            vmax=kwargs['vmax']
        else:
            vmax=1350
        
        pp = ax.scatter(XXreg, YYreg, s=s, c=color,
             marker=marker, cmap=cmap,linewidth=linewidth,
            edgecolor=edgecolor, alpha=alpha, norm=COLORS.LogNorm(vmin=vmin, vmax=vmax))
        

        if len(XXupp.values)>0:
            print(len(XXupp.values))
            ax.quiver(XXupp.values, YYupp.values, u,v,  uppcolor, headwidth=4, headaxislength=3,headlength=3, width=.005, scale=30, alpha=alpha, norm=COLORS.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, edgecolor=edgecolor, linewidth=linewidth)
    return pp
        
    
def plot_graylines(ax, df_caspar):
    df_casparCOM = df_caspar.loc[(df_caspar['Mass']<=0.075) & (df_caspar['Companion']=='COM')]
    NAMES = df_casparCOM['Unique Name'].unique()
    for name in NAMES:
        dfname = df_casparCOM.loc[df_casparCOM['Unique Name']==name]
        XX = dfname['log Mass'].values
        XX = np.append(XX, XX[0])
        YY = dfname['log Mdot'].values
        YY = np.append(YY,YY[0])
        ax.plot(XX, YY, '--', color='gray', lw=0.5, zorder=-1)
    
    df_caspar2 = df_caspar.loc[(df_caspar['Companion']!='COM')]
    NAMES = df_caspar2['Unique Name'].unique()
    for name in NAMES:
        dfname = df_caspar2.loc[df_caspar2['Unique Name']==name]
        
        ax.plot(dfname['log Mass'], dfname['log Mdot'], '--', color='gray', lw=0.5, zorder=-1)
        
def plot_litlines(ax):
    x = np.linspace(-5.5,-1.12,10)
    y=0.12*x-10.5
    line1, = ax.plot(x,y,color='k',ls='--', zorder=4, linewidth=1)
    ax.fill_between(x, y-0.3, y+0.3, color='gray', alpha=0.4, zorder=-1)
    
    x = np.arange(np.log10(0.075),1,0.1)
    y2=2.1*x-7.9
    line2, = ax.plot(x,y2,color='k', zorder=4, linewidth=1)
    ax.fill_between(x, y2-0.5, y2+0.5, color='gray', alpha=0.2, zorder=100)
    
    x = np.arange(-5.5,np.log10(0.075)+0.1,0.1)
    y2=2.1*x-7.9
    ax.plot(x,y2,color='gray', alpha=0.8, zorder=-3, linewidth=1)
    ax.fill_between(x, y2-0.5, y2+0.5, color='gray', alpha=0.1, zorder=100)
    return line1, line2
    

    
    
    
    
    
