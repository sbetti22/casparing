from scipy import interpolate
import matplotlib.patheffects as pe
import matplotlib.gridspec as gridspec
from sklearn.metrics import r2_score
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import .CASPAR_util as cutil
import .CASPAR_sortdata as csort
import .CASPAR_plotting as cplot
import .CASPAR_fitMMdot as cfit
# get dictionary of accretion diagnosticas
accdiaglines = cutil.AccDiagLines()

colors = cutil.AccDiagColor()
darkcolors = cutil.AccDiagColorDark()
linfit = cutil.linfitMMdot()


def markerlegend(ax, cS, cBD, cP, fs = 10, loc='upper left'):
    s1= ax.scatter([],[], color=cS, marker='o', edgecolor='k', linewidth=0.5)
    s2= ax.scatter([],[], color=cBD,  marker='D', edgecolor='k', linewidth=0.5)
    s3= ax.scatter([],[], color=cP,  marker='s', edgecolor='k', linewidth=0.75)
    second_leg = ax.legend([s1, s2, s3], ['Star', 'Brown Dwarf', 'Planetary Mass\nCompanion'], 
    ncol=1, fontsize=fs, loc=loc,title_fontsize=fs, frameon=True, columnspacing=0.3,
    handletextpad=0.1)  
    second_leg.set_zorder(1000)
    
def plotfeatures(ax, linetype='CASPAR all', **kwargs):
    ax.annotate('Hydrogen\nBurning Limit', (-1.23, -13.5), color = 'dimgray', rotation = 90, fontsize = 2, ha='left')
    ax.annotate('Deuterium\nBurning Limit', (-2.008, -13.5), color = 'dimgray', rotation = 90, fontsize = 2, ha='left')
    ax.axvline(x = np.log10(0.075), color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
    ax.axvline(x = np.log10(13 * .0009543), color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)

    x = np.arange(np.log10(0.075),1,0.1)

    y2 = linfit.at[linetype, 'slope']*x + linfit.at[linetype, 'intercept']
    C1 = kwargs.get('C', 'k')
    line2, = ax.plot(x,y2,color=C1, label="collapsing prestellar cores", zorder=4, linewidth=1)
    x = np.arange(-5.5,np.log10(0.075)+0.1,0.1)
    y2 = linfit.at[linetype, 'slope']*x + linfit.at[linetype, 'intercept']

    C2 = kwargs.get('C', 'gray')
    line2, = ax.plot(x,y2,color=C2, alpha=1, zorder=-3, linewidth=1)

    ax.set_xlim(-2.5, 0.5)    
    ax.set_ylim(-13.95, -5)

    ax.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True)
    ax.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True)
    ax.minorticks_on()
    
def MdotMplot(ax, a, b, diag, **kwargs):
    Areg, Aupp = csort.CASPAR_separateByUppLimit(a)
    Breg, Bupp = csort.CASPAR_separateByUppLimit(b)
    cplot.plot_MMdot(ax, Areg['log Mass'], Areg['log Mdot'], Aupp['log Mass'], Aupp['log Mdot'], **kwargs)
    cplot.plot_MMdot(ax, Breg['log Mass'], Breg['log Mdot'], Bupp['log Mass'], Bupp['log Mdot'], **kwargs)
        
def MdotbyDiagnostic(df_caspar, cols=3,fit_data=True, fig=None, ax=None,title=True, **kwargs):
    '''
    kwargs: labelC - boolean
    marker
    s
    edgecolor
    linewidth
    alpha
    lines
    binn
    fit_mass_regions
    fitname
    '''
    fit_mass_regions= kwargs.get('fit_mass_regions')
    
    fitname=kwargs.get('fitname','doLinearFit')
    if fig is None:
        fig, (ax1, ax2, ax3,ax4) = plt.subplots(figsize=(8,2.4), nrows=1,ncols=4,dpi=150, sharey=True)
        ax = [ax1, ax2, ax3, ax4]
    else:
        ax1, ax2, ax3, ax4 = ax
    
    if kwargs.get('labelC'):
        fig.text(0.01, 0.94, 'c)', fontsize=20)

    for i, diag in enumerate(np.sort(list(colors.keys()))):
        print(diag)
        plotfeatures(ax[i])
        
        
        a = df_caspar.loc[df_caspar['AD']==diag]
        reg, upp = csort.CASPAR_separateByUppLimit(a)
        cplot.plot_MMdot(ax[i], reg['log Mass'], reg['log Mdot'], upp['log Mass'], upp['log Mdot'], 
                             color=colors[diag], **kwargs)
    
        ALL_MASSES = a['log Mass'].to_numpy() 
        ALL_MASSESerr = a['log Mass err'].to_numpy() 
        ALL_MDOTS = a['log Mdot'].to_numpy()
        ALL_MDOTSerr = a['log Mdot err'].to_numpy() 
        ALL_UPPLIM = a['Upper Limit bool'].to_numpy() 
        
        if diag=='Hα 10% width':  
            line = accdiaglines[diag]
            b = df_caspar.loc[df_caspar[f'{line} Accretion Rate'].notnull()]
            breg, bupp = csort.CASPAR_separateByUppLimit(b)
            b[f'{line} log Mdot err'] = np.nanmedian(a['log Mdot err'])
            b[f'{line} Upper Limit bool'] = b[f'{line} Upper Limit'] != 'UPP' # true is detected, false is upper limit
            ALL_MASSES = np.append(ALL_MASSES, b['log Mass'].values)
            print(ALL_MASSES)
            ALL_MASSESerr = np.append(ALL_MASSESerr, b['log Mass err'].values)
            ALL_MDOTS = np.append(ALL_MDOTS, np.log10(b[f'{line} Accretion Rate']).values)
            ALL_MDOTSerr = np.append(ALL_MDOTSerr, b[f'{line} log Mdot err'].values)
            ALL_UPPLIM = np.append(ALL_UPPLIM, b[f'{line} Upper Limit bool'].values)
            cplot.plot_MMdot(ax[i], breg['log Mass'], np.log10(breg[f'{line} Accretion Rate']), bupp['log Mass'],
                                 np.log10(bupp[f'{line} Accretion Rate']),
                                 color=colors[diag], **kwargs)
        CC=colors[diag]
        dCC= darkcolors[diag]
        if kwargs.get('lines'):
            if fitname == 'original':
                if 'α' in diag:
                    diagname = diag.replace('α', 'a')
                else: diagname=diag
                finfitname = 'original ' + diagname
            else:
                finfitname = fitname
            XX, YY, XXerr, YYerr, UPP = ALL_MASSES, ALL_MDOTS, ALL_MASSESerr, ALL_MDOTSerr, ALL_UPPLIM
            
            if fit_mass_regions:
                mass_bounds = [-2.5, -1.1249, 0.5]
            else:
                mass_bounds = [-2.5, 0.5]

            
            for k in np.arange(len(mass_bounds)-1):
                ind = np.where((XX > mass_bounds[k]) & (XX <= mass_bounds[k+1]))
                if len(ind[0]) >0:
                    X2, Y2, Xerr2, Yerr2, upp2 = XX[ind], YY[ind], XXerr[ind], YYerr[ind], UPP[ind]
                    idx = np.isfinite(X2) & np.isfinite(Y2)
                    Xfin, Yfin, Xerrfin, Yerrfin, uppfin = X2[idx], Y2[idx], Xerr2[idx], Yerr2[idx], upp2[idx]
                    df_vals = pd.DataFrame({'log Mass':Xfin, 'log Mass err':Xerrfin, 'log Mdot':Yfin, 'log Mdot err':Yerrfin, 'Upper Limit bool':uppfin})
                    
                    result = cfit.linmix_fitting(df_vals, fitname=finfitname, xval='logMass')
                    N, m, slopeerrlow, slopeerrupp, b, intercepterrlow, intercepterrupp, stdfit, r2=result
                    print(diag, mass_bounds[k], result)
                    
                    x = np.linspace(mass_bounds[k], mass_bounds[k+1], 50)
                    if (diag == 'Veiling') and (k == 0):
                        X = x[20]
                    else:
                        X = x[3]
                    ax[i].text(X, (m*X+b)+0.3, 'log$\dot M$='+str(round(m,1))+'log$M$'+str(round(b,1)), fontsize=4, zorder=101,
                                color='k', rotation=np.rad2deg(np.arctan(m)), rotation_mode='anchor', transform_rotates_text=True,)


                    ax[i].plot(x, m*x + b, '-', color=dCC, zorder=100, linewidth=2)
        else: ## binns
            if fit_mass_regions:
                XX, YY = a['log Mass'].values, a['log Mdot'].values
                XXcom, YYcom = a['log Mass'].values, a['log Mdot'].values
                mass_bounds = [-2.5, -1.1249, 0.5]
                mark= ['D', 'o']
                binns = [3,7]
            else:
                XX, YY = ALL_MASSES, ALL_MDOTS
                mass_bounds = [-5, 0.5]
                mark =['o']
                binns=[7]

            for k in np.arange(len(mass_bounds)-1):
                print(mark[k])
                ind = np.where((XX > mass_bounds[k]) & (XX <= mass_bounds[k+1]))
                if len(ind[0]) >0:
                    X2, Y2 = XX[ind], YY[ind]
                    idx = np.isfinite(X2) & np.isfinite(Y2)
                    cutil.binn(X2[idx], Y2[idx], nbins=binns[k], plot=True, ax=ax[i], color=CC, marker=mark[k])

        if title:
            if diag == 'Hα Photometric Luminosity':
                diag = 'Hα Photometric\nLuminosity'
            fs = kwargs.get('fs', 11)
            ax[i].text(0.05, 0.95, diag, fontsize=fs, ha='left', va='top', transform=ax[i].transAxes)
    ax1.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)')
    for axi in ax:
        axi.set_xlabel('log $M\ (M_{\odot})$')
    if cols==2:
        ax4.set_xlabel('')
        ax3.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)')
        ax5.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)')
    elif cols==3:
        ax4.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)')
    else:
        ax1.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)')
    plt.subplots_adjust(wspace=0.03)
    plt.tight_layout()
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    if fig is None:
        plt.show()
    

def MdotVsAccDiagresiduals(df_caspar, **kwargs):
    fs = kwargs.get('fs', 8)
    fig, ax = cplot.plot_createSinglePlot(figsize=(4,2))
    ax.axvline(x = -1.124938736, color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
    ax.axvline(x = -1.906371, color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
    df = df_caspar.sort_values(by='AD', ascending=True)
    AD = np.flip(df_caspar['AD'].unique())
        
    mark = 'o'
    bins = 8

    DF = df_caspar
    VALS = np.sort(DF['AD'].unique())
    for AD in VALS:
        
        if (AD =='Continuum Excess') | (AD == 'Line Luminosity'):
            if AD == 'Line Luminosity':
                col = 'darkred'
            else:
                col=colors[AD]
            x = DF['log Mass'].loc[DF['AD']==AD].values
            y = DF['log Mdot'].loc[DF['AD']==AD].values
 
            XX, YY, YYerr = cutil.binn(x, y, nbins=bins, plot=False)

            s = np.linspace(5,40,len(XX))
            x = np.linspace(np.min(XX), np.max(XX), 500)
            y = linfit.at['CASPAR all', 'slope']*x + linfit.at['CASPAR all', 'intercept']
            f = interpolate.interp1d(x, y)
            ax.errorbar(XX, YY-f(XX), yerr=YYerr, color=col, ecolor=col,lw=0.7, zorder=-1, fmt='None')
            ax.scatter(XX, YY-f(XX), c=col, marker=mark,  s=50, edgecolor='k',linewidth=0.7, alpha=1)

    ax.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
    ax.tick_params(which='minor', direction='in', top=False,bottom=False, right=True, left=True, labelsize=fs)
    ax.minorticks_on()
    ax.set_ylim(-1.5, 3.1)
    ax.axhline(0, color='k', zorder=-1)

    plt.subplots_adjust(bottom=0.38, top=0.95, left=0.18)
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    plt.show()    
    
def MdotVsAccDiagAll(df_caspar,  fig=None, ax=None,**kwargs):
    fs = kwargs.get('fs',8)
    if fig is None:
        fig, ax = cplot.plot_createSinglePlot()
    
    df = df_caspar.sort_values(by='AD', ascending=True)
    AD = np.flip(df_caspar['AD'].unique())
        
    mark = 'o'
    bins = 7

    valpos = {}
    labels = list(colors.keys())
    for i, name in enumerate(labels):
        valpos[name] = i+1

    VALS = np.sort(df_caspar['AD'].unique())
    XXfin, YYfin, YYerrfin = [], [],[]
    for AD in VALS:
        x = df_caspar['log Mass'].loc[df_caspar['AD']==AD].values
        y = df_caspar['log Mdot'].loc[df_caspar['AD']==AD].values
        XX, YY, YYerr = cutil.binn(x, y, nbins=bins, plot=False)
        XXfin.append(XX)
        YYfin.append(YY)
        YYerrfin.append(YYerr)
    
    for i in np.arange(len(VALS)):
        s = np.linspace(5,40,len(XXfin[i]))
        x = np.linspace(np.min(XXfin[i]), np.max(XXfin[i]), 500)
        y = linfit.at['CASPAR all', 'slope']*x + linfit.at['CASPAR all', 'intercept']

        f = interpolate.interp1d(x, y)
        offset = np.linspace(valpos[VALS[i]]-0.3, valpos[VALS[i]]+0.3, len(YYfin[i]))
        ax.errorbar(offset, YYfin[i]-f(XXfin[i]), yerr=YYerrfin[i], color=colors[VALS[i]], ecolor=colors[VALS[i]],lw=0.7, zorder=-1, fmt='None')
        ax.scatter(offset, YYfin[i]-f(XXfin[i]), c=colors[VALS[i]], marker=mark,  s=10+s, edgecolor='k',linewidth=0.7, alpha=1)
        print(VALS[i],YYfin[i]-f(XXfin[i]))

    for j in np.arange(len(XXfin[i])):
        ax.scatter([],[], c='k', marker=mark,  s=10+s[j], edgecolor='k',linewidth=0.7, alpha=1, 
                   label=round(10**XXfin[i][j],3))
    ax.legend(fontsize=5, ncol=1, labelspacing=0.6, handletextpad=0, columnspacing=0, framealpha=1, title='Mass ($M_\odot$)', title_fontsize=5)
    labels[1] = 'Hα Photometric\nLuminosity'
    ax.set_xticks(np.arange(1,len(labels)+1))
    ax.set_xticklabels(labels, rotation=90, fontsize=fs)
    ax.set_xlim(0.5,5.5)
    ax.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
    ax.tick_params(which='minor', direction='in', top=False,bottom=False, right=True, left=True, labelsize=fs)
    ax.minorticks_on()
    ax.set_ylim(-2.75, 3.5)
    ax.axhline(0, color='k', zorder=-1)
    ax.set_ylabel('$\dot M$ residual\n($\dot M_\mathrm{bin}-\dot M_\mathrm{fit}$; dex)', fontsize=fs)
#    fig.tight_layout()
    plt.subplots_adjust(bottom=0.38, top=0.95, left=0.18)
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    if fig is None:
        plt.show()
    
    


def Ha10VSFlux(df_caspar, **kwargs):
    IND = (~df_caspar['Ha 10% Accretion Rate'].isna()) & (~df_caspar['Ha Accretion Rate'].isna())

    Ha10 = df_caspar['Ha 10% Accretion Rate'][IND].values
    Ha = df_caspar['Ha Accretion Rate'][IND].values
    Haerr = df_caspar['Ha Accretion Rate err'][IND].values
    d = df_caspar['Distance'][IND].values
    M = df_caspar['Mass'][IND].values
    age = df_caspar['Age'][IND].values
    COM = df_caspar['Companion'][IND].values
    yerr = 0.434 * Ha/Haerr
   
    fig, (ax1, ax2) = plt.subplots(figsize=(6,8), ncols=1, nrows=2)
    for i in np.arange(len(Ha10)):
        if M[i] > 0.075:
            marker='o'
        elif (M[i] <= 0.075) & (COM[i] == 'COM'):
            marker='s'
        else:
            marker='D'

        ax2.scatter(np.log10(M[i]), np.log10(Ha10[i]/Ha[i]), c=np.log10(M[i]), marker=marker,vmin=-2, vmax=-0.5, ec='k', s=age[i]*10,zorder=100, cmap='turbo')

        ax1.errorbar(np.log10(Ha10[i]), np.log10(Ha[i]), yerr=yerr[i], ecolor='k', capsize=3)
        im = ax1.scatter(np.log10(Ha10[i]), np.log10(Ha[i]), c=np.log10(M[i]), s=age[i]*10, marker=marker,vmin=-2, vmax=-0.5, zorder=100, ec='k', cmap='turbo')
        
    cplot.plot_MMdotFeatures(None, ax1, 
                             fs=16, fs_BL=10, fs_ticks=16, 
                             burning_lines=False, tickparams=True,
                             xlabel=r'log[$\dot{M}$ (Hα 10%)]',
                             ylabel=r'log[$\dot{M}$ (Hα Lum)]', 
                             xlim=(-13.5,-6), ylim=(-13.5,-6))
    
    cplot.plot_MMdotFeatures(None, ax2, 
                             fs=16, fs_BL=10, fs_ticks=16, 
                             burning_lines=True, tickparams=True,
                             xlabel='log Mass (M$_\odot$)',
                             ylabel=r'log[$\dot{M}$ (Hα 10%) / $\dot{M}$ (Hα Lum)] (dex)', 
                             xlim=(-2.5, 0.3), ylim=(-2,2.5))
    
    ax1.plot([-13.5,-6], [-13.5,-6], color='k')
    for ax in [ax1, ax2]:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('log Mass (M$_\odot$)',fontsize=16)
        cbar.ax.tick_params(labelsize=16, direction='in', length=7)
    plt.tight_layout()
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    plt.show()
    
def Ha10VSexcess(df_caspar, **kwargs):
    
    IND = (~df_caspar['Ha 10%'].isna()) & (df_caspar['Accretion Diagnostic']=='Continuum Excess')
    
    Ha10 = df_caspar['Ha 10%'][IND].values
    UV = df_caspar['Accretion Rate'][IND].values
    UVerr = df_caspar['Accretion Rate err'][IND].values
    UVlogerr = 0.434 * UV/UVerr
    
    d = df_caspar['Distance'][IND].values
    M = df_caspar['Mass'][IND].values
    age = df_caspar['Age'][IND].values
    COM = df_caspar['Companion'][IND].values

    fig, (ax1, ax2) = plt.subplots(figsize=(6,8), ncols=1, nrows=2)
    for i in np.arange(len(Ha10)):
        if M[i] > 0.075:
            marker='o'
        elif (M[i] <= 0.075) & (COM[i] == 'COM'):
            marker='s'
        else:
            marker='D'
        ax2.scatter(np.log10(M[i]), np.log10(Ha10[i]/UV[i]), c=np.log10(M[i]), marker=marker,vmin=-2, vmax=-0.5, ec='k', s=age[i]*10,zorder=100, cmap='turbo')

        ax1.errorbar(Ha10[i], np.log10(UV[i]), yerr=UVlogerr[i], ecolor='k', capsize=3)
        im = ax1.scatter(Ha10[i], np.log10(UV[i]), c=np.log10(M[i]), s=age[i]*10, marker=marker,vmin=-2, vmax=-0.5, zorder=100, ec='k', cmap='turbo')
        
    cplot.plot_MMdotFeatures(None, ax2, 
                             fs=16, fs_BL=10, fs_ticks=16, 
                             burning_lines=True, tickparams=True,
                             ylabel=r'log[$\dot{M}$ (Hα) / $\dot{M}$ (UV)] (dex)', 
                             xlim=(-2.5, 0.3), ylim=(10, 13))
    
    cplot.plot_MMdotFeatures(None, ax1, 
                             fs=16, fs_BL=10, fs_ticks=16, 
                             burning_lines=False, tickparams=True,
                             xlabel=r'log[$\dot{M}$ (Hα)]', ylabel=r'log[$\dot{M}$ (UV)]', 
                             xlim=(150, 600), ylim=(-13.5,-6))

    ax1.plot(np.array([150, 600]), 9.7e-3*np.array([150, 600]) - 12.89)
    for ax in [ax1, ax2]:
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('log Mass (M$_\odot$)',fontsize=16)
        cbar.ax.tick_params(labelsize=16, direction='in', length=7)

    plt.tight_layout()
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    plt.show()
    
    
def MdotbyNIRLineDiagnostic(df_caspar, fit_data=True, **kwargs):
    fs = kwargs.get('fs', 14)
    fit_mass_regions=kwargs.get('fit_mass_regions')
    fitname=kwargs.get('fitname', 'doLinearFit')
    
    fig = plt.figure(figsize=(10,8))
    gs0 = gridspec.GridSpec(2, 2, figure=fig)

    gs00 = gs0[0].subgridspec(3,2,wspace=0, hspace=0, height_ratios=[1,4,2], width_ratios=[4,1])
    ax1top = fig.add_subplot(gs00[0,0])
    ax1 = fig.add_subplot(gs00[1,0])
    ax2 = fig.add_subplot(gs00[2,0])
    ax1side = fig.add_subplot(gs00[1,1])
    ax2side = fig.add_subplot(gs00[2,1])

    gs01 = gs0[1].subgridspec(3,2,wspace=0, hspace=0, height_ratios=[1,4,2], width_ratios=[4,1])
    ax3top = fig.add_subplot(gs01[0,0])
    ax3 = fig.add_subplot(gs01[1,0])
    ax4 = fig.add_subplot(gs01[2,0])
    ax3side = fig.add_subplot(gs01[1,1])
    ax4side = fig.add_subplot(gs01[2,1])

    gs02 = gs0[2].subgridspec(3,2,wspace=0, hspace=0, height_ratios=[1,4,2], width_ratios=[4,1])
    ax5top = fig.add_subplot(gs02[0,0])
    ax5 = fig.add_subplot(gs02[1,0])
    ax6 = fig.add_subplot(gs02[2,0])
    ax5side = fig.add_subplot(gs02[1,1])
    ax6side = fig.add_subplot(gs02[2,1])
    
    gs03 = gs0[3].subgridspec(3,2,wspace=0, hspace=0, height_ratios=[1,4,2], width_ratios=[4,1])
    ax7top = fig.add_subplot(gs03[0,0])
    ax7 = fig.add_subplot(gs03[1,0])
    ax8 = fig.add_subplot(gs03[2,0])
    ax7side = fig.add_subplot(gs03[1,1])
    ax8side = fig.add_subplot(gs03[2,1])
    
    def plot(ax):
        cplot.plot_MMdotFeatures(fig, ax, fs=fs, fs_BL=5)

        x = np.arange(np.log10(0.075),1,0.1)
        y2 = linfit.at['CASPAR all', 'slope']*x + linfit.at['CASPAR all', 'intercept']
        line2, = ax.plot(x,y2,color='k', label="collapsing prestellar cores", zorder=4, linewidth=1)
        x = np.arange(-5.5,np.log10(0.075)+0.1,0.1)
        y2 = linfit.at['CASPAR all', 'slope']*x + linfit.at['CASPAR all', 'intercept']
        line2, = ax.plot(x,y2,color='gray', alpha=1, zorder=-3, linewidth=1)

    ax = [ax1, ax3, ax5, ax7]
    axres = [ax2, ax4, ax6, ax8]
    axtop = [ax1top, ax3top, ax5top, ax7top]
    axside = [ax1side, ax3side, ax5side, ax7side]
    axsideres= [ax2side, ax4side, ax6side, ax8side]
    CC = ['red', 'orange', 'green','blue']
    letter_label=['a) ', 'b) ', 'c) ', 'd) ']
        
    for i, diag in enumerate(['Balmer', 'Paschen/Brackett/Pfund', r'He $\mathtt{I}$', r'Ca $\mathtt{II}$']): #
        if fitname == 'original':
            if 'α' in diag:
                diagname = diag.replace('α', 'a')
            else: diagname=diag
            finfitname = 'original ' + diagname
        elif fitname == 'CASPAR':
            if '$\mathtt' in diag:
                diagname = diag.replace(' $\mathtt{', '').replace('}$', '')
            else: diagname=diag
            finfitname = 'CASPAR ' + diagname
        else:
            finfitname = fitname
                
        ALL_MASSES = np.array([])
        ALL_MDOTS = np.array([])
        ALL_MASSESerr = np.array([])
        ALL_MDOTSerr = np.array([])
        ALL_UPP = np.array([])
        for line in accdiaglines[diag]:
            df_caspar[f'{line} log Mdot'] = np.log10(df_caspar[f'{line} Accretion Rate'])
            df_caspar[f'{line} log Mdot err'] = 0.434 * df_caspar[f'{line} Accretion Rate err']/df_caspar[f'{line} Accretion Rate']
            df_caspar[f'{line} Upper Limit bool'] = df_caspar[f'{line} Upper Limit'] != 'UPP'
            
            a = df_caspar.loc[(df_caspar[f'{line} Accretion Rate'].notnull()&(df_caspar[f'{line} Upper Limit']!='UPP'))]
            reg, upp = csort.CASPAR_separateByUppLimit(a)
            
            cplot.plot_MMdot(ax[i], reg['log Mass'], reg[f'{line} log Mdot'], 
                                 upp['log Mass'], upp[f'{line} log Mdot'],
                                 color=CC[i], alpha=0.15, s=25)
            
            ALL_MASSES = np.append(ALL_MASSES, a['log Mass'].values)
            ALL_MASSESerr = np.append(ALL_MASSESerr, a['log Mass err'].values)
            ALL_MDOTS = np.append(ALL_MDOTS, a[f'{line} log Mdot'].values)
            ALL_MDOTSerr = np.append(ALL_MDOTSerr, a[f'{line} log Mdot err'].values)
            ALL_UPP = np.append(ALL_UPP, a[f'{line} Upper Limit bool'].values)

        medmdoterrorig = np.nanmedian(ALL_MDOTSerr)
        bad = np.where(np.isnan(ALL_MDOTSerr))
        ALL_MDOTSerr[bad] = medmdoterrorig
        
        axtop[i].hist(ALL_MASSES, color=CC[i], bins=20, histtype='stepfilled', edgecolor=CC[i], linewidth=2, alpha=0.5)
        
        axside[i].hist(ALL_MDOTS, color=CC[i], orientation='horizontal', bins=20, histtype='stepfilled', edgecolor=CC[i], linewidth=2, alpha=0.5)

        if fit_data:
            if fit_mass_regions:
                mass_bounds = [-2.5, -1.1249, 0.5]
            else:
                mass_bounds = [-2.5, 0.5]
            for k in np.arange(len(mass_bounds)-1):
                ind = np.where((ALL_MASSES > mass_bounds[k]) & (ALL_MASSES <= mass_bounds[k+1]))
                if len(ind[0]) >0:
                    X2, Y2, X2err, Y2err, upp2 = ALL_MASSES[ind], ALL_MDOTS[ind],ALL_MASSESerr[ind], ALL_MDOTSerr[ind], ALL_UPP[ind]
                    idx = np.isfinite(X2) & np.isfinite(Y2)
                    Xfin, Yfin, Xfinerr, Yfinerr, uppfin = X2[idx], Y2[idx],X2err[idx], Y2err[idx], upp2[idx]
                    df_vals = pd.DataFrame({'log Mass':Xfin, 'log Mass err':Xfinerr, 'log Mdot':Yfin, 'log Mdot err':Yfinerr, 'Upper Limit bool':uppfin})


                    result = cfit.linmix_fitting(df_vals, fitname=finfitname, xval='logMass')
                    N, m, slopeerrlow, slopeerrupp, b, intercepterrlow, intercepterrupp, stdfit, r2=result
                    print(diag, mass_bounds[k], result)

                    x = np.linspace(mass_bounds[k], mass_bounds[k+1], 50)

                    ax[i].plot(x, m*x + b, '-', color='dark'+CC[i], zorder=100, linewidth=3)

                    axres[i].scatter(X2, Y2-(m*X2 + b), s=25, color=CC[i], alpha=0.15)
                    axsideres[i].hist(np.array(Y2-(m*X2 + b)), color=CC[i], orientation='horizontal', bins=20, alpha=0.5, linewidth=2, histtype='stepfilled', edgecolor=CC[i])

                    x_main1 = np.arange(-5.5,0.5,0.1)
                    y_main1 = linfit.at['CASPAR all', 'slope']*x_main1 + linfit.at['CASPAR all', 'intercept']
                    axres[i].plot(x_main1, y_main1-(m*x_main1 + b),color='k', alpha=1, zorder=-3, linewidth=1)
                    axres[i].axhline(0, color='dark'+CC[i], zorder=-1, linewidth=3)
            plot(ax[i])
            ax[i].text(0.02, 0.95, letter_label[i]+diag, fontsize=fs+4, ha='left', va='top', transform=ax[i].transAxes)

            ind = np.where((ALL_MASSES > -1.1249))
            X2, Y2, X2err, Y2err, upp2 = ALL_MASSES[ind], ALL_MDOTS[ind],ALL_MASSESerr[ind], ALL_MDOTSerr[ind], ALL_UPP[ind]
            idx = np.isfinite(X2) & np.isfinite(Y2)
            Xfin, Yfin, Xfinerr, Yfinerr, uppfin = X2[idx], Y2[idx],X2err[idx], Y2err[idx], upp2[idx]

            y_main1 = linfit.at['CASPAR all', 'slope']*Xfin + linfit.at['CASPAR all', 'intercept']
            y_BF = m*Xfin + b

            r2 = r2_score(Yfin, y_main1)
            MAE = (1/len(y_main1)) * np.sum(abs(Yfin-y_main1)) # lower = better
            MSE = (1/len(y_main1)) * np.sum((Yfin-y_main1)**2.) # lower = better
            MSAE = np.sqrt((1/len(y_main1)) * np.sum((Yfin-y_main1)**2.)) # lower= better
            rss = np.sum((Yfin-y_main1)**2.)
            aic = len(y_main1) * np.log10(rss/len(y_main1)) + 2*2 # more negative = better
            chi2 = (1/(len(y_main1)-2)) * np.nansum((Yfin-y_main1)**2. / Yfinerr)
            print('stars CASPAR aic:', aic)

            r2 = r2_score(Yfin, y_BF)
            MAE = (1/len(y_BF)) * np.sum(abs(Yfin-y_BF)) # lower = better
            MSE = (1/len(y_BF)) * np.sum((Yfin-y_BF)**2.) # lower = better
            MSAE = np.sqrt((1/len(y_BF)) * np.sum((Yfin-y_BF)**2.)) # lower= better
            rss = np.sum((Yfin-y_BF)**2.)
            aic = len(y_BF) * np.log10(rss/len(y_BF)) + 2*2 # more negative = better
            chi2 = (1/(len(y_BF)-2)) * np.nansum((Yfin-y_BF)**2. / Yfinerr)
            print('stars BF aic:', aic)

            ##########
            ind = np.where((ALL_MASSES <= -1.1249))
            X2, Y2, X2err, Y2err, upp2 = ALL_MASSES[ind], ALL_MDOTS[ind],ALL_MASSESerr[ind], ALL_MDOTSerr[ind], ALL_UPP[ind]
            idx = np.isfinite(X2) & np.isfinite(Y2)
            Xfin, Yfin, Xfinerr, Yfinerr, uppfin = X2[idx], Y2[idx],X2err[idx], Y2err[idx], upp2[idx]

            y_main1 = linfit.at['CASPAR all', 'slope']*Xfin + linfit.at['CASPAR all', 'intercept']
            y_BF = m*Xfin + b

            r2 = r2_score(Yfin, y_main1)
            MAE = (1/len(y_main1)) * np.sum(abs(Yfin-y_main1)) # lower = better
            MSE = (1/len(y_main1)) * np.sum((Yfin-y_main1)**2.) # lower = better
            MSAE = np.sqrt((1/len(y_main1)) * np.sum((Yfin-y_main1)**2.)) # lower= better
            rss = np.sum((Yfin-y_main1)**2.)
            aic = len(y_main1) * np.log10(rss/len(y_main1)) + 2*2 # more negative = better
            chi2 = (1/(len(y_main1)-2)) * np.nansum((Yfin-y_main1)**2. / Yfinerr)
            print('bds CASPAR aic:', aic)

            r2 = r2_score(Yfin, y_BF)
            MAE = (1/len(y_BF)) * np.sum(abs(Yfin-y_BF)) # lower = better
            MSE = (1/len(y_BF)) * np.sum((Yfin-y_BF)**2.) # lower = better
            MSAE = np.sqrt((1/len(y_BF)) * np.sum((Yfin-y_BF)**2.)) # lower= better
            rss = np.sum((Yfin-y_BF)**2.)
            aic = len(y_BF) * np.log10(rss/len(y_BF)) + 2*2 # more negative = better
            chi2 = (1/(len(y_BF)-2)) * np.nansum((Yfin-y_BF)**2. / Yfinerr)
            print('bds BF aic:', aic)


    
    
    for axi in ax:
        axi.set_xlim(-2.5, 0.3)  
        plt.setp(axi.get_xticklabels(), visible=False)
        axi.set_xlabel('')
        if axi in [ax1, ax5]:
            axi.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)', fontsize=fs)
        else:
            axi.set_ylabel('')
            plt.setp(axi.get_yticklabels(), visible=False)
            
    for axi in axtop:
        axi.set_xlim(-2.5, 0.3)  
        axi.set_xticklabels('')
        axi.set_ylim(2,300)
        axi.set_yticks([10,100])
        axi.set_yscale('log')
        cplot.plot_tickparams(axi, fs=fs)
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
        

    for axi in axres:
        axi.set_xlim(-2.5, 0.3)  
        axi.set_ylim(-3,3)
        axi.axvline(x = -1.124938736, color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
        axi.axvline(x = -1.906371, color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
        cplot.plot_tickparams(axi, fs=fs)
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
        
    for axi in axsideres:
        axi.set_xlim(1,800)  
        axi.set_ylim(-3,3)
        cplot.plot_tickparams(axi, fs=fs)
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
        plt.setp(axi.get_yticklabels(), visible=False)
        axi.set_xscale('log')
        
    for axi in axside:
        axi.set_xlim(1,800) 
        axi.set_ylim(-13.95, -5)
        cplot.plot_tickparams(axi, fs=fs)
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
        plt.setp(axi.get_yticklabels(), visible=False)
        plt.setp(axi.get_xticklabels(), visible=False)
        axi.set_xscale('log')
        
    for axi in [ax6, ax8]:
        axi.set_xlabel('log $M\ (M_{\odot})$', fontsize=fs)
    for axi in [ax2, ax6]:
        axi.set_ylabel('residual\n(dex)', fontsize=fs)
    for axi in [ax4, ax8, ax3top, ax7top]:
        plt.setp(axi.get_yticklabels(), visible=False)
    for axi in [ax2, ax4, ax2side, ax4side]:
        axi.set_xticklabels([])
        plt.setp(axi.get_xticklabels(), visible=False)
    for axi in [ax1top, ax5top]:
        axi.set_ylabel('N', fontsize=fs)
    for axi in [ax6side, ax8side]:
        axi.set_xlabel('N', fontsize=fs)

    plt.subplots_adjust()
    plt.tight_layout()
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    plt.show()
    
    
def Ha10_UVexcess_combo(df_caspar, **kwargs):
    fig, (ax2, ax1) = plt.subplots(figsize=(5,8), ncols=1, nrows=2)
    
    
    # HA 10 vs HA
    IND = (~df_caspar['Ha 10% Accretion Rate'].isna()) & (~df_caspar['Ha Accretion Rate'].isna())

    Ha10 = df_caspar['Ha 10% Accretion Rate'][IND].values
    Ha = df_caspar['Ha Accretion Rate'][IND].values
    Haerr = df_caspar['Ha Accretion Rate err'][IND].values
    d = df_caspar['Distance'][IND].values
    M = df_caspar['Mass'][IND].values
    age = df_caspar['Age'][IND].values
    COM = df_caspar['Companion'][IND].values
    yerr = 0.434 * Ha/Haerr
   
    for i in np.arange(len(Ha10)):
        
        if M[i] > 0.075:
            col='dimgray'
            marker='o'
        elif (M[i] <= 0.075) & (COM[i] == 'COM'):
            col='magenta'
            marker='s'
        else:
            col='cyan'
            marker='D'
            
        ax1.scatter(np.log10(M[i]), np.log10(Ha10[i]/Ha[i]), c=col, marker=marker, ec='k', s=50,zorder=100)
        
    cplot.plot_MMdotFeatures(None, ax1, fs=16, fs_BL=10, fs_ticks=16,
                             burning_lines=True, tick_params=True,
                            xlabel='log Mass (M$_\odot$)', ylabel=r'log[$\frac{ \dot{M}\ (\mathrm{H}\alpha \ 10\%)}{\dot{M} (H\alpha)}$] (dex)',
                            xlim=(-2.5, 0.3), ylim=(-2.5, 2.5))
    
    ### UV excess 
    IND = (~df_caspar['Ha Accretion Rate'].isna()) & (df_caspar['Accretion Diagnostic']=='Continuum Excess')

    Ha = df_caspar['Ha Accretion Rate'][IND].values
    UV = df_caspar['Accretion Rate'][IND].values
    UVerr = df_caspar['Accretion Rate err'][IND].values
    UVlogerr = 0.434 * UV/UVerr
    
    d = df_caspar['Distance'][IND].values
    M = df_caspar['Mass'][IND].values
    age = df_caspar['Age'][IND].values
    COM = df_caspar['Companion'][IND].values
    
    for i in np.arange(len(Ha10)):
        
        if M[i] > 0.075:
            col='dimgray'
            marker='o'
        elif (M[i] <= 0.075) & (COM[i] == 'COM'):
            col='magenta'
            marker='s'
        else:
            col='cyan'
            marker='D'
        im=ax2.scatter(np.log10(M[i]), np.log10(UV[i]/Ha[i]), c=col, marker=marker, ec='k', s=50,zorder=100)

    cplot.plot_MMdotFeatures(None, ax2, fs=16, fs_BL=10, fs_ticks=16,
                         burning_lines=True, tick_params=True,
                        xlabel='log Mass (M$_\odot$)', 
                             ylabel=r'log[$\frac{\dot{M}\ (\mathrm{Continuum\ Excess})}{\dot{M}\ (\mathrm{H}\alpha)}$] (dex)',
                        xlim=(-2.5, 0.3), ylim=(-2,2))

    ax1.text(-2.4, -1.8, 'Hα 10% width', fontsize= 16)
    ax2.text(-2.4, -1.8, 'continuum excess', fontsize=16)
    

    for ax in [ax1, ax2]:
        ax.axhline(0, color='k', linewidth=3)
#        cbar = plt.colorbar(im, ax=ax)
#        cbar.set_label('log Mass (M$_\odot$)',fontsize=16)
#        cbar.ax.tick_params(labelsize=16, direction='in', length=7)

#        ax.tick_params(labelsize=16)
#        ax.tick_params(which='major', length=7, direction='in', top=True, right=True, bottom=True)
#        ax.tick_params(which='minor', length=4, direction='in', top=True, right=True, bottom=True)
#        ax.minorticks_on()
        
    markerlegend(ax2, 'k', 'cyan', 'magenta', loc='upper right')
    
    plt.tight_layout()
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    plt.show()
    

