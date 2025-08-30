import CASPAR_plotting as cplot
import CASPAR_util as cutil
import CASPAR_sortdata as csort
import CASPAR_fitMMdot as cfit

import os
import numpy as np
import pandas as pd
import matplotlib.colors as mc
import matplotlib.pyplot as plt
from scipy import interpolate

def _colorbySFR_cbar(fig, ax, plot, vmin, vmax,cmap, fs):
    cbaxes = fig.add_axes([0.98, 0.19, 0.03, 0.76])
    cbar = fig.colorbar(plot, ax=ax, cax=cbaxes)
    clu_vals = {'σ Ori':[400, 'd'],'25 Ori':[354, 'o'],'118 Tau':[100, '1'], 'Argus':[80, 'v'],
                'β Pic':[40, 'H'], 'Cha I':[190, '^'],'IC 348':[321, '2'], 'Lagoon Nebula':[1326, 'p'], 
                'NGC 2024':[420, '*'],'η Cha':[94, 'h'], '\n    Tuc-Hor, 50 pc\n    TW Hya':[50, ['D', '4']], 'Lupus':[158, 'P'], 'CrA': [149, 's'],  'US':[146,'<'], 
                ' ρ Oph & Taurus':[140, ['3', 'x']],'UCL':[130, 'X'], }
    clu_cbar = pd.DataFrame.from_dict(clu_vals, orient='index', columns=['dist', 'marker']).reset_index(names='sfr')

    pos = cbar.ax.get_position()
    ax2 = cbar.ax.twinx()
    ax2.set_ylim([vmin, vmax])
    ax2.set_yscale('log')
    pos.x0 +=0.05
    cbar.ax.set_position(pos)
    ax2.set_position(pos)
    cbar.ax.tick_params(axis='y',which='both', labelsize=0, left=False)
    ax2.tick_params(axis='y',which='both', labelsize=fs/2, direction='in', )

    vals = clu_cbar['dist'].to_list()
    names = clu_cbar['sfr'].to_list()
    markers = clu_cbar['marker'].to_list()
    c = clu_cbar['dist'].astype(str).to_list()

    c[0] = '406'
    lab = ['     ' + s1 + ', ' + s2 + ' pc' for s1, s2 in zip(names, c)]
    lab[-4] = '     CrA, 149 pc; US, 146 pc'
    lab[-3] = ''
    ax2.minorticks_on()
    ax2.yaxis.set_ticks([])
    ax2.yaxis.set_ticks(vals, minor=True)
    ax2.yaxis.set_ticklabels(lab, minor=True)
    ax2.tick_params(axis='y',which='minor', labelsize=3.) 
    ax2.tick_params(axis='y',which='major', labelsize=fs/2.2, direction='in', length = 2, pad=-8, ) 
    for i in np.arange(len(names)):
        if isinstance(markers[i], str):
            ax2.scatter(1.7, vals[i], marker=markers[i], norm=mc.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, c=vals[i], s=4, alpha=0.7, clip_on=False, edgecolor='k', linewidth=0.5)
        else:
            ax2.scatter(1.7, vals[i], marker=markers[i][0], norm=mc.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, c=vals[i], s=4, alpha=0.7, clip_on=False, edgecolor='k', linewidth=0.5)
            if 'Taurus' in names[i]:
                xx, yy_off = 2., 0
            elif 'TW Hya' in names[i]:
                xx, yy_off = 1.7, 3.5
            else:
                xx, yy_off = 1.7, 7
            ax2.scatter(xx, vals[i]-yy_off, marker=markers[i][1], norm=mc.LogNorm(vmin=vmin, vmax=vmax), cmap=cmap, c=vals[i], s=4, alpha=0.7, clip_on=False, edgecolor='k', linewidth=0.5)
    return cbaxes


def _colorbySFR(fig, ax, df_caspar, fs):
    cluster, marker, distance = cplot.plot_cluster(dist=True)
    cmap = plt.cm.jet
    cluster_color = plt.cm.jet(np.linspace(0.2, 1, len(cluster)-1))
    vmin, vmax= 35, 1350

    df_caspar['Main Association'] = df_caspar['Main Association'].fillna('none')
    bd_P = df_caspar.loc[(df_caspar['Companion']=='COM')&(df_caspar['Mass']<=0.075)]
    for i, name in enumerate(cluster):
        bd_cluster = df_caspar[df_caspar['Main Association'].str.contains(name)]
        bd_clusterP = bd_P[bd_P['Main Association'].str.contains(name)]

        REG, UPP = csort.CASPAR_separateByUppLimit(bd_cluster)
        REGP, UPPP = csort.CASPAR_separateByUppLimit(bd_clusterP)

        p = cplot.plot_MMdot(ax, REG['log Mass'], REG['log Mdot'], UPP['log Mass'], UPP['log Mdot'], color=REG['Distance'], uppcolor=UPP['Distance'], color_cmap=True, marker=marker[i], cmap=cmap, edgecolor=None, alpha=0.7, vmin=vmin, vmax=vmax)
        
        cplot.plot_MMdot(ax, REGP['log Mass'], REGP['log Mdot'], UPPP['log Mass'], UPPP['log Mdot'], color=REGP['Distance'], uppcolor=UPPP['Distance'], color_cmap=True, marker=marker[i], cmap=cmap, edgecolor='k', alpha=0.7, vmin=vmin, vmax=vmax)
        
    ec=None
    cbaxes= _colorbySFR_cbar(fig, ax, p, vmin, vmax,cmap, fs)
    return cbaxes, ec

                
def _colorbyAccDiag(main_ax, Sreg, BDreg, Preg, Supp, BDupp, Pupp, fs, legend=True, **kwargs):

    if 'color' not in list(Sreg.columns):
        Sreg['color'] = csort.AccDiag_AssD(Sreg)[3]
        BDreg['color'] = csort.AccDiag_AssD(BDreg)[3]
        Preg['color'] = csort.AccDiag_AssD(Preg)[3]
        Supp['color'] = csort.AccDiag_AssD(Supp)[3]
        BDupp['color'] = csort.AccDiag_AssD(BDupp)[3]
        Pupp['color'] = csort.AccDiag_AssD(Pupp)[3]
        
    cS, cBD, cP = Sreg['color'], BDreg['color'], Preg['color']
    cUPPS, cUPPBD, cUPPP = Supp['color'], BDupp['color'], Pupp['color']
    cSdark, cBDdark, cPdark = 'k', 'k', 'k'
    alpha=0.7

    if legend:
        color = cutil.AccDiagColor()

        fin_col, lab, co = [], [], []
        for col in np.sort(list(color.keys())):
            s = main_ax.scatter([],[], color=color[col],s=0)
            fin_col.append(s)
            lab.append(col)
            co.append(color[col])

        first_legend = main_ax.legend(fin_col,lab, loc='lower right', ncol=1, fontsize=fs/2,title='Accretion Diagnostic',
                       title_fontsize=fs/2, markerscale=0.5, labelcolor=co, handletextpad=-1.55, frameon=True) 
        first_legend._legend_box.align = "left"
        first_legend.set_zorder(1000)
        main_ax.add_artist(first_legend)
    edgecolor='none'
    return cS, cBD, cP, cUPPS, cUPPBD, cUPPP, cSdark, cBDdark, cPdark, alpha, edgecolor

def MdotVSmass(df_caspar, colorby='object', fit_data=True, **kwargs):
    mainlegend = kwargs.get('mainlegend')
    mainlegendloc = kwargs.get('mainlegendloc', 'lower right')
    s = kwargs.get('s', 12)
    fs = kwargs.get('fs', 12)
    
    if 'ax' in kwargs:
        fig = kwargs['fig']
        ax = kwargs['ax']
    else:
        fig, ax = cplot.plot_createSinglePlot()
    
    stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_caspar)

    cplot.plot_graylines(ax, df_caspar)
    
    if colorby == 'sfr':
        cbaxes, ec = _colorbySFR(fig, ax, df_caspar, fs)
        cS, cBD, cP, cUPPS, cUPPBD, cUPPP, cSdark, cBDdark, cPdark, alpha = 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 1
    
    else:
        if colorby=='object':
            cS, cBD, cP, cUPPS, cUPPBD, cUPPP, cSdark, cBDdark, cPdark, alpha,ec = cutil.ObjectColor()
        elif colorby == 'accretion diagnostic':
            cS, cBD, cP, cUPPS, cUPPBD, cUPPP, cSdark, cBDdark, cPdark, alpha,ec = _colorbyAccDiag(ax, Sreg, BDreg, Preg, Supp, BDupp, Pupp,fs,  legend=mainlegend,**kwargs)
        else:
            cS, cBD, cP, cUPPS, cUPPBD, cUPPP, cSdark, cBDdark, cPdark, alpha,ec = 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 1, 'k'

        if 'color' in kwargs:
            C = kwargs['color']
            cS, cBD, cP,cUPPS, cUPPBD, cUPPP = C, C, C, C, C, C
        
            
        cplot.plot_MMdot(ax, Sreg['log Mass'], Sreg['log Mdot'], Supp['log Mass'], Supp['log Mdot'], color_cmap=False, s=s, color=cS, marker='o', edgecolor=ec, linewidth=0.5,alpha=alpha, uppcolor=cUPPS, zorder=100)

        cplot.plot_MMdot(ax, BDreg['log Mass'], BDreg['log Mdot'], BDupp['log Mass'], BDupp['log Mdot'], color_cmap=False, s=s, color=cBD, marker='D', edgecolor=ec, linewidth=0.5,alpha=alpha, uppcolor=cUPPBD, zorder=100)

        cplot.plot_MMdot(ax, Preg['log Mass'], Preg['log Mdot'], Pupp['log Mass'], Pupp['log Mdot'], color_cmap=False, s=s+2, color=cP, marker='s', edgecolor='k', alpha=alpha, linewidth=0.75, uppcolor=cUPPP, zorder=100)

    cplot.plot_MMdotFeatures(fig, ax, fs=fs)
    
    # Line fitting here
    _line_fitting(ax, df_caspar, fit_data=fit_data, colorby=colorby, **kwargs)
    
    ############### legend ################
    if colorby == 'accretion diagnostic':
        cS, cBD, cP = 'k', 'k', 'k'
        if 'mainlegendloc' not in kwargs:
            mainlegendloc = 'upper left'
    cplot.plot_objtype_legend(ax, cS, cBD, cP, fs=fs, loc=mainlegendloc)
    #######################################

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    if 'savefig' in kwargs:
        if colorby=='SFR':
            plt.savefig(kwargs['savefig'], dpi=300, transparent=False, bbox_extra_artists=(cbaxes,), bbox_inches='tight')
        else:
            plt.savefig(kwargs['savefig'], dpi=300, transparent=False)
    plt.show()
    
def _line_fitting(ax, df_caspar, fit_data=True, colorby='object', **kwargs):
    ############## LINE FITTING ######################
    if colorby=='object':
        _, _, _, _, _, _, cSdark, cBDdark, cPdark, _,_ = cutil.ObjectColor()
    elif 'color' in kwargs:
        cSdark, cBDdark, cPdark = kwargs['color'], kwargs['color'], kwargs['color']
    else:
        cSdark, cBDdark, cPdark = 'k', 'k', 'k'
             
    if fit_data:
        fitname = kwargs.get('fitname', 'doLinearFit')
        print(fitname)
        if kwargs.get('fit_mass_regions'):
            if fitname == 'doLinearFit':
                eS, eBD, eP = '','',''
            else:
                eS, eBD, eP = ' stars', ' BDs', ' planet'
            print(fitname + eS)
            N, slope, slopeerrlow,slopeerrupp, intercept, intercepterrlow,intercepterrupp, stdfit, r2 = cfit.linmix_fitting(stars, fitname=fitname + eS, xval='logMass', plot=True, as_list=False, ax=ax, color=cSdark, extend=True)
            print('stars fitting information')
            print('#', N, 'slope: ', slope, slopeerrlow, slopeerrupp, 'intercept: ', intercept, intercepterrlow, intercepterrupp, 'stdfit: ', stdfit, 'R2: ', r2)
            
            N, slope, slopeerrlow,slopeerrupp, intercept, intercepterrlow,intercepterrupp, stdfit, r2 = cfit.linmix_fitting(bds, fitname=fitname + eBD, xval='logMass', plot=True, as_list=False, ax=ax, color=cBDdark, extend=True)
            print('BDs fitting information')
            print('#', N, 'slope: ', slope, slopeerrlow, slopeerrupp, 'intercept: ', intercept, intercepterrlow, intercepterrupp, 'stdfit: ', stdfit, 'R2: ', r2)
            
            N, slope, slopeerrlow,slopeerrupp, intercept, intercepterrlow,intercepterrupp, stdfit, r2 = cfit.linmix_fitting(planets, fitname=fitname + eP, xval='logMass', plot=True, as_list=False, ax=ax, color=cPdark, extend=True)
            print('PMCs fitting information')
            print('#', N, 'slope: ', slope, slopeerrlow, slopeerrupp, 'intercept: ', intercept, intercepterrlow, intercepterrupp, 'stdfit: ', stdfit, 'R2: ', r2)
            
        else:
            if fitname != 'doLinearFit':
                fitname  = fitname + ' all'
            N, slope, slopeerrlow,slopeerrupp, intercept, intercepterrlow,intercepterrupp, stdfit, r2 = cfit.linmix_fitting(df_caspar, fitname=fitname, xval='logMass', plot=True, as_list=False, ax=ax, color=cSdark, extend=False)
            print('stars fitting information')
            print('#', N, 'slope: ', slope, slopeerrlow, slopeerrupp, 'intercept: ', intercept, intercepterrlow, intercepterrupp, 'stdfit: ', stdfit, 'r: ', r2)
        
    else:
        cplot.plot_litlines(ax)


def _plot_hist(ax, data, color, alpha, bins=None, bottom=0):
    if bins is None:
        nbins = int(np.sqrt(len(data)))
    else:
        nbins=bins
    _,B,_=ax.hist(data, color=color, histtype='stepfilled', edgecolor='k', linewidth=1.5,alpha=alpha, bottom=bottom,
           bins=nbins, orientation='horizontal', density=True, zorder=10)
    ax.set_xlim(0,2.8)
    ax.axvline(bottom, color='k')
    ax.plot([bottom, bottom+1], [np.mean(data), np.mean(data)], color='dimgray',
           linewidth=0.75, linestyle='-', zorder=9)
    return B

def MdotVSmass_residuals(df_caspar,**kwargs):
    fs = kwargs.get('fs', 12)
    fig, ((ax, ax11),(ax2,ax22), (ax3,ax33), (ax4,ax44), (ax5, ax55)) = plt.subplots(figsize=(6,6), nrows=5, ncols=2, dpi=150,  gridspec_kw={'height_ratios':[3,1,1,1,1], 'width_ratios':[3,1]},)
    
    cplot.plot_graylines(ax, df_caspar)

    x = np.linspace(-5.5,5,10)
    yZhou=2.43*x-7.9 #Zhou2014
    yHartmann=2.1*x-7.9 #Hartmann2015
    
    _, slopestars, _, _, interceptstars, _, _, stdfitstars, _= cfit.linmix_fitting(df_caspar, fitname='CASPAR stars')
    yCASPARstars = slopestars*x + interceptstars#CASPAR stars
    
    _, slope, _, _, intercept, _, _, stdfit, _= cfit.linmix_fitting(df_caspar, fitname='CASPAR all')
    yCASPAR = slope*x + intercept#CASPAR ALL

    fZhou = interpolate.interp1d(x, yZhou)
    fHartmann = interpolate.interp1d(x, yHartmann)
    fCASPAR = interpolate.interp1d(x, yCASPAR)
    fCASPARstars = interpolate.interp1d(x, yCASPARstars)
    axi = [ax2, ax3, ax4, ax5]
    axh = [ax22, ax33, ax44, ax55]
    fits = [fZhou, fHartmann, fCASPARstars, fCASPAR]
    fitname = ['Zhou+2014', 'Hartmann+2016', 'CASPAR stars','CASPAR']
    combocolor = ['k', 'cyan', 'magenta']
    facecolor = ['k', 'cyan', 'magenta']
    edgecolor = ['gray', 'k', 'k']
    ls=['-.', ':', '--', '-']
    color=['brown', 'gray', 'k', 'k']
    spreadsheet = ['CASPAR']
    lw = 0.75
    stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_caspar)
    ### M-Mdot plot
    cplot.plot_MMdot(ax, Sreg['log Mass'], Sreg['log Mdot'], Supp['log Mass'], Supp['log Mdot'], color_cmap=False, s=15, color=facecolor[0], marker='o', linewidth=lw, edgecolor=edgecolor[0], alpha=0.4, zorder=10)

    cplot.plot_MMdot(ax, BDreg['log Mass'], BDreg['log Mdot'], BDupp['log Mass'], BDupp['log Mdot'], color_cmap=False, s=15, color=facecolor[1], marker='D',linewidth=lw, edgecolor=edgecolor[1], alpha=0.4, zorder=10)

    cplot.plot_MMdot(ax, Preg['log Mass'], Preg['log Mdot'], Pupp['log Mass'], Pupp['log Mdot'], color_cmap=False, s=15, color=facecolor[2], marker='s', edgecolor=edgecolor[2], alpha=0.4, linewidth=lw, zorder=10)
      
    B = _plot_hist(ax11, stars['log Mdot'].values, 'k', 1)
    _plot_hist(ax11, bds['log Mdot'].values, 'cyan', 1,bins=B, bottom=1)
    _plot_hist(ax11, planets['log Mdot'].values, 'magenta', 1,bins=B, bottom=2)
    ax11.set_ylim(-13.95, -5)
    ax11.set_yticklabels([])
    ax11.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
    ax11.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
    ax11.minorticks_on()
    for val, f in enumerate(fits):
        #### RESIDUALS
        cplot.plot_MMdot(axi[val], Sreg['log Mass'], Sreg['log Mdot']-f(Sreg['log Mass']), Supp['log Mass'], Supp['log Mdot']-f(Supp['log Mass']), color_cmap=False, s=15, color=facecolor[0], marker='o', linewidth=0.5, edgecolor=edgecolor[0], alpha=0.4, zorder=10)

        cplot.plot_MMdot(axi[val], BDreg['log Mass'], BDreg['log Mdot']-f(BDreg['log Mass']), BDupp['log Mass'], BDupp['log Mdot']-f(BDupp['log Mass']), color_cmap=False, s=15, color=facecolor[1], marker='D',linewidth=0.5, edgecolor=edgecolor[1], alpha=0.4, zorder=10)

        cplot.plot_MMdot(axi[val], Preg['log Mass'], Preg['log Mdot']-f(Preg['log Mass']), Pupp['log Mass'], Pupp['log Mdot']-f(Pupp['log Mass']), color_cmap=False, s=15, color=facecolor[2], marker='s', edgecolor=edgecolor[2], alpha=0.4, linewidth=0.5, zorder=10)
        
        Bh = _plot_hist(axh[val], (bds['log Mdot']-f(bds['log Mass'])).values, 'cyan',1, bottom=1)
        _plot_hist(axh[val], (stars['log Mdot']-f(stars['log Mass'])).values, 'k',1,bins=Bh )
        
        _plot_hist(axh[val], (planets['log Mdot']-f(planets['log Mass'])).values, 'magenta', 1,bins=Bh, bottom=2)
        
        for c, obj in enumerate([stars, bds, planets]):
            YFIT = f(obj['log Mass'])
            VAL = obj['log Mdot']-YFIT
        
        axi[val].text(0.02, 0.1, fitname[val], transform=axi[val].transAxes)
        axi[val].set_ylabel('\nresidual\n(dex)', fontsize=fs)
        axi[val].set_ylim(-4.75,4.75)
        axh[val].set_ylim(-4.75, 4.75)
        axi[val].set_yticks([-3,0,3])
        axh[val].set_yticks([-3,0,3])
        axh[val].set_yticklabels([])
        axi[val].axhline(0, color=color[val], linestyle=ls[val], linewidth=2, zorder=10)
        axh[val].axhline(0, color=color[val], linestyle=ls[val], linewidth=2, zorder=10)
        axi[val].axvline(x = np.log10(0.075), color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
        axi[val].axvline(x = np.log10(13 * .0009543), color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
        axi[val].tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
        axi[val].tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
        axi[val].minorticks_on()
        axh[val].tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
        axh[val].tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
        axh[val].minorticks_on()

    cplot.plot_MMdotFeatures(fig, ax, fs=fs)

    lw = 2
    ax.plot(x, yZhou, linestyle=ls[0], color=color[0], label='Zhou+2014', linewidth=lw, zorder=10) 
    
    ax.plot(x, yCASPARstars, linestyle=ls[2], color=color[2], label='CASPAR stars', linewidth=lw, zorder=10)
    
    
    ax.plot(x, yCASPAR, linestyle=ls[3], color=color[3], label='CASPAR', linewidth=lw+1, zorder=10) 
    ax.plot(x, yHartmann, linestyle=ls[1], color=color[1], label='Hartmann+2016', linewidth=lw, zorder=10) 
    
    
    
    ax.fill_between(x, yCASPAR-stdfit, 
                    yCASPAR+stdfit, color=color[3], alpha=0.1, zorder=-1)
    ax.fill_between(x, yCASPARstars-stdfitstars, 
                    yCASPARstars+stdfitstars, color=color[2], alpha=0.1, zorder=-1)
    ax.fill_between(x, yHartmann-0.75, yHartmann+0.75, color=color[1], alpha=0.1, zorder=-2)
    ax.fill_between(x, yZhou-0.5, yZhou+0.5, color=color[0], alpha=0.1, zorder=-3)

    ############### legend ################
    cplot.plot_objtype_legend(ax, 'k', 'cyan', 'magenta', fs=fs/1.25, loc='lower right')
    ##########################################################################
    ax5.set_xlabel('log Mass ($M_{\odot}$)', fontsize = fs)
    ax55.set_xlabel('Normalized N\n+offset', fontsize=fs)
    ax.set_xticklabels([])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0, wspace=0)
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=150, transparent=False) 
    plt.show()

