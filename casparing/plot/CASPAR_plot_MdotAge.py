import CASPAR_util as cutil
import CASPAR_fitMMdot as cfit
import CASPAR_plotting as cplot
import CASPAR_sortdata as csort

from scipy.stats import norm
from scipy.optimize import curve_fit
from scipy import interpolate
import matplotlib.colors as COLORS
import pandas as pd
import numpy as np
import matplotlib.cm as cm
from scipy import odr as odr
import statsmodels.api as sm

import matplotlib.pyplot as plt
import os

HERE = os.path.dirname(os.path.abspath(__file__))
fil = os.path.join(HERE, '../data/linearfits_data/betti2023_capar_linearfits.csv')
linfit = pd.read_csv(fil, delimiter=',', comment='#', skiprows=[0])

linfit = linfit.set_index('obj')
linfit = linfit.replace(-1000, np.nan)

def _AccDiaglegend(ax, loc, fs):

    color = cutil.AccDiagColor()

    fin_col=[]
    lab, co = [], []
    for col in np.sort(list(color.keys())):
        if 'Ha Photometric Line Strength' in col:
            pass
        else:
            if 'Ha' in col:
                label = 'Hα' + col.split('Ha')[-1]
            else:
                label=col
            s = ax.scatter([],[], color=color[col],s=0)
            fin_col.append(s)
            lab.append(label)
            co.append(color[col])
    first_legend = ax.legend(fin_col,lab,  loc=loc, ncol=1, fontsize=fs/2,title='Accretion Diagnostic',
               title_fontsize=fs/2, markerscale=0.5, labelcolor=co, handletextpad=-1.55, frameon=True) 
    first_legend._legend_box.align = "left"
    ax.add_artist(first_legend)
        
    
def _run_scaleByMass(objs, REG, UPP, fitname, doUPP=True, db='CASPAR'):
    REG['log Mdot orig'] = REG['log Mdot'] 
    fitparams = cfit.linmix_fitting(objs, fitname=db+' ' +fitname,
                                   plot=False, as_list=False)

    N, slope, slopeerrlow, slopeerrupp, intercept, intercepterrlow, intercepterrupp, stdfit, r2 = fitparams       
    C = REG['log Mdot'] - slope*REG['log Mass']
    NEWMdot = slope*np.log10(0.7) + C

    REG['log Mdot'] = NEWMdot


    if doUPP:
        UPP['log Mdot orig'] = UPP['log Mdot']
        C = UPP['log Mdot'] - slope*UPP['log Mass']
        NEWMdot = slope*np.log10(0.7) + C
        UPP['log Mdot'] = NEWMdot

        return REG, UPP
    else:
        return REG, None





def MdotVSage(df_caspar, colorby='mass', scaleByMass=True, fit_data=True, **kwargs):
    
    cmap=kwargs.get('cmap', 'turbo')
    s = kwargs.get('s', 20)
    fs = kwargs.get('fs', 12)
    mainlegend= kwargs.get('mainlegend'), 
    markerlegend=kwargs.get('markerlegend')
    mainloc=kwargs.get('mainloc', 'lower left')
    markerlegendloc = kwargs.get('markerlegendloc','lower left')
    dbname = kwargs.get('dbname', 'CASPAR')
    fitname=kwargs.get('fitname','doLinearFit')
    xlim = kwargs.get('xlim', [5.5, 7.9])
    
    fig, ax = cplot.plot_createSinglePlot()

    stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_caspar)
    print('std log age = ', round(np.nanstd(df_caspar['log Age']),3), 'yr ' )
    
    if colorby == 'accretion diagnostic':
        cS = Sreg['color']
        cB = BDreg['color']
        cP = Preg['color'] 
        cmap=None
        _AccDiaglegend(ax, mainloc, fs)
    elif colorby == 'mass':
        cS = Sreg['log Mass']
        cB = BDreg['log Mass']
        cP = Preg['log Mass']
    else:
        cS, cB, cP ='k', 'cyan', 'magenta'
        
    if scaleByMass:
        Sreg, Supp = _run_scaleByMass(stars.copy(), Sreg.copy(), Supp.copy(), 'all', db=dbname)
        BDreg, BDupp = _run_scaleByMass(bds.copy(), BDreg.copy(), BDupp.copy(), 'all', db=dbname)
        Preg, Pupp = _run_scaleByMass(planets.copy(), Preg.copy(), Pupp.copy(), 'all', db=dbname)
        df_caspar, _ = _run_scaleByMass(df_caspar.copy(), df_caspar.copy(), df_caspar.copy(), 'all', doUPP=False, db=dbname)
        cmap=cmap
        
    im = _scatterplot(ax, Sreg['log Age'], Sreg['log Mdot'],Supp['log Age'], Supp['log Mdot'], c=cS, s=s,  alpha=0.65, marker='o', edgecolor='none',cmap=cmap)
    
    if colorby == 'mass':
        cbaxes = fig.add_axes([0.91, 0.125, 0.03, 0.755])
        cbar = fig.colorbar(im, ax=ax, cax=cbaxes)
        cbar.set_label('log Mass ($M_\odot$)', fontsize=fs)
        cbar.ax.tick_params(which='both',direction='in', labelsize=fs)
        cbar.ax.minorticks_on()
    
    im = _scatterplot(ax, BDreg['log Age'], BDreg['log Mdot'], BDupp['log Age'], BDupp['log Mdot'], cB, s, cmap, 0.65, 'D', 'none')    
    
    im = _scatterplot(ax, Preg['log Age'], Preg['log Mdot'], Pupp['log Age'], Pupp['log Mdot'], cP, s+3, cmap, 0.8, 's', 'k', lw=1)   
    
    if fit_data:
        if fitname == 'lit_lines':
            x = np.linspace(3.2, 8, 50)
            y = -1.32 - (1.07 * x)
            line1, = ax.plot(x,y,color='k',ls='-', zorder=4, linewidth=1)
            ax.fill_between(x, y-0.5, y+0.5, color='gray', alpha=0.4, zorder=-1)
        else:
            if fitname != 'doLinearFit':
                fitname = fitname + ' MdotVSAge'
            print('max log Mdot = ',round(df_caspar['log Mdot'].max(), 3), 'Msun/yr' )
            fitparams = cfit.linmix_fitting(df_caspar, fitname=fitname, xval='Age', 
                                       plot=False, ax=ax, extend=False, as_list=False, color='k')
            N, slope, slopeerrlow, slopeerrupp, intercept, intercepterrlow, intercepterrupp, stdfit, r2 = fitparams
            print('best fit parameters')
            print(f'N = {N}, std= {stdfit:.3}, R2 = {r2:.3}')
            print(f'slope = {slope:.3}_{slopeerrlow:.3}^{slopeerrupp:.3}, intercept = {intercept:.3}_{intercepterrlow:.3}^{intercepterrupp:.3}')

    #### LEGEND
    if markerlegend:
        if colorby != 'object':
            cS, cB, cP = 'k', 'k', 'k'
        cplot.plot_objtype_legend(ax, cS, cB, cP, fs=fs,loc=markerlegendloc)
    
    ax.set_xlabel('log Age (yr)', fontsize=fs)
    ax.set_ylabel(r'log $\dot{M}\ (M_\odot)$', fontsize = fs)
    
    cplot.plot_tickparams(ax, fs=fs, minor='on')
    ax.set_xlim(xlim[0], xlim[1])

    if 'savefig' in kwargs:
        if colorby == 'Mass':
            plt.savefig(kwargs['savefig'], dpi=150,transparent=False, 
                   bbox_extra_artists=(cbaxes,), bbox_inches='tight')
        else:
            plt.savefig(kwargs['savefig'], dpi=150,transparent=False)
    plt.show()
    


def _scatterplot(ax, X, Y, Xupp, Yupp, c, s, cmap, alpha, marker, edgecolor, lw=0.5):
    
    im = ax.scatter(X, Y, c=c, s=s,
            cmap=cmap, alpha=alpha, vmin=-2.5, vmax=0.6, 
            marker=marker, edgecolor=edgecolor, linewidth=lw)
    
    u = np.zeros_like(Xupp)
    v = -np.ones_like(Xupp)
    if cmap is not None:
        norm = COLORS.Normalize(vmin=-2.5, vmax=0.6)
        colmap = cm.get_cmap(cmap)
        ax.quiver(Xupp, Yupp, u,v,  headwidth=4, headaxislength=3,headlength=3, width=.005, scale=30, alpha=alpha, color=colmap(norm(c)))
    else:
        ax.quiver(Xupp, Yupp, u,v,  headwidth=4, headaxislength=3,headlength=3, width=.005, scale=30, alpha=alpha, color=c)
    
    return im
    


def _age_grouping(df_caspar, age, i, group=True):
    if (i == 0) and group==True:
        
        df_ageupdateS = csort.CASPAR_separateByAge(df_caspar, age[i], age[i+1])
        stars, _, planets, [Sreg, Supp], _, [Preg, Pupp] = csort.CASPAR_separateByMass(df_ageupdateS)     
        df_ageupdate = csort.CASPAR_separateByAge(df_caspar, age[i], age[i+2])
        _, bds, _, _, [BDreg, BDupp],_ = csort.CASPAR_separateByMass(df_ageupdate)

    elif (i == 1) and group==True:
        
        df_ageupdateS = csort.CASPAR_separateByAge(df_caspar, age[i], age[i+1])
        stars, _, planets, [Sreg, Supp], _, [Preg, Pupp] = csort.CASPAR_separateByMass(df_ageupdateS)
        
        df_ageupdate = csort.CASPAR_separateByAge(df_caspar, age[i-1], age[i+1])
        _, bds, _, _, [BDreg, BDupp],_ = csort.CASPAR_separateByMass(df_ageupdate)

    else:
      
        df_ageupdate = csort.CASPAR_separateByAge(df_caspar, age[i], age[i+1])
        stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_ageupdate)                                                  
    return df_ageupdate, stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp]
    
    

def _fit_mass_regions(df, object_type):
    if len(df) > 0:
        object_type = object_type.lower()
        if (object_type == 'bd') or (object_type == 'bds'):
            ot_fit = 'bd'
            ot_linfit = 'BDs'
            

        if (object_type == 'star') or (object_type == 'stars'):
            ot_fit = 'stars'
            ot_linfit = 'stars'
            
            
        if object_type == 'all':
            ot_fit = 'all'
            ot_linfit = 'all'
            
        intercept, intercept_err = cfit.run_odrI(df['log Mass'], df['log Mdot'],
             df['log Mass err'], df['log Mdot err'], ot_fit)

        slope = linfit.at[f'CASPAR {ot_linfit}', 'slope']  

        DBL = slope*(-1.906)+intercept[0]
        DBLerr = np.sqrt(intercept_err[0]**2. + linfit.at[f'CASPAR BDs', 'slope err low']**2.)

        HBL = slope*(-1.12)+intercept[0]
        HBLerr = np.sqrt(intercept_err[0]**2. + linfit.at[f'CASPAR BDs', 'slope err low']**2.)
        return DBL, DBLerr, HBL, HBLerr, slope, intercept[0]
    else:
        return [np.nan], [np.nan], [np.nan], [np.nan], None,None

                        
def MMdotVSage(df_caspar, fit_mass_regions=True, group=True):
    age, _, _, _ = cplot.plot_age()
    AGE = []
    AGEerr = []
    HBL = []
    HBLerr = []
    DBL = []
    DBLerr = []
    for i in np.arange(len(age)-1):
        df_ageupdate, stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = _age_grouping(df_caspar, age, i, group)
        
        AGE.append(np.mean([age[i], age[i+1]]))
        AGEerr.append(df_ageupdate['Age err'].median())
        
        if fit_mass_regions:
            bd_DBL, bd_DBLerr, bd_HBL, bd_HBLerr, bdslope, bdintercept = _fit_mass_regions(bds, 'bds')
            s_DBL, s_DBLerr, s_HBL, s_HBLerr, sslope, sintercept = _fit_mass_regions(stars, 'stars')
            DBL_fin = [s_DBL, bd_DBL]
            DBLerr_fin = [s_DBLerr, bd_DBLerr]
            HBL_fin = [s_HBL, bd_HBL]
            HBLerr_fin = [s_HBLerr, bd_HBLerr]
            
        else:

            df_caspar2 = df_caspar.loc[df_caspar['Companion']!= 'COM']
            DBL_fin, DBLerr_fin, HBL_fin, HBLerr_fin, allslope, allintercept = _fit_mass_regions(df_caspar2, 'all')
        HBL.append(HBL_fin)
        HBLerr.append(HBLerr_fin)
        DBL.append(DBL_fin)
        DBLerr.append(DBLerr_fin)

    return AGE, AGEerr, [HBL, DBL], [HBLerr, DBLerr]

def plot_MMdotVSage(df_caspar, group=True, fit_mass_regions=True, plot_residuals=True, plot_pdf=True, **kwargs):
    fs = kwargs.get('fs', 12)
    if plot_residuals:
        if plot_pdf:
            fig, ((ax, No_ax), (axN, axP)) = plt.subplots(dpi=150,nrows=2, ncols=2,figsize=(5,4), gridspec_kw={'width_ratios':[3,1], 'height_ratios':[1.5,1]})
            No_ax.axis('off')
            axP.axhline(0, color='k', linewidth=1, zorder=-1)
        else:
            fig, (ax, axN) = plt.subplots(figsize=(5,4), dpi=300, nrows=2, ncols=1, sharey=False, gridspec_kw={'height_ratios':[3,1], 'wspace':0, 'hspace':0})
        axN.axhline(0, color='k', linewidth=1, zorder=-1)

    else:
        fig, ax = cplot.plot_createSinglePlot()
    cplot.plot_graylines(ax, df_caspar)
    
    bdx = np.linspace(-1.9, -1.12, 50) # BD
    sx = np.linspace(-1.12, 5, 50) #star
    allx = np.linspace(-2.5, 5, 50)
    
    age, marker, c, darkc = cplot.plot_age()
    combo_age_residuals = []
    for i in np.arange(len(age)-1):    

        cc='peru' if (i <= 1) and (group==True) else c[i]
        dc = 'saddlebrown' if (i <= 1) and (group==True) else darkc[i]
        mm='p' if (i <= 1) and (group==True) else marker[i]
        s = 60 if mm == '*' else 25
        
        if age[i]==1:
            zorder=-10
            label='(' + str(age[i]) + '-' + str(age[i+1]) + ']'
        elif age[i] == 8:
            zorder=-age[i]
            label = '> 8'
        elif age[i+1] == 1:
            zorder=-age[i]
            label = '$\leq$ 1' 
        else:
            zorder=-age[i]
            label='(' + str(age[i]) + '-' + str(age[i+1]) + ']'
            
        age_residuals = []
            
        df_ageupdate, stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = _age_grouping(df_caspar, age, i, group=group)
            
        ax.scatter([],[], marker=mm, color=cc, label=label, s=20, zorder=i)
        if len(stars)>0:
            cplot.plot_MMdot(ax, Sreg['log Mass'], Sreg['log Mdot'], Supp['log Mass'], Supp['log Mdot'], s=s, color=c[i], marker=mm, edgecolor='none', alpha=0.5, zorder=zorder)
        if len(bds)>0:
            
            cplot.plot_MMdot(ax, BDreg['log Mass'], BDreg['log Mdot'], BDupp['log Mass'], BDupp['log Mdot'], s=s, color=cc, marker=mm, edgecolor='none', alpha=0.5, zorder=zorder)
        if len(planets)>0:
            cplot.plot_MMdot(ax, Preg['log Mass'], Preg['log Mdot'], Pupp['log Mass'], Pupp['log Mdot'], s=s, color=c[i], marker=mm, edgecolor='k', facecolor=c[i], alpha=0.6, zorder=zorder, linewidth=0.5)
        
        
        if fit_mass_regions:
            if len(bds)>0:
                _, _, _, _, bdslope, bdintercept = _fit_mass_regions(bds, 'bds')
                ax.plot(bdx, bdslope*bdx + bdintercept, '-', color=dc, linewidth=2)
            if len(stars)>0:
                _, _, _, _, sslope, sintercept = _fit_mass_regions(stars, 'stars')
                ax.plot(sx, sslope*sx + sintercept, '-', color=darkc[i], linewidth=2)
                ax.plot(bdx, sslope*bdx + sintercept, '--', color=darkc[i], linewidth=1)
            
            if group:
                if i == 0:
                    label2='$\leq$' + str(age[i+2]) 
                    ax.scatter([],[], marker=mm, color=cc, s=20, label=label2)
                
            if plot_residuals:
                if len(stars)>0:
                    f = interpolate.interp1d(sx, sslope*sx + sintercept)              
                    cplot.plot_MMdot(axN, Sreg['log Mass'], Sreg['log Mdot']-f(Sreg['log Mass']), Supp['log Mass'], Supp['log Mdot']-f(Supp['log Mass']), s=s, color=c[i], marker=mm, edgecolor='none', alpha=0.6, zorder=zorder)
                    age_residuals.extend((stars['log Mdot']-f(stars['log Mass'])).values)
                if len(bds)>0:
                    f = interpolate.interp1d(bdx, bdslope*bdx + bdintercept) 
                    cplot.plot_MMdot(axN, BDreg['log Mass'], BDreg['log Mdot']-f(BDreg['log Mass']), BDupp['log Mass'], BDupp['log Mdot']-f(BDupp['log Mass']), s=s, color=cc, marker=mm, edgecolor='none', alpha=0.6, zorder=zorder)
                    age_residuals.extend((bds['log Mdot']-f(bds['log Mass'])).values)
                                
            
        else:
 
            df_caspar2 = df_caspar.loc[df_caspar['Companion']!= 'COM']
            _, _, _, _, slope, intercept = _fit_mass_regions(df_caspar2, 'all')
            
            ax.plot(allx, slope*allx + intercept, '-', color=darkc[i], linewidth=2)
               
            if plot_residuals:
                f = interpolate.interp1d(allx,slope*allx+intercept)
                df_casparreg, df_casparupp = csort.CASPAR_separateByUppLimit(df_caspar2)
                cplot.plot_MMdot(axN, df_casparreg['log Mass'], df_casparreg['log Mdot']-f(df_casparreg['log Mass']), df_casparupp['log Mass'], df_casparupp['log Mdot']-f(df_casparupp['log Mass']), s=s, color=cc, marker=mm, edgecolor='none', alpha=0.6, zorder=zorder)
                age_residuals.extend((df_caspar2['log Mdot']-f(df_caspar2['log Mass'])).values)
                
        if plot_residuals:
            combo_age_residuals.extend(age_residuals) 
            if plot_pdf:
                axP.hist(age_residuals, color=cc, orientation='horizontal', histtype='step', bins=12, density=True, linewidth=1.5)

    combo_age_residuals = np.sort(np.array(combo_age_residuals))
    if plot_residuals and plot_pdf:
        _, bins, _ = axP.hist(combo_age_residuals, color='k', histtype='step', orientation='horizontal', bins=12, density=True, label=r'$M$-$\dot{M}$ age fit')
        if fit_mass_regions:
            stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_caspar)
            
            sy = linfit.at['CASPAR stars', 'slope']*sx + linfit.at['CASPAR stars', 'intercept']
            fstar = interpolate.interp1d(sx,sy)

            allstar_residuals = stars['log Mdot']-fstar(stars['log Mass'])
            
            bdy = linfit.at['CASPAR BDs', 'slope']*bdx + linfit.at['CASPAR BDs', 'intercept']
            fbd = interpolate.interp1d(bdx,bdy)
            allbd_residuals = bds['log Mdot']-fbd(bds['log Mass'])
            combo_all_residuals = np.append(allbd_residuals, allstar_residuals)
            all_pdf_label = r'$M$-$\dot{M}$ star/BD fit'
            
        else:
  
            ally = linfit.at['CASPAR all', 'slope']*allx + linfit.at['CASPAR all', 'intercept']
            fall = interpolate.interp1d(allx,ally)
            combo_all_residuals = df_caspar['log Mdot']-fall(df_caspar['log Mass'])
            all_pdf_label = r'$M$-$\dot{M}$ fit'
            
        axP.hist(combo_all_residuals, color='gray', orientation='horizontal', histtype='stepfilled', bins=12, density=True, zorder=0, alpha=0.6,label=all_pdf_label)
        axP.legend(fontsize=fs/2.2, loc='lower right',)
        
        cplot.plot_MMdotFeatures(fig, axP, fs=fs,ylabel=None, fs_ticks=fs,  burning_lines=True,burning_lines_labels=False, tick_params=True,xlabel='PDF', xlim=(0, 0.8), ylim=(-3.5, 3.5))
        axP.set_yticklabels([])
            
            
    cplot.plot_MMdotFeatures(fig, ax, fs=fs,  fs_ticks=fs, burning_lines=True, burning_lines_labels=True, tick_params=True, 
                            xlim=(-2.5, 0.3), ylim=(-13, -5))
    ax.set_xlabel('')
    ax.set_xticklabels([])
    
    if plot_residuals:
        cplot.plot_MMdotFeatures(fig, axN, fs=fs, fs_ticks=fs, 
                                 xlabel='log Mass ($M_{\odot}$)', 
                                 ylabel='log $\dot{M}/\dot{M}_\mathrm{age\ fit}$',  burning_lines=True,
                                 burning_lines_labels=False, 
                                 tick_params=True,
                                 xlim=(-2.5, 0.3),ylim=(-3.5, 3.5))

    ax.legend(ncol=1, fontsize=fs/1.5, title='Age (Myr)', title_fontsize=fs/1.5, loc='upper left')
    if group:
        ax.get_legend().remove()
        ph = [plt.plot([],marker="", ls="")[0]]*3
        l2 = ax.scatter([],[], color=c[0], s=s, marker=marker[0])
        l3 = ax.scatter([],[] , color=c[1], s=s, marker=marker[1])
        l4 = ax.scatter([],[] ,color=c[2], s=s, marker=marker[2])
        l5 = ax.scatter([],[], color=c[3],s=s,  marker=marker[3])
        l7 = ax.scatter([],[], color='peru', s=s,  marker='p')
        l8 = ax.scatter([],[], color=c[2], s=s, marker=marker[2])
        l9 = ax.scatter([],[], color=c[3], s=s, marker=marker[3])
        labels = ['Star Ages (Myr)', '$\leq$ 1','(1-3]','(3-8])', '> 8', '',
                 'BD Ages (Myr)', '$\leq$ 3','(3-8]', '> 8' ]
        leg = ax.legend([ph[0], l2, l3, l4, l5, ph[1], ph[2], l7, l8, l9], 
                   labels,
                   loc='upper left', ncol=1, fontsize=fs/2, framealpha=1)
        for vpack in leg._legend_handle_box.get_children():
            for i, hpack in enumerate(vpack.get_children()):
                if i == 0:

                    hpack.get_children()[0].set_width(0)
                elif i==6:
                    hpack.get_children()[0].set_width(0)
        
    
    plt.subplots_adjust(hspace=0, wspace=0)


    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)

        
    if plot_residuals:
        if plot_pdf:
            return fig, (ax, axN, axP)
        else:
            return fig, (ax, axN)
    else:
        return fig, ax
    
    


    


    
def _interceptVSAge(ax, x,xerr,y, yerr, fit_mass_regions=True, **kwargs):
    if fit_mass_regions:
        y = np.array(y)
        yerr = np.array(yerr)

        ybd = y[1::,1]
        y = y[:,0]
        
        yerrbd = yerr[1::,1]
        yerr = yerr[:,0]
        
        xbd = x[1::]
        xerrbd = xerr[1::]
        xbd[-1]=16

    x[-1]=16
        
    Smfc = kwargs.get('star_mfc', 'k')
    Bmfc = kwargs.get('bd_mfc', 'darkcyan')
    
    ls = kwargs.get('ls', '-')
    
    ax.errorbar(x, y, xerr = xerr, yerr= yerr, fmt='ko', mfc=Smfc, mec='k')
    
    if fit_mass_regions:
        ax.errorbar(xbd, ybd,xerr = xerrbd, yerr=yerrbd,fmt= 'D', color='darkcyan', mfc=Bmfc, mec='darkcyan')

    def func(B, x):
        return B[0] * np.exp(-B[1] * x) + B[2]

    data = odr.RealData(x, y, sx = xerr, sy = yerr)
    quad_model = odr.Model(func)
    Odr = odr.ODR(data, quad_model,beta0 =[1.1, 0.15, -8.6])
    Odr.set_job(fit_type=0)
    out =Odr.run()
    vals = out.beta
    #print('stars:', vals)
    xx = np.linspace(0,20, 100)
    label=kwargs.get('label', 'intercept')
    ax.plot(xx, func(vals, xx), 'k',ls=ls, linewidth=1, label=label)

    if fit_mass_regions: 
        data = odr.RealData(xbd, ybd, sx = xerrbd, sy = yerrbd)
        Odr = odr.ODR(data, quad_model,beta0 =[1.1, 0.15, -7.02])
        Odr.set_job(fit_type=0)
        out =Odr.run()
        vals = out.beta
        #print('BD:', vals)
        xx = np.linspace(0,20, 100)
        ax.plot(xx, func(vals, xx), linestyle=ls,color='darkcyan', linewidth=1)
        ax.plot([],[], color='k', linestyle=ls,  linewidth=1)
        
def interceptVSAge(df_caspar, fit_mass_regions=True, **kwargs):
    x, xerr, y, yerr = MMdotVSage(df_caspar, fit_mass_regions=fit_mass_regions, group=True) 

    fig, ax = cplot.plot_createSinglePlot()
    
    if fit_mass_regions:
        s1=ax.scatter([],[], label='Brown Dwarf', color='darkcyan', marker='D', )
        s2=ax.scatter([],[], label='Star', color='k', marker='o', )
    
    _interceptVSAge(ax, x, xerr,y[0], yerr[0], fit_mass_regions=fit_mass_regions, ls='-', label='HBL $M = 0.075\ M_\odot$', **kwargs) #HBL
    _interceptVSAge(ax, x, xerr,y[1], yerr[1], fit_mass_regions=fit_mass_regions, ls='--', label='DBL $M = 0.012\ M_\odot$', Smfc='none', Bmfc='none',**kwargs) #DBL
      
    fs=kwargs.get('fs', 10)
    
    cplot.plot_tickparams(ax, fs=fs, minor='on')
    
    ax.set_xlabel('Age (Myr)', fontsize = fs)
    ax.set_ylabel('$\log \dot{M}$ at Burning Limit ($M_\odot$ yr$^{-1}$)', fontsize = fs)
    
    ax.set_ylim(kwargs.get('ylim',(-13.4, -8)))
    ax.legend(fontsize=8, ncol=2, loc='upper right')

    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
    
    return fig, ax
 

def BDMdotvsAgeInteractionModel(df_caspar, **kwargs):
    fs = kwargs.get('fs', 8)
    fig, ax = cplot.plot_createSinglePlot()

    marker='D'

    stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_caspar)
    df_caspar = df_caspar.loc[df_caspar['Companion']!='COM']
    
    reg, upp = csort.CASPAR_separateByUppLimit(df_caspar)


    def mult_lin_fit(data, reg, upp, datatype, extend=True, inleg=True):
        if datatype == 'all':
            mass = data['log Mass'].values
            age = data['Age'].values
            Mdot = data['log Mdot'].values
        else:
            mass = reg['log Mass'].values
            age = reg['Age'].values
            Mdot = reg['log Mdot'].values

        d = {'mass':mass, 'age':age, 'Mdot':Mdot}
        df = pd.DataFrame(d)

        model = sm.OLS.from_formula('Mdot ~ mass + age + mass:age',df)
        result = model.fit()
        stdfit = np.sqrt(result.scale)
        intercept = result.params['Intercept']
        Smass = result.params['mass']
        Sage = result.params['age']
        combo = result.params['mass:age']
        print(result.params)
        def newmod(mass, mdot, age, c, extend=True):
            newmodel = intercept + Smass*mass + Sage*age + combo*mass*age
            RES = newmodel - mdot
            stds = np.nanstd(RES)
            
            mass_big = np.linspace(mass.min(), mass.max(), 100)
            newmodel = intercept + Smass*mass_big + Sage*age + combo*mass_big*age
            
            print(age, (Smass+combo*age), 'x  ', intercept+Sage*age)
#            print(stds)
            ax.plot(mass_big, newmodel, color=c, zorder=20)
            ax.fill_between(mass_big, newmodel-stds, newmodel+stds, color=c, alpha=0.1, edgecolor='none', zorder=4)
            mass_ext = np.linspace(-1.95, -1.09, 20)
            newmodel = intercept + Smass*mass_ext + Sage*age + combo*mass_ext*age
            if extend:
                ax.plot(mass_ext, newmodel, color=c, linestyle='--', zorder=19)
            


        age_list = [1, 2, 4, 10, 16]
        
        zorder=[10,9,8,7,6]
        c = cm.magma_r(np.linspace(0.3,0.9,len(age_list)))
        for i in np.arange(len(age_list)):
            if inleg:
                ax.scatter([],[], color=c[i], s=0, label=f'{age_list[i]} Myr')
            newmod(mass, Mdot, age_list[i], c=c[i], extend=extend)

            if datatype == 'all':
                massreg = reg['log Mass'].values
                agereg = reg['Age'].values
                Mdotreg = reg['log Mdot'].values
                massupp = upp['log Mass'].values
                ageupp = upp['Age'].values
                Mdotupp = upp['log Mdot'].values
                if age_list[i] == 4:
                    offset=1.5
                elif age_list[i] == 2:
                    offset = 0.015
                elif age_list[i]==10:
                    offset=2
                else:
                    offset=0.5

                INDreg = np.where((agereg >= age_list[i]-offset)& (agereg < age_list[i]+offset))
                INDupp = np.where((ageupp >= age_list[i]-offset)& (ageupp < age_list[i]+offset))
                ax.scatter(massreg[INDreg], Mdotreg[INDreg], marker=marker, color=c[i], s=15, zorder=zorder[i])
                u = np.zeros_like(massupp[INDupp])
                v = -np.ones_like(massupp[INDupp])
                ax.quiver(massupp[INDupp], Mdotupp[INDupp], u, v, headwidth=4, headaxislength=3, 
                          headlength=3, width=0.005, scale=30, alpha=1, color=c[i], zorder=zorder[i])

            else:
                IND = np.where(age == age_list[i])
                ax.scatter(mass[IND], Mdot[IND], marker=marker, color=c[i], s=15, zorder=zorder[i])
                
        return c

    c = mult_lin_fit(stars, Sreg, Supp, 'all', extend=False)
    c = mult_lin_fit(bds, BDreg, BDupp, 'all', extend=True, inleg=False)
    cplot.plot_MMdotFeatures(fig,ax, fs=fs)

    t = ax.text(-1.8, -6.1, r'$\log \dot{M}_\mathrm{star} = 2.22 \log M_\mathrm{star} - 0.14 t -0.11 (\log M_\mathrm{star}\times t) -7.50$'+'\n' + r'$\log \dot{M}_\mathrm{BD} = 1.28 \log M_\mathrm{BD} + 0.32 t + 0.31 (\log M_\mathrm{BD}\times t) -8.29$', fontsize=6,zorder=100)
    t.set_bbox(dict(facecolor='white', alpha=0.5, edgecolor='none'))
    ax.legend(ncol=1,loc = 'upper left', fontsize=fs, handlelength=0, handletextpad=0, labelcolor=c)

    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)




def BDMdotvsAgeInteractionModelRES(df_caspar, **kwargs):
    fs = kwargs.get('fs', 6)
    fig, (ax, ax2) = plt.subplots(figsize=(5,4), dpi=300, nrows=2, ncols=1, sharey=False, gridspec_kw={'height_ratios':[3,1], 'wspace':0, 'hspace':0})
    
    marker='D'

    stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(df_caspar)
    df_caspar = df_caspar.loc[df_caspar['Companion']!='COM']
    
    reg, upp = csort.CASPAR_separateByUppLimit(df_caspar)


    def mult_lin_fit(data, reg, upp, datatype, extend=True, inleg=True):
        if datatype == 'all':
            mass = data['log Mass'].values
            age = data['Age'].values
            Mdot = data['log Mdot'].values
        else:
            mass = reg['log Mass'].values
            age = reg['Age'].values
            Mdot = reg['log Mdot'].values

        d = {'mass':mass, 'age':age, 'Mdot':Mdot}
        df = pd.DataFrame(d)

        model = sm.OLS.from_formula('Mdot ~ mass + age + mass:age',df)
        result = model.fit()
        stdfit = np.sqrt(result.scale)
        intercept = result.params['Intercept']
        Smass = result.params['mass']
        Sage = result.params['age']
        combo = result.params['mass:age']
        def newmod(mass, mdot, age, c, extend=True):
            newmodel = intercept + Smass*mass + Sage*age + combo*mass*age
            RES = newmodel - mdot
            stds = np.nanstd(RES)
            
            mass_big = np.linspace(mass.min(), mass.max(), 100)
            newmodel = intercept + Smass*mass_big + Sage*age + combo*mass_big*age
            print(age, (Smass+combo*age), 'x  ', intercept+Sage*age)
            ax.plot(mass_big, newmodel, color=c, zorder=20)
            ax.fill_between(mass_big, newmodel-stds, newmodel+stds, color=c, alpha=0.1, edgecolor='none', zorder=4)
            mass_ext = np.linspace(-1.95, -1.09, 20)
            newmodel_ext = intercept + Smass*mass_ext + Sage*age + combo*mass_ext*age
            if extend:
                ax.plot(mass_ext, newmodel_ext, color=c, linestyle='--', zorder=19)
                
            return mass_big, newmodel
            

            
        RESes = []    
        age_list = [1, 2, 4, 10,16]
        offset = [0.5, 0.015, 1.5, 0.5, 2]
        zorder=[10,9,8,7,6]
        c = cm.magma_r(np.linspace(0.3,0.9,len(age_list)))
        for i in np.arange(len(age_list)):
            if inleg:
                ax.scatter([],[], color=c[i], s=0, label=f'{age_list[i]} Myr')
            MODmass, MODmdot = newmod(mass, Mdot, age_list[i], c=c[i], extend=extend)
 
            if datatype == 'all':
                massreg = reg['log Mass'].values
                agereg = reg['Age'].values
                Mdotreg = reg['log Mdot'].values
                massupp = upp['log Mass'].values
                ageupp = upp['Age'].values
                Mdotupp = upp['log Mdot'].values
                

                INDreg = np.where((agereg >= age_list[i]-offset[i])& (agereg < age_list[i]+offset[i]))
                INDupp = np.where((ageupp >= age_list[i]-offset[i])& (ageupp < age_list[i]+offset[i]))
                ax.scatter(massreg[INDreg], Mdotreg[INDreg], marker=marker, color=c[i], s=15, zorder=zorder[i])
                u = np.zeros_like(massupp[INDupp])
                v = -np.ones_like(massupp[INDupp])
                ax.quiver(massupp[INDupp], Mdotupp[INDupp], u, v, headwidth=4, headaxislength=3, 
                          headlength=3, width=0.005, scale=30, alpha=1, color=c[i], zorder=zorder[i])
                
                
                allIND = np.where((agereg >= age_list[i]-offset[i])& (agereg < age_list[i]+offset[i]))
                allmass = massreg[allIND]
                allmdot = Mdotreg[allIND]
                f = interpolate.interp1d(MODmass, MODmdot)
                newMODmdot = f(allmass)
                allRES = allmdot - newMODmdot
                ax2.scatter(allmass, allRES, color=c[i], s=15)
                RESes.extend(list(allRES))
            

            else:
                IND = np.where((age >= age_list[i]-offset)& (age < age_list[i]+offset))
                ax.scatter(mass[IND], Mdot[IND], marker=marker, color=c[i], s=15, zorder=zorder[i])

        print(np.nanstd(RESes))
        return c, RESes
    print('stars')
    c, RESesS = mult_lin_fit(stars, Sreg, Supp, 'all', extend=False)
    print('BDs')
    c, RESesBD = mult_lin_fit(bds, BDreg, BDupp, 'all', extend=True, inleg=False)
    print('all')
    allRES = RESesS +RESesBD
    print(np.nanstd(allRES))
    cplot.plot_MMdotFeatures(fig,ax, fs=fs, xlabel=None)
    cplot.plot_MMdotFeatures(fig,ax2, fs=fs, ylabel='residuals', burning_lines=False, burning_lines_labels=False, ylim=(-3,3))
    
    ax2.axhline(0, color='k', zorder=-1)
    ax.set_xticklabels([])
    ax.legend(ncol=1,loc = 'upper left', fontsize=8, handlelength=0, handletextpad=0, labelcolor=c)

    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
#    plt.show()