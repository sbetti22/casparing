import seaborn as sns
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import .CASPAR_plotting as cplot
import .CASPAR_util as cutil
import .CASPAR_sortdata as csort
import .CASPAR_fitMMdot as cfit

def _plot_refVSnew(ax, df_x, df_y, colname, AD, label, legend=False, **kwargs):
        axlims= kwargs.get('axlims', None)
        fs = kwargs.get('fs', 10)
        s = kwargs.get('ms', 30)
        color = cutil.AccDiagColor()
        div = 1e3 if colname == 'Teff' else 1
        XX = df_x[colname]/div
        YY = df_y[colname]/div
        sns.scatterplot(XX, YY,
                        hue=df_y['AD'], s=s, ax=ax, legend=legend, 
                        zorder=10, palette=color, hue_order=AD, edgecolor='k',alpha=1)         
        if legend:
            ax.get_legend().remove()   

        if kwargs.get('verbose'):
            per_dif = abs(YY-XX)/XX
            IND = np.where(per_dif < .5)
            print(f'percent of values within a factor of 2 between lit/caspar for {label}:',len(IND[0])/len(XX))

        if axlims is None:
            axlims = [min(np.min(XX), np.min(YY)), max(np.max(XX), np.max(YY))]

        ax.set_xlim(axlims)
        ax.set_ylim(axlims)
        x = np.arange(axlims[0], axlims[1]*2)
        ax.plot(x,x,'k')

        ax.set_xlabel(f'literature {label}', fontsize=fs)
        ax.set_ylabel(f'CASPAR {label}', fontsize=fs)  
        
def _axis_label(key):
    d =  {'Age': 'age (Myr)', 
           'Distance':'distance (pc)',
           'Mass': 'mass ($M_\odot$)',
           'Teff':'effective\ntemperature ($10^3$ K)' ,
           'Radius': 'radius ($R_\odot$)'}
    return d[key]
        
def refVSnew(df_lit, df_caspar, **kwargs):
    fs = kwargs.get('fs', 16)
    s = kwargs.get('ms', 30)
    
    df_lit = df_lit.sort_values(by=['Unique ID'], ignore_index=True)
    df_caspar = df_caspar.sort_values(by=['Unique ID'], ignore_index=True)
    
    AD = sorted(df_caspar.AD.unique())

    fig, ((ax0, ax1, ax3), (ax2, ax4, ax5)) = plt.subplots(figsize=(12,8),dpi=150, nrows=2, ncols=3)
    ax5.axis('off')
    
    _plot_refVSnew(ax0, df_lit, df_caspar, 'Age',AD, _axis_label('Age'), fs=fs, s=s, axlims=[0.08, 60])
    _plot_refVSnew(ax1, df_lit, df_caspar, 'Distance', AD, _axis_label('Distance'), fs=fs, s=s, axlims=[20, 8000])
    _plot_refVSnew(ax2, df_lit, df_caspar, 'Mass', AD, _axis_label('Mass'), fs=fs, s=s, axlms=[3e-3, 4])
    _plot_refVSnew(ax3, df_lit, df_caspar, 'Teff', AD, _axis_label('Teff'), fs=fs, s=s, axlims=[1.7, 5.5])
    _plot_refVSnew(ax4, df_lit, df_caspar, 'Radius', AD,  _axis_label('Radius'), fs=fs, s=s, axlims=[6e-2, 10], legend=True)  
        
    for i, ax in enumerate([ax0, ax1, ax2,ax3, ax4]):
        if i != 3:
            ax.set_xscale('log')
            ax.set_yscale('log')
        ax.tick_params(which='major', direction='in', top=True, right=True, length=5,labelsize=fs-2)
        ax.tick_params(which='minor', direction='in', top=True, right=True, length=3,labelsize=fs-2)
        ax.minorticks_on()
        ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    handles, labels = ax4.get_legend_handles_labels()
    leg=fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.7,0.2),
               fontsize=fs-3, ncol=1, title='Accretion Diagnostic', facecolor='none', title_fontsize=15)

    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    if 'savefig' in kwargs:
         plt.savefig(kwargs['savefig'], dpi=300, transparent=True, 
                    bbox_extra_artists=(leg,), bbox_inches='tight')
    plt.show()
    
    
def refVSnew_single(df_lit, df_caspar, phyParam, **kwargs):
    xyscale = kwargs.get('xyscale', 'linear')
    axlims = kwargs.get('axlims', None)
    fs = kwargs.get('fs', 10)
    leg_fs = fs-5
    if leg_fs  <= 0 :
        leg_fs = 1
    
    
    df_lit = df_lit.sort_values(by=['Unique ID'], ignore_index=True)
    df_caspar = df_caspar.sort_values(by=['Unique ID'], ignore_index=True)
    
    AD = sorted(df_caspar.AD.unique())

    if 'ax' in kwargs:
        fig = kwargs['fig']
        ax = kwargs['ax']
    else:
        fig, ax = cplot.plot_createSinglePlot()
    
    _plot_refVSnew(ax, df_lit, df_caspar, phyParam, AD, _axis_label(phyParam), legend=True, **kwargs)  
        
    if xyscale == 'log':
        ax.set_xscale('log')
        ax.set_yscale('log')
    ax.tick_params(which='major', direction='in', top=True, right=True, length=5,labelsize=fs-2)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=3,labelsize=fs-2)
    ax.minorticks_on()
    
    handles, labels = ax.get_legend_handles_labels()
    leg=ax.legend(handles, labels, loc='lower right',
               fontsize=leg_fs, ncol=1, title='Accretion Diagnostic', facecolor='none', title_fontsize=leg_fs)

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    if 'savefig' in kwargs:
         plt.savefig(kwargs['savefig'], dpi=300)
    plt.show()
    
    
def refVSnew_mdotresidual(df_lit, df_caspar, **kwargs):
    fs = kwargs.get('fs', 12)
    s = kwargs.get('ms', 20)
    
    df_lit = df_lit.sort_values(by=['Unique ID'], ignore_index=True)
    df_caspar = df_caspar.sort_values(by=['Unique ID'], ignore_index=True)

    colname = 'log Mdot'
    AD = sorted(df_caspar.AD.unique())

    fig, ((ax1, ax11),(ax2,ax22), (ax3,ax33))= plt.subplots(nrows=3, ncols=2, 
                                 figsize=(5,6), dpi=150, sharey=True, gridspec_kw={'width_ratios':[4,1]})
    axes = [ax1, ax2, ax3]
    axhist = [ax11, ax22, ax33]
    for i, accdiag in enumerate(AD[0:-1]):
        df_caspar_accdiag = df_caspar[colname].loc[df_caspar['AD'] == accdiag]
        df_lit_accdiag = df_lit[colname].loc[df_caspar['AD'] == accdiag]
        df_residual = df_caspar_accdiag-df_lit_accdiag
        hue = df_caspar['AD'].loc[df_caspar['AD'] == accdiag]
        
        axes[i].scatter(df_caspar_accdiag, df_residual, edgecolor='k', s=s, facecolor='none')

        sns.scatterplot(x=df_caspar_accdiag, y=df_residual,
                        hue=hue, s=s, ax=axes[i],
                        zorder=10, palette=cutil.AccDiagColor(), hue_order=AD, edgecolor='k',alpha=0.4, legend=False)  

        axes[i].axhline(0, color='k')
        
        axes[i].set_ylabel(r'log $\frac{\dot M_{CASPAR}}{\dot M_{lit}}$ (dex)', fontsize=fs)  
        
        t=axes[i].text(0.05, 0.8, AD[i], fontsize=fs+2, transform=axes[i].transAxes, zorder=1000)
        t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
        
        ################################
    
        # best fit of data
        FWHM = 2.355 * np.std(df_residual)
        (mu, sigma) = norm.fit(df_residual)

        # the histogram of the data
        binns = int(np.sqrt(len(df_residual)))
        n, bins, patches = axhist[i].hist(df_residual, binns, orientation="horizontal",density=True, 
                      color=cutil.AccDiagColor()[accdiag], histtype='stepfilled', edgecolor='k', align='mid')

        axhist[i].axhline(mu-2.355*sigma/2, color='gray', linewidth=0.75, )
        axhist[i].axhline(mu+2.355*sigma/2, color='gray',linewidth=0.75, )
        if  i==2:
            axhist[i].text(1,mu+2.355*sigma/2+0.08, 'FWHM', color='gray', fontsize=6)
        ################################
        axes[i].set_xlim(-14, -5.5)
        axes[i].set_ylim(-2.5, 2.5)
        axes[i].tick_params(which='major', direction='in', top=True, right=True, length=5,labelsize=fs)
        axes[i].tick_params(which='minor', direction='in', top=True, right=True, length=3,labelsize=fs)
        axes[i].minorticks_on()
        
        axhist[i].axhline(0, color='k')
        axhist[i].set_xlim(0, 1.55)
        axhist[i].set_ylim(-2.5,2.5)
        axhist[i].set_xticks([0.2, 0.8, 1.4])
        axhist[i].tick_params(which='major', direction='in', top=True, right=True, length=5,labelsize=fs)
        axhist[i].tick_params(which='minor', direction='in', top=True, right=True, length=3,labelsize=fs)
        axhist[i].minorticks_on()
        
        if i == 2:
            axhist[i].set_xlabel('PDF', fontsize=fs)
            axes[i].set_xlabel(r'CASPAR log $\dot{M}\ (M_\odot$/yr)', fontsize=fs)
        else:
            axes[i].set_xlabel('') 
            axes[i].set_xticklabels('')
            axhist[i].set_xticklabels('')
             
    XX = 10**df_lit_accdiag
    YY = 10**df_caspar_accdiag
    # % difference
    per_dif = abs(YY-XX)/XX
    IND = np.where(per_dif <= .5) # the values that changed by 50% or less
    print(len(IND[0]), len(XX))
    print(f'percent of values that changed by 50% or less between lit/caspar for Mdot:',len(IND[0])/len(XX))
                                      
                                                                         
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.1, wspace=0)
    
    if 'savefig' in kwargs:
         plt.savefig(kwargs['savefig'], dpi=300, transparent=False)
    plt.show()

    
