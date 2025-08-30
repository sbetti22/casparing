import .CASPAR_plotting as cplot
import .CASPAR_fitMMdot as cfitmmdot
import .CASPAR_sortdata as csort
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

def MdotVSmassVSsfr(df_caspar, lit_lines=True, **kwargs):
    fs = kwargs.get('fs', 8)
    fig, axes = plt.subplots(figsize=(6,8), nrows=5, ncols=2, dpi=150, sharex=True, sharey=True)
    axi = axes.flatten()

    cluster = {'Lagoon Nebula':[1326, 'p'], 'Chamaeleon I':[190, '^'],'Taurus':[135, 'x'],'σ Ori':[406, 'd'],'p Oph':[150, '3'],'Lupus':[158, 'P'],  'TW Hya':[51, '4'],'Upper Scorpius':[141,'<'],'n Chamaeleontis':[94, 'h'], 
           'Upper Centaurus Lupus':[130, 'X'], 
               }
    marker = ['^',  'p', 'P',  'h','3', 'x', '4',
                 'X', '<', 'd',  '8', '>']
    cluster_all, marker_all = cplot.plot_cluster(dist=False)
    col = cm.jet(np.linspace(0.2, 1, len(cluster_all)-1))
    for i, name in enumerate(list(cluster.keys())):
        ax = axi[i]
        ind = np.where(np.array(cluster_all)==name)[0][0]
        color=col[ind]
        
        CLUSTER_DB = df_caspar.loc[df_caspar['Main Association'].str.contains(name)]

        meanage = int(CLUSTER_DB['Age'].median())
        if name == 'Lagoon Nebula':
            meanage = 0.5
        ax.text(0.95, 0.05, name+', '+str(meanage)+' Myr', ha='right',transform=ax.transAxes, fontsize=6)
        
        NAMES = CLUSTER_DB['Unique Name'].unique()
        
        cplot.plot_graylines(ax, CLUSTER_DB)
         
        stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] = csort.CASPAR_separateByMass(CLUSTER_DB)
        
        cplot.plot_MMdot(ax, Sreg['log Mass'], Sreg['log Mdot'], Supp['log Mass'], Supp['log Mdot'], color_cmap=False, s=16, color=[color], marker=marker[i], edgecolor='none', alpha=0.7, zorder=10)

        cplot.plot_MMdot(ax, BDreg['log Mass'], BDreg['log Mdot'], BDupp['log Mass'], BDupp['log Mdot'], color_cmap=False, s=16, color=[color], marker=marker[i], edgecolor='none', alpha=0.7, zorder=10)

        cplot.plot_MMdot(ax, Preg['log Mass'], Preg['log Mdot'], Pupp['log Mass'], Pupp['log Mdot'], color_cmap=False, s=16, color=[color], marker=marker[i], edgecolor='k', alpha=0.7, zorder=10)
        
        cplot.plot_MMdotFeatures(fig, ax, fs=fs,fs_BL=4, **kwargs)
        

        ############## LINE FITTING ######################
        if lit_lines:
            cplot.plot_litlines(ax)

        else:
            if len(stars['log Mass'].values)>1:
                cfitmmdot.run_statsmodels(stars['log Mass'].values, stars['log Mdot'].values, plot=True, ax=ax, extend=True, color='k')
            if len(bds['log Mass'].values)>1:
                cfitmmdot.run_statsmodels(bds['log Mass'].values, bds['log Mdot'].values, plot=True, ax=ax, extend=True, color='k')
        ##################################################################
        
        if i%2==0:
            ax.set_ylabel('log $\dot{M}$ ($M_{\odot}\ \mathrm{yr}^{-1}$)', fontsize = fs)
        else:
            ax.set_ylabel('')
            
    axi[8].set_xlabel('log Mass ($M_{\odot}$)', fontsize = fs)
    axi[9].set_xlabel('log Mass ($M_{\odot}$)', fontsize = fs)
    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=150, transparent=False)
    plt.show()
