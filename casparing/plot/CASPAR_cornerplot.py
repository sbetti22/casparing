'''
plotting corner plot of physical parameters
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns

from casparing.plot.CASPAR_util import AccDiagColor
color = AccDiagColor()

def corner_map(df_caspar, fs=20):
    '''
    Function to plot a corner plot of the physical parameters in CASPAR (Association, Spectral Type, Age, Mass, Luminosity, Accretion Luminosity, Accretion Rate, and Accretion Diagnostic)
    
    Parameters
    -------
    df_caspar: pandas.DataFrame
        pandas Dataframe of CASPAR or Literature Database after opened with casparing.CASPAR_sortdata.CASPAR_loaddata() 
    fs: int
        fontsize of legend labels
        
    
    Returns
    -------
    NoneType
        None
    
    '''
    ###############
    name = ['Main Association', 'Sp Type Num', 'Age', 'log Mass', 'log Lum',  'Log Accretion Luminosity', 'log Mdot',  'AD']
    data = df_caspar[name]

    data = data.set_index('AD', drop=True)
    sorter = list(color.keys())
    data = data.loc[sorter[::-1]]
    data = data.reset_index()
    #################
    AD = data['AD'].unique()

    pp = sns.pairplot(data, corner=True, hue='AD', diag_kind='kde', palette=color, 
                      plot_kws=dict(alpha=0.4, linewidth=1, s=60, edgecolor='k'), kind='scatter', 
                      hue_order=AD, height=3)

    def plot_extra(x, y, **kwargs):
        sns.scatterplot(x=df_caspar[x.name], y=df_caspar[y.name], color='gray', zorder=-1, alpha=0.3, marker='x')
    pp.map_offdiag(plot_extra)


    Name = ['Association', 'Spectral Type', 'Age (Myr)', r'$\log(M/M_\odot$)', r'log $L\ (L_\odot$)', r'log($L_\mathrm{acc}/L_\odot$)', r'log($\dot{M}/(M_\odot\ \mathrm{yr}^{-1}))$']
    print('AXES LABELS ARE OFF IN CODE!!')
#    for i in range(6):
#        pp.axes[5,i].set_xlabel(Name[i], fontsize=fs)
#        pp.axes[5,i].tick_params(which='both', direction='in', labelsize=fs)
#        if i > 0:
#            pp.axes[i,0].set_ylabel(Name[i], fontsize=fs)
#            pp.axes[i,0].tick_params(which='both', direction='in', labelsize=fs)

#    # Age
#    pp.axes[5,2].set_xticks([0,20,40,60])
#    pp.axes[2,0].set_yticks([0,20,40,60])
#    pp.axes[5,2].set_xlim([-5, 62])
#    pp.axes[2,0].set_ylim([-5, 62])
#
#    # MASS
#    pp.axes[5,3].set_xticks([-2,-1,0,1])
#    pp.axes[3,0].set_yticks([-2,-1,0,1])
#    pp.axes[5,3].set_xlim([-2.5, 1.2])
#    pp.axes[3,0].set_ylim([-2.5, 1.2])
#
#    # Luminosity
#    # pp.axes[5,4].set_xticks([0,1,2,3,4,5])
#    # pp.axes[4,0].set_yticks([0,1,2,3,4,5])
#    pp.axes[5,4].set_xlim([-3.5,1.255])
#    pp.axes[4,0].set_ylim([-3.5,1.255])
#
#    # Lacc
#    pp.axes[5,5].set_xticks([-7, -5, -3, -1, 1 ])
#    pp.axes[5,0].set_yticks([-7, -5, -3, -1, 1])
#    pp.axes[5,5].set_xlim([-7.5,1.5])
#    pp.axes[5,0].set_ylim([-7.5,1.5])
#
#    # Spectral Type
#    pp.axes[5,1].set_xticks([60,70,80])
#    pp.axes[5,1].set_xticklabels(['K0', 'M0', 'L0'])
#    pp.axes[1,0].set_yticks([60,70,80])
#    pp.axes[1,0].set_yticklabels(['K0', 'M0', 'L0'])
#    pp.axes[5,1].set_xlim([57, 83])
#    pp.axes[1,0].set_ylim([57, 83])
#
#    # Association
#    pp.axes[5,0].set_xlim(0,20.5)
#    pp.axes[5,0].set_xticks(np.arange(1,20))
#    pp.axes[5,0].set_xticklabels(['Sh 2-284','Lagoon Nebula','Upper Centaurus Lupus','Tucana-Horologium','IC 348','Taurus','118 Tau','σ Ori','NGC 2024','Chamaeleon I','TW Hya','Upper Scorpius','Lupus','Argus','p Oph','Corona-Australis', 'n Chamaeleontis','25 Orionis', 'β Pictoris'], rotation='vertical', fontsize=11)
#    pp.axes[5,2].set_xticks([0,20,40,60])
#    pp.axes[2,0].set_yticks([0,20,40,60])
#    pp.axes[5,2].set_xlim([-5, 62])
#    pp.axes[2,0].set_ylim([-5, 62])
#
#    # MASS
#    pp.axes[5,3].set_xticks([-2,-1,0,1])
#    pp.axes[3,0].set_yticks([-2,-1,0,1])
#    pp.axes[5,3].set_xlim([-2.5, 1.2])
#    pp.axes[3,0].set_ylim([-2.5, 1.2])
#
#    # Luminosity
#    # pp.axes[5,4].set_xticks([0,1,2,3,4,5])
#    # pp.axes[4,0].set_yticks([0,1,2,3,4,5])
#    pp.axes[5,4].set_xlim([-3.5,1.255])
#    pp.axes[4,0].set_ylim([-3.5,1.255])
#
#    # Lacc
#    pp.axes[5,5].set_xticks([-7, -5, -3, -1, 1 ])
#    pp.axes[5,0].set_yticks([-7, -5, -3, -1, 1])
#    pp.axes[5,5].set_xlim([-7.5,1.5])
#    pp.axes[5,0].set_ylim([-7.5,1.5])
#
#    # Spectral Type
#    pp.axes[5,1].set_xticks([60,70,80])
#    pp.axes[5,1].set_xticklabels(['K0', 'M0', 'L0'])
#    pp.axes[1,0].set_yticks([60,70,80])
#    pp.axes[1,0].set_yticklabels(['K0', 'M0', 'L0'])
#    pp.axes[5,1].set_xlim([57, 83])
#    pp.axes[1,0].set_ylim([57, 83])
#
#    # Association
#    pp.axes[5,0].set_xlim(0,20.5)
#    pp.axes[5,0].set_xticks(np.arange(1,20))
#    pp.axes[5,0].set_xticklabels(['Sh 2-284','Lagoon Nebula','Upper Centaurus Lupus','Tucana-Horologium','IC 348','Taurus','118 Tau','σ Ori','NGC 2024','Chamaeleon I','TW Hya','Upper Scorpius','Lupus','Argus','p Oph','Corona-Australis', 'n Chamaeleontis','25 Orionis', 'β Pictoris'], rotation='vertical', fontsize=11)

    sns.move_legend(pp, 'upper right', bbox_to_anchor=(0.85, 0.98), 
                    title='Accretion Diagnostic', fontsize=fs, frameon=True, 
                    title_fontsize=20)
    plt.savefig('test.png', dpi=300)
    plt.show()

##############################################################################



    















