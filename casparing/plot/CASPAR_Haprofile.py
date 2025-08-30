import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def Haprofile(df_lit, df_caspar, **kwargs):
    
    bd_orig03 = df_lit.loc[(df_caspar['Original Reference']=='Muzerolle 2003') & (df_caspar['Accretion Diagnostic']=='Line Profile')]
    bd_update03 = df_caspar.loc[(df_caspar['Original Reference']=='Muzerolle 2003')& (df_caspar['Accretion Diagnostic']=='Line Profile')]
    
    bd_orig05 = df_lit.loc[(df_caspar['Original Reference']=='Muzerolle 2005') & (df_caspar['Accretion Diagnostic']=='Line Profile')]
    bd_update05 = df_caspar.loc[(df_caspar['Original Reference']=='Muzerolle 2005')& (df_caspar['Accretion Diagnostic']=='Line Profile')]

    bd_origNatta = df_lit.loc[(df_caspar['Original Reference']=='Natta 2004') & (df_caspar['Accretion Diagnostic']=='Line Profile')]
    bd_updateNatta = df_caspar.loc[(df_caspar['Original Reference']=='Natta 2004')& (df_caspar['Accretion Diagnostic']=='Line Profile')]

    bd_origEsp = df_lit.loc[(df_caspar['Original Reference']=='Espaillat2008') & (df_caspar['Accretion Diagnostic']=='Line Profile')]
    bd_updateEsp = df_caspar.loc[(df_caspar['Original Reference']=='Espaillat2008')& (df_caspar['Accretion Diagnostic']=='Line Profile')]

    fig, (ax3, ax1) = plt.subplots(figsize=(4,6.5),dpi=150, nrows=2, ncols=1)
    ax2 = ax3.twinx().twiny()
    ax4 = ax3.twinx()

    s = 70
    mass_ratio1 = bd_update03['Mass'] / bd_orig03['Mass']
    ax1.scatter(bd_update03['Mass'], mass_ratio1, color='k', edgecolor='white', s=s, marker='D',zorder=-1,label='Muzerolle+2003')

    sptype_ratio1 = abs(bd_update03['Sp Type Num'] - bd_orig03['Sp Type Num'])
    ax2.scatter(bd_update03['Sp Type Num'], sptype_ratio1, color='k',edgecolor='white', s=s,zorder=-1, marker='D',label='Muzerolle+2003')

    mass_ratio2 = bd_update05['Mass'] / bd_orig05['Mass']
    ax1.scatter(bd_update05['Mass'], mass_ratio2, color='k', s=s-30,  zorder=10, marker='o', label='Muzerolle+2005')

    temp_ratio1 = abs(bd_update05['Teff'] - bd_orig05['Teff'])   
    ax3.scatter(bd_update05['Teff'], temp_ratio1, color='k', marker='o', s=s-30,zorder=10, label='Muzerolle+2005')

    mass_ratio3 = bd_updateNatta['Mass'] / bd_origNatta['Mass']
    ax1.scatter(bd_updateNatta['Mass'], mass_ratio3, color='k', s=s-30,  zorder=10, marker='s', label='Natta+2004')

    temp_ratio2 = abs(bd_updateNatta['Teff'] - bd_origNatta['Teff'])   
    ax3.scatter(bd_updateNatta['Teff'], temp_ratio2, color='k', marker='s', s=s-30,zorder=10, label='Natta+2004')

    mass_ratio4 = bd_updateEsp['Mass'] / bd_origEsp['Mass']
    ax1.scatter(bd_updateEsp['Mass'], mass_ratio4, color='k', edgecolor='white', s=s, zorder=10, marker='^',label='Espaillat+2008')

    temp_ratio3 = abs(bd_updateEsp['Teff'] - bd_origEsp['Teff'])   
    ax3.scatter(bd_updateEsp['Teff'], temp_ratio3, color='k', marker='^', s=s-30,zorder=10,label='Espaillat+2008')
    
    allvals = []
    for val in [mass_ratio1.values, mass_ratio2.values, mass_ratio3.values, mass_ratio4.values]:
        for i in val:
            allvals.append(i)
    
    allvals = []
    for val in [temp_ratio1.values, temp_ratio2.values, temp_ratio3.values]:
        for i in val:
            allvals.append(i)
    for i in [200, 0,0,0,0,0]:
        allvals.append(0)
    

    ax1.set_xlim(0.009, 0.5)
    ax1.set_ylim(0, 2)
    ax1.set_xlabel('log($M_\mathrm{CASPAR}$/$M_\odot$)')
    ax1.set_ylabel('$M_\mathrm{CASPAR}/M_\mathrm{literature}$') # 
    ax1.set_xscale('log')
    ax1.set_xticks([0.01, 0.1])
    ax1.set_xticklabels(['0.01', '0.1'])

    ax2.set_xlim(2.5, 9.5)
    ax2.set_xticks([3,4,5, 6, 7,8,9])
    ax2.set_ylim(-0.1,1.2)
    ax2.set_xticklabels(['M3', 'M4','M5', 'M6', 'M7', 'M8', 'M9'])
    ax2.set_xlabel('CASPAR spectral type', color='k')


    ax3.set_xlim(3485, 2520)
    ax3.set_ylim(-20, 240)
    ax3.set_xlabel('T$_\mathrm{CASPAR}$ (K)')
    ax3.set_ylabel('ΔT (CASPAR - literature)')

    ls = 10
    ax1.tick_params(which='major', direction='in', top=True, right=True, length=5,labelsize=ls)
    ax1.tick_params(which='minor', direction='in', top=True, right=True, length=3,labelsize=ls)
    ax1.minorticks_on()

    ax2.tick_params(which='major', direction='in', top=True, right=False, length=5,labelsize=ls,  labelcolor='k', color='k')
    ax2.tick_params(which='minor', direction='in', top=True, right=False, length=3,labelsize=ls,  labelcolor='k', color='k')
    ax2.minorticks_on()

    ax3.tick_params(axis='x',which='major', direction='in', top=False, right=False, left=True, bottom=True, length=5,labelsize=ls)
    ax3.tick_params(axis='x',which='minor', direction='in', top=False, right=False, left=True, bottom=True, length=3,labelsize=ls)
    ax3.tick_params(axis='y',which='major', direction='in', top=False, right=False, left=True, bottom=True, length=5,labelsize=ls)
    ax3.tick_params(axis='y',which='minor', direction='in', top=False, right=False, left=True, bottom=True, length=3,labelsize=ls)
    ax3.minorticks_on()

    ax4.set_ylim(-0.1,1.2)
    ax4.tick_params(which='major', direction='in', axis='y', color='k', labelcolor='k', length=5)
    ax4.tick_params(which='minor', direction='in', axis='y', color='k', labelcolor='k', length=3)
    ax4.minorticks_on()
    ax4.set_ylabel('Δsptype (CASPAR - literature)', color='k')

    ax2.set_yticks([])
    
    ax1.legend(loc='lower right', ncol=2, fontsize=7,facecolor='white', framealpha=1)
    handles, labels = ax1.get_legend_handles_labels()
    
    leg = ax2.legend(handles, labels,loc='upper left', ncol=1, fontsize=7, frameon=True, handletextpad=0.5)
    leg.set_zorder(1000)

    plt.tight_layout()

    if 'savefig' in kwargs:
        plt.savefig(kwargs['savefig'], dpi=150, transparent=False) 
    plt.show()
