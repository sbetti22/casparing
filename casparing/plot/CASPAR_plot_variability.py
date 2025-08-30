import scipy.stats as stats
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

import casparing.plot.CASPAR_util as cutil
import casparing.plot.CASPAR_plotting as cplot


def plot_allvariability(df_caspar, **kwargs):
    fs = kwargs.get('fs', 14)
    color = cutil.ObjectColor()
    fig, (ax0, ax00)= plt.subplots(figsize=(4,3), ncols=2,nrows=1, dpi=150, sharey=True,
                                    gridspec_kw={'width_ratios':[3,1]})
    
    lines = cutil.AccDiagLines()['optical-NIR Line Flux']
    lineAccRate = [line + ' Accretion Rate' for line in lines]
    colnames = ['Simbad-Resolvable Name', 'Accretion Rate']  + lineAccRate
    
    dupobjs_df = df_caspar[colnames]
    dupobjs_df[colnames[1:]] = np.log10(dupobjs_df[colnames[1:]])

    dupobjs_df_range = dupobjs_df.set_index('Simbad-Resolvable Name').stack().groupby(level=0).agg([np.ptp])
    dupobjs_df_range = dupobjs_df_range[dupobjs_df_range.ptp != 0]
    
    dupobjs_mass = df_caspar[['Simbad-Resolvable Name', 'log Mass', 'Companion']].groupby('Simbad-Resolvable Name').agg({'log Mass':'mean', 'Companion':pd.Series.mode}).reset_index()

    dupobs_fin = pd.merge(dupobjs_df_range, dupobjs_mass, 
                           on='Simbad-Resolvable Name', how='left')
    
    dupobs_fin['Obj'] = np.where(dupobs_fin['log Mass'] > np.log10(0.075), 'S', np.where(dupobs_fin['Companion'] == 'COM', 'P', 'BD')) 
    
    sns.scatterplot(x="log Mass", y="ptp", data=dupobs_fin, style='Obj', hue='Obj', markers=['D', 'o', 's'], palette = {'S':color[0], 'BD':color[1], 'P':color[2]}, ax=ax0, legend=False)

    cplot.plot_objtype_legend(ax0, color[0], color[1], color[2], fs=fs, loc='upper left', labels_color=False)

    ax00.hist(dupobs_fin['ptp'].loc[dupobs_fin['Obj']=='S'], orientation='horizontal', color='k', alpha=0.8, 
            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='k')
    ax00.hist(dupobs_fin['ptp'].loc[dupobs_fin['Obj']=='P'], orientation='horizontal', color='magenta', alpha=0.4, 
            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='magenta')

    ax00.hist(dupobs_fin['ptp'].loc[dupobs_fin['Obj']=='BD'], orientation='horizontal', color='cyan', alpha=0.6, 
            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='cyan')    
    
    cplot.plot_MMdotFeatures(fig, ax0, fs=fs, burning_lines=True, burning_lines_labels=True, tickparams=True, xlim=(-2.5, 0.3), ylim=(-0.2, 6), xlabel= 'log Mass ($M_\odot$)', ylabel='$\Delta \dot{M}$ (dex)')
    
    cplot.plot_MMdotFeatures(fig, ax00, fs=fs, burning_lines=False, burning_lines_labels=False, tickparams=True, xlim=None, ylim=(-0.2, 6), xlabel= 'PDF', ylabel=None)

    plt.subplots_adjust(wspace=0)
    
    
#def plot_emissionline_variability(df_caspar, **kwargs):
#    fs = kwargs.get('fs', 14)
#    color = cutil.ObjectColor()
#    fig, (ax0, ax00)= plt.subplots(figsize=(4,3), ncols=2,nrows=1, dpi=150, sharey=True,
#                                    gridspec_kw={'width_ratios':[3,1]})
    
    

#def Fig14(bd_update, fs=12):
#    all_names = bd_update['Unique Name'].unique()
#    lines = ['Ha', 'Hb', 'Hg', 'Hd', 'Hep', 'H8', 'H9', 'H10', 'H11', 'H12',
#                    'H13', 'H14', 'H15', 'PaB', 'PaG', 'PaD', 'Pa8', 'Pa9', 'Pa10', 'BrG', 'Br8', 'PfB','HeI 402.6', 'HeI 447.1', 'HeI 471.3', 'HeI 501.6', 'HeI 587.6','HeI 667.8', 'HeI 706.5',
#                    'HeII 468.6', 'CaII K 393.4', 'CaII H 396.9', 'CaII 854.2', 'CaII 866.2', 'NaI 588.9','NaI 589.6', 
#                    'OI 844.6']
#
#    fig, ((ax0, ax00), (ax1,ax11),(ax2, ax22))= plt.subplots(figsize=(4,10), ncols=2, 
#                                                            nrows=3, dpi=150, 
#                                                            gridspec_kw={'width_ratios':[3,1]})
#
#    uniquedf = bd_update[['Unique Name', 'Mass']].loc[(bd_update['Duplicate #']==2)]
#    uniquedf = uniquedf.sort_values(by='Mass', ascending=False)
#    unique = uniquedf['Unique Name']
#    Srange = []
#    BDrange = []
#    Prange=[]
#    for i in np.arange(len(unique.values)):
#        name = unique.values[i]
#        udf = bd_update.loc[bd_update['Unique Name']==name]
#        all_lines = []
#        for line in lines:
#            a = udf[f'{line} Accretion Rate']
#            all_lines.append(np.log10(a.values))
#                
#        mdot = np.log10(udf['Accretion Rate']) 
#        all_lines.append(mdot)
#        mass = udf['Mass'].values[0]
#        COM = udf['Companion'].values[0]
#        rang = np.nanmax(all_lines) - np.nanmin(all_lines)
#        if mass> 0.075:
#            col = 'k'
#            m = 'o'
#            Srange.append(rang)
#        elif COM == 'COM':
#            col='magenta'
#            m='s'
#            Prange.append(rang)
#        else:
#            col='cyan'
#            m='D'
#            BDrange.append(rang)
#        ax0.plot(np.log10(mass), rang, marker=m, color=col, ms=6, alpha=0.6)
#
#    print(np.nanmedian(Srange))
#    subrange = BDrange+Prange
#    print(np.nanmedian(subrange))
#
##    print('scipy t = ', stats.ttest_ind(Srange, subrange, equal_var = False))
#
#    t = (np.nanmean(subrange)-np.nanmean(Srange)) / np.sqrt( (np.var(subrange)/len(subrange)) 
#                                                            + (np.var(Srange)/len(Srange)))
#
#    print('stellar', np.nanmean(Srange), np.nanstd(Srange))
#    print('sub', np.nanmean(subrange), np.nanstd(subrange))
#    print('t test', t)
#
#    N = ( (np.var(subrange)**2 / len(subrange)) + (np.var(Srange)**2./len(Srange))  )**2.
#    D = (((np.var(subrange)**2 / len(subrange)) ** 2.)/(len(subrange)-1)) + (((np.var(Srange)**2 / len(Srange)) ** 2.)/(len(Srange)-1))
#    dof = N/D
#    print('dof', dof)
#            
#    ax00.hist(Srange, orientation='horizontal', color='k', alpha=0.8, 
#            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='k')
#    ax00.hist(Prange, orientation='horizontal', color='magenta', alpha=0.4, 
#            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='magenta')
#
#    ax00.hist(BDrange, orientation='horizontal', color='cyan', alpha=0.6, 
#            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='cyan')
#
#
#
#    finlines = []
#    Srange = []
#    BDrange = []
#    Prange=[]
#    m = ['o', 'D', 's']
#    col = ['k', 'cyan', 'magenta']
#    lab = ['Star', 'Brown Dwarf', 'Planetary Mass\nCompanion']
#    for i in np.arange(3):
#        ax1.plot([],[], marker=m[i], color=col[i], ms=4, alpha=0.6, label=lab[i], linestyle='none')
#        ax0.plot([],[], marker=m[i], color=col[i], ms=4, alpha=0.6, label=lab[i], linestyle='none')
#        
#    num = 0
#    for index, row in bd_update.iterrows():
#        all_lines = []
#        for line in lines:
#            a = row[f'{line} Accretion Rate']
#            if np.isfinite(a):
#                all_lines.append(np.log10(a))
#        if len(all_lines) > 1:
#            rang = np.nanmax(all_lines) - np.nanmin(all_lines)
#            mdot2 = all_lines - np.nanmedian(all_lines)
#            if row['Mass']> 0.075:
#                col = 'k'
#                m = 'o'
#                Srange.append(rang)
#            elif row['Companion'] == 'COM':
#                col='magenta'
#                m='s'
#                Prange.append(rang)
#            else:
#                col='cyan'
#                m='D'
#                BDrange.append(rang)
#            ax1.plot(np.log10(row['Mass']), rang, marker=m, color=col, ms=6, alpha=0.6)
#            if rang > 5:
#                print(row['Unique ID'])
#            finlines.append(rang)
#            num+=1
#        else:
#            finlines.append(np.nan)
#    print('fin num', num)
#    forVAR = bd_update.copy()
#    forVAR['varrange'] = finlines
#    for name in all_names:
#        dfname = forVAR.loc[forVAR['Unique Name']==name]
#        XX = dfname['log Mass'].values
#        XX = np.append(XX, XX[0])
#        YY = dfname['varrange'].values
#        YY = np.append(YY,YY[0])
#        ax1.plot(XX, YY, '--', color='gray', lw=0.5, zorder=-1)
#        
#    print(np.nanmedian(Srange))
#    subrange = BDrange+Prange
#    print(np.nanmedian(subrange))
#    ax11.hist(Prange, orientation='horizontal', color='magenta', alpha=0.4, 
#            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='magenta')
#
#    ax11.hist(Srange, orientation='horizontal', color='k', alpha=0.8, 
#            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='k')
#    ax11.hist(BDrange, orientation='horizontal', color='cyan', alpha=0.6, 
#            histtype='stepfilled', bins=np.arange(0,6.5,0.3), density=True, edgecolor='cyan')
#
#
#
#    COL = ['darkred', 'lightgreen', 'darkgreen', 'lightblue', 'purple', 'violet']
#    sym = ['o', 's', 'D', 'p', '^', '>']
#    # plot 2 = 1 object, 1 emission line, all epochs
#
#    for k, line in enumerate(['Hα', 'Paβ', 'Paγ', 'Brγ', 'CaII K', 'CaII H']):
#        ax2.plot([],[],marker=sym[k], color=COL[k], ms=4, label=line, linestyle='none')
#
#    LL = ['Ha', 'PaB', 'PaG', 'BrG', 'CaII K 393.4', 'CaII H 396.9']
#    harange, pab, pag, brg, cak, cah = [],[],[],[],[],[]
#    hist_vals = [harange, pab, pag, brg, cak, cah]
#    mass_vals = [[],[],[],[],[],[]]
#    for i in np.arange(len(all_names)):
#        a = bd_update.loc[bd_update['Unique Name']==all_names[i]]
#        labels = []
#        masses = []
#        
#        for k, line in enumerate(['Ha', 'PaB', 'PaG', 'BrG', 'CaII K 393.4', 'CaII H 396.9']):
#            b = a[f'{line} Accretion Rate'].loc[a[f'{line} Accretion Rate'].notnull()]
#            m = a['Mass'].loc[a[f'{line} Accretion Rate'].notnull()]
#            if (len(b) > 1) & (len(b.unique())>1):
#                bb = b.unique()
#                mdot = np.log10(bb)
#                
#                mdot2 = mdot - np.nanmedian(mdot)
#                mdotrange = mdot.max() - mdot.min()
#                hist_vals[k].append(mdotrange)
#                mass_vals[k].append(np.log10(m.values[0]))
#                ax2.plot(np.log10(m.values[0]), mdotrange, marker=sym[k], color=COL[k], 
#                        zorder=10, ms=6)
#                
#    #             violin_parts = ax2.violinplot(mdot2, positions=[pos],
#    #                                           widths=1,showmedians=True)
#                labels.append(mdotrange)
#                masses.append(np.log10(m.values[0]))
#        
#        ax2.plot(masses, labels, '--', color='gray', zorder=1, linewidth=0.5)
#
#    for num, val in enumerate(hist_vals): 
#        print(LL[num],len(val), len(mass_vals[num]))
#        IND = np.where(np.array(mass_vals[num]) > -1.12)
#        IND2 = np.where(np.array(mass_vals[num]) <= -1.12)
#        
#        print(np.nanmedian(np.array(val)), np.nanmedian(np.array(val)[IND]), np.nanmedian(np.array(val)[IND2]))
#        
#        if LL[num] == 'PaB':
#            from scipy import stats
#            ste = np.array(val)[IND]
#            sub = np.array(val)[IND2]
#            print('scipy t = ', stats.ttest_ind(ste, sub, equal_var = False))
#            
#            t = (np.nanmean(sub)-np.nanmean(ste)) / np.sqrt( (np.var(sub)/len(sub)) + (np.var(ste)/len(ste)))
#            print('t test', t)
#            
#            N = ( (np.var(sub)**2 / len(sub)) + (np.var(ste)**2./len(ste))  )**2.
#            D = (((np.var(sub)**2 / len(sub)) ** 2.)/(len(sub)-1)) + (((np.var(ste)**2 / len(ste)) ** 2.)/(len(ste)-1))
#            dof = N/D
#            print('dof', dof)
#            
#        ax22.hist(val, orientation='horizontal', color=COL[num], alpha=0.8, 
#                histtype='stepfilled', bins=np.arange(0,1.5,0.08), density=True, edgecolor=COL[num])
#
#
#
#    ax0.set_title('           a) range of $\dot M$ for\n              measured accretion flux in all epochs')
#    ax1.set_title('            b) range of $\dot M$ for all\n               measured emission lines in 1 epoch')
#    ax2.set_title('              c) range of $\dot M$ for 1 emission\n                line measured in multiple epochs')
#    
#    for axi in [ax00,ax11, ax22]:
#        axi.set_yticklabels([])
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
#        axi.set_xlabel('PDF', fontsize=fs)
#        
#    for axi in [ax0, ax1, ax2]:
#        axi.set_xlim(-2.5, 0.3) 
#        axi.set_xlabel('log $M$ ($M_\odot$)', fontsize=fs)
#        axi.set_ylabel(r'$\Delta \dot M$ (dex)', fontsize=fs)    
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
#        fs_BL=fs-6
#        axi.text(-1.124938736, 0.95, 'Hydrogen\nBurning Limit', color='dimgray', rotation=90, 
#                fontsize=fs_BL, ha='center', transform=axi.get_xaxis_transform(), va='top')
#        
#        axi.text(-1.906371, 0.95, 'Deuterium\nBurning Limit', color='dimgray', rotation=90, 
#                fontsize=fs_BL, ha='center', transform=axi.get_xaxis_transform(), va='top')
#
#        axi.axvline(x = -1.124938736, color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
#        axi.axvline(x = -1.906371, color = 'dimgray', linestyle='-', alpha=0.4, linewidth=0.5)
#
#    ax0.set_ylim(-0.2, 6)
#    ax00.set_ylim(-0.2, 6)
#    ax1.set_ylim(0,7.)
#    ax11.set_ylim(0,7.)
#    ax2.set_ylim(-0.1,1.8)
#    ax22.set_ylim(-0.1,1.8)
#
#
#
#    ax0.legend(loc='upper left', ncol=1, fontsize=fs-4)
#    ax1.legend(loc='upper left', ncol=1, fontsize=fs-4)
#    ax2.legend(loc='upper left', ncol=2, fontsize=fs-4)
#    plt.tight_layout()
#    plt.subplots_adjust(wspace=0)
#    plt.savefig('/Users/sbetti/Desktop/CASPAR_plotslagoon/variability_V2.pdf', dpi=150)
#    plt.show()


#def variability(bd_update, **kwargs):
#    uniquedf = bd_update[['Unique Name', 'Mass']].loc[(bd_update['Duplicate #']==2)&(bd_update['Companion']!='COM')]
#    uniquedf = uniquedf.sort_values(by='Mass', ascending=False)
#    unique = uniquedf['Unique Name']
#
#    stellar = []
#    substellar = []
#    fs = 10
#    binsize=0.15
#    red = '#882255'
#    blue = '#88CCEE'
#    fig, (ax1, ax2) = plt.subplots(figsize=(4,6), ncols=1, nrows=2, dpi=150)
#
#
#    for i in np.arange(len(unique.values)):
#        name = unique.values[i]
#        udf = bd_update.loc[bd_update['Unique Name']==name]
#        mdot = np.log10(udf['Accretion Rate']) 
#
#        rang = np.nanmax(mdot) - np.nanmin(mdot)
#        mdot2 = mdot - np.nanmedian(np.log10(udf['Accretion Rate']))
#        violin_parts = ax2.violinplot(mdot2, positions=[np.log10(udf['Mass'].values[0])], widths=0.07,
#                showmedians=True)
#        if udf['Mass'].values[0] > 0.075:
#            stellar.append(rang)
#            for partname in ('cbars','cmins','cmaxes','cmedians'):
#                vp = violin_parts[partname]
#                vp.set_edgecolor('darkblue')
#                vp.set_linewidth(1)
#
#            # Make the violin body blue with a red border:
#            for vp in violin_parts['bodies']:
#                vp.set_facecolor(blue)
#                vp.set_edgecolor(blue)
#                vp.set_linewidth(1)
#                vp.set_alpha(0.6)
#        else:
#            substellar.append(rang)
#            for partname in ('cbars','cmins','cmaxes','cmedians'):
#                vp = violin_parts[partname]
#                vp.set_edgecolor('darkred')
#                vp.set_linewidth(1)
#
#            # Make the violin body blue with a red border:
#            for vp in violin_parts['bodies']:
#                vp.set_facecolor(red)
#                vp.set_edgecolor(red)
#                vp.set_linewidth(1)
#                vp.set_alpha(0.3)
#
#    ax2.set_xlim(-1.85, 0.3)
#    ax2.set_xlabel('log $M$ ($M_\odot$/yr)', fontsize=fs)
#    ax2.set_ylabel(r'log $\dot M/ \overline{\dot M}$ (dex)', fontsize=fs)
#
#    print(np.nanmedian(stellar), np.nanstd(stellar))
#    print(np.nanmedian(substellar), np.nanstd(substellar))
#
#    bins=np.arange(min(stellar), max(stellar) + binsize, binsize)
#    entries, bin_edges, patches = ax1.hist(stellar, bins=bins, density=True, color=blue, alpha=0.5, label='stellar')
#    ax1.hist(stellar, bins=bins, density=True, color=blue, histtype='step', linewidth=2)
#
#    bins=np.arange(min(substellar), max(substellar) + binsize, binsize)
#    entries, bin_edges, patches = ax1.hist(substellar, bins=bins, density=True, color=red, alpha=0.3, label='substellar')
#    ax1.hist(substellar, bins=bins, density=True, color=red,histtype='step', linewidth=2)
#    ax1.set_xlabel('range of multi-epoch $\dot M$ (dex)', fontsize=fs)
#    ax1.set_ylabel('N', fontsize=fs)
#
#    ax1.legend(fontsize=fs)
#
#    for axi in [ax1, ax2]:
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
#
#    plt.tight_layout()
#    if 'savefig' in kwargs:
#        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
#    plt.show()
#    
#    
#    
#    
#def variability_V2(bd_update, **kwargs):
#    
#    # plot 1 = 1 object, 1 epoch, all lines
#    # plot 2 = 1 object, 1 emission line, all epochs
#    
#    all_names = bd_update['Unique Name'].unique()
##    for i in np.arange(len(all_names)):
#        
#    uniquedf = bd_update[['Unique Name', 'Mass']].loc[(bd_update['Duplicate #']==2)&(bd_update['Companion']!='COM')]
#    uniquedf = uniquedf.sort_values(by='Mass', ascending=False)
#    unique = uniquedf['Unique Name']
#
#    stellar = []
#    substellar = []
#    fs = 10
#    binsize=0.15
#    red = '#882255'
#    blue = '#88CCEE'
#    fig, (ax1, ax2) = plt.subplots(figsize=(4,6), ncols=1, nrows=2, dpi=150)
#
#
#    for i in np.arange(len(unique.values)):
#        name = unique.values[i]
#        udf = bd_update.loc[bd_update['Unique Name']==name]
#        mdot = np.log10(udf['Accretion Rate']) 
#
#        rang = np.nanmax(mdot) - np.nanmin(mdot)
#        mdot2 = mdot - np.nanmedian(np.log10(udf['Accretion Rate']))
#        violin_parts = ax2.violinplot(mdot2, positions=[np.log10(udf['Mass'].values[0])], widths=0.07,
#                showmedians=True)
#        if udf['Mass'].values[0] > 0.075:
#            stellar.append(rang)
#            for partname in ('cbars','cmins','cmaxes','cmedians'):
#                vp = violin_parts[partname]
#                vp.set_edgecolor('darkblue')
#                vp.set_linewidth(1)
#
#            # Make the violin body blue with a red border:
#            for vp in violin_parts['bodies']:
#                vp.set_facecolor(blue)
#                vp.set_edgecolor(blue)
#                vp.set_linewidth(1)
#                vp.set_alpha(0.6)
#        else:
#            substellar.append(rang)
#            for partname in ('cbars','cmins','cmaxes','cmedians'):
#                vp = violin_parts[partname]
#                vp.set_edgecolor('darkred')
#                vp.set_linewidth(1)
#
#            # Make the violin body blue with a red border:
#            for vp in violin_parts['bodies']:
#                vp.set_facecolor(red)
#                vp.set_edgecolor(red)
#                vp.set_linewidth(1)
#                vp.set_alpha(0.3)
#
#    ax2.set_xlim(-1.85, 0.3)
#    ax2.set_xlabel('log $M$ ($M_\odot$/yr)', fontsize=fs)
#    ax2.set_ylabel(r'log $\dot M/ \overline{\dot M}$ (dex)', fontsize=fs)
#
#    print(np.nanmedian(stellar), np.nanstd(stellar))
#    print(np.nanmedian(substellar), np.nanstd(substellar))
#
#    bins=np.arange(min(stellar), max(stellar) + binsize, binsize)
#    entries, bin_edges, patches = ax1.hist(stellar, bins=bins, density=True, color=blue, alpha=0.5, label='stellar')
#    ax1.hist(stellar, bins=bins, density=True, color=blue, histtype='step', linewidth=2)
#
#    bins=np.arange(min(substellar), max(substellar) + binsize, binsize)
#    entries, bin_edges, patches = ax1.hist(substellar, bins=bins, density=True, color=red, alpha=0.3, label='substellar')
#    ax1.hist(substellar, bins=bins, density=True, color=red,histtype='step', linewidth=2)
#    ax1.set_xlabel('range of multi-epoch $\dot M$ (dex)', fontsize=fs)
#    ax1.set_ylabel('N', fontsize=fs)
#
#    ax1.legend(fontsize=fs)
#
#    for axi in [ax1, ax2]:
#        axi.tick_params(which='major', direction='in', top=True, right=True,left=True, bottom=True, labelsize=fs)
#        axi.tick_params(which='minor', direction='in', top=True,bottom=True, right=True, left=True, labelsize=fs)
#        axi.minorticks_on()
#
#    plt.tight_layout()
#    if 'savefig' in kwargs:
#        plt.savefig(kwargs['savefig'], dpi=300,transparent=False)
#    plt.show()
#    
