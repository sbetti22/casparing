import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
from scipy import stats

def linfitMMdot():
    HERE = os.path.dirname(os.path.abspath(__file__))
    fil = os.path.join(HERE, 'Linear_fitting_ALLValues_LINMAX.csv')
    linfit = pd.read_csv(fil, delimiter=',', comment='#', skiprows=[0])

    linfit = linfit.set_index('obj')
    linfit = linfit.replace(-1000, np.nan)
    return linfit

def AccDiagColor():
    magenta = '#CC6677'  #line luminosity
    darkmagenta = '#882255'  # ha photometry
    gold = '#DDCC77' # profile modeling
    forestgreen = '#117733' # 10% width
    lightblue = '#88CCEE' # continuum

    color = {'Line Luminosity':magenta, 'Hα Photometric Luminosity':darkmagenta, 'Line Profile':forestgreen, 'Continuum Excess':lightblue}
    return color

def AccDiagColorDark():
    darkmagenta = '#7a3d47' #'#CC6677'
    superdarkmagenta = '#541535' #'#882255'
    darkgold = '#857a47'# '#DDCC77'
    darkforestgreen = '#0a471f' #'#117733'
    darkblue = '#446677' #'#88CCEE'

    color = {'Line Luminosity':darkmagenta, 'Hα Photometric Luminosity':superdarkmagenta, 'Line Profile':darkforestgreen, 'Continuum Excess':darkblue}
    return color

def AccDiagLines():
    d = {'Hα 10% width':['Ha 10%'],
         'Line Luminosity':['Ha', 'Hb', 'Hg', 'Hd', 'Hep', 'H8', 'H9', 'H10', 'H11', 'H12',
             'H13', 'H14', 'H15', 'PaB', 'PaG', 'PaD', 'Pa8', 'Pa9', 'Pa10', 'BrG', 'Br8', 'PfB','HeI 402.6', 'HeI 447.1', 'HeI 471.3', 'HeI 501.6', 'HeI 587.6','HeI 667.8', 'HeI 706.5',
            'HeII 468.6', 'CaII K 393.4', 'CaII H 396.9', 'CaII 854.2', 'CaII 866.2', 'NaI 588.9','NaI 589.6', 
            'OI 844.6'],
         'optical-NIR Line Flux':['Ha', 'Hb', 'Hg', 'Hd', 'Hep', 'H8', 'H9', 'H10', 'H11', 'H12',
                 'H13', 'H14', 'H15', 'PaB', 'PaG', 'PaD', 'Pa8', 'Pa9', 'Pa10', 'BrG', 'Br8', 'PfB','HeI 402.6', 'HeI 447.1', 'HeI 471.3', 'HeI 501.6', 'HeI 587.6','HeI 667.8', 'HeI 706.5',
                'HeII 468.6', 'CaII K 393.4', 'CaII H 396.9', 'CaII 854.2', 'CaII 866.2', 'NaI 588.9','NaI 589.6', 
                'OI 844.6'],
         'Balmer':['Ha', 'Hb', 'Hg', 'Hd', 'Hep', 'H8', 'H9', 'H10', 'H11', 'H12',
                 'H13', 'H14', 'H15'], 
         'Paschen/Brackett/Pfund':['PaB', 'PaG', 'PaD', 'Pa8', 'Pa9', 'Pa10','BrG', 'Br8', 'PfB'],
         r'He $\mathtt{I}$':['HeI 402.6', 'HeI 447.1', 'HeI 471.3', 'HeI 501.6', 'HeI 587.6','HeI 667.8', 'HeI 706.5'],
         r'Ca $\mathtt{II}$':['CaII K 393.4', 'CaII H 396.9', 'CaII 854.2', 'CaII 866.2']}
    return d


ramp = lambda u: np.maximum( u, 0 )
step = lambda u: ( u > 0 ).astype(float)
def SegmentedLinearReg( X, Y, breakpoints ):
    nIterationMax = 10

    breakpoints = np.sort( np.array(breakpoints) )

    dt = np.min( np.diff(X) )
    ones = np.ones_like(X)

    for i in range( nIterationMax ):
        # Linear regression:  solve A*p = Y
        Rk = [ramp( X - xk ) for xk in breakpoints ]
        Sk = [step( X - xk ) for xk in breakpoints ]
        A = np.array([ ones, X ] + Rk + Sk )
        p =  lstsq(A.transpose(), Y, rcond=None)[0] 

        # Parameters identification:
        a, b = p[0:2]
        ck = p[ 2:2+len(breakpoints) ]
        dk = p[ 2+len(breakpoints): ]

        # Estimation of the next break-points:
        newBreakpoints = breakpoints - dk/ck 

        # Stop condition
        if np.max(np.abs(newBreakpoints - breakpoints)) < dt/5:
            break

        breakpoints = newBreakpoints
    else:
        print( 'maximum iteration reached' )

    # Compute the final segmented fit:
    Xsolution = np.insert( np.append( breakpoints, max(X) ), 0, min(X) )
    ones =  np.ones_like(Xsolution) 
    Rk = [ c*ramp( Xsolution - x0 ) for x0, c in zip(breakpoints, ck) ]

    Ysolution = a*ones + b*Xsolution + np.sum( Rk, axis=0 )

    return Xsolution, Ysolution


def to_SpTyNum(SPTs):
        '''
        This function takes a spectral type in the form - 'letternumber', i.e. 'M5'.
        It then translates that spec type into a spectral type number identifier used for interpolation.
        '''
        spty_dict = {'O' : 1,'B' : 2, 'A' : 3, 'F' : 4, 'G' : 5,'K' : 6, 'M' : 7,'L' : 8}
        sptynum = []
        for SPT in SPTs:
            if isinstance(SPT, str):
                if 'V' in SPT:
                    SPT = SPT.split('V')[0]
                    dec = 0
                if 'e' in SPT:
                    SPT = SPT.split('e')[0]
                    dec = 0
                elif '-' in SPT:
                    SPT = SPT.split('-')[0]
                    dec = 0
                elif '.' in SPT:
                    dec = float('0.' + SPT.split('.')[-1])
                    SPT = SPT.split('.')[0]
                elif '+' in SPT:
                    SPT = SPT.split('+')[0]
                    dec = 0

                else:
                    dec = 0
                letter = SPT[0]
                number = SPT[1:]
                sptynum.append(spty_dict[letter]*10 + int(number) + dec)
            else:
                sptynum.append(np.nan)
        return sptynum


def binn(x, y, nbins=10, plot=True,ax=None, marker='o', color='k' ):
    bin_means, bin_edges, binnumber = stats.binned_statistic(x, y, 'mean', bins=nbins)
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centers = bin_edges[1:] - bin_width/2
    XX = bin_centers
    YY = bin_means
    print('Y', YY)
    bin_std, bin_edges, binnumber = stats.binned_statistic(x, y, 'std', bins=nbins)
    YYerr = bin_std

    if plot:
        if ax is None:
            raise ValueError('must supply ax if plot=True')
        else:
            ax.errorbar(XX, YY, yerr=YYerr, color=color, capsize=0, markersize=6, mec='k', ecolor=color,
                   mew=0.7, lw=0.7, fmt=marker, zorder=100)

    return XX, YY, YYerr
       