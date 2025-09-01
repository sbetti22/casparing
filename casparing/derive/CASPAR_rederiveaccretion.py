import pandas as pd
import numpy as np
from casparing.derive.rederiveaccretion import CASPAR_derive_newMdot

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

def derive_newMdot(row, unc):

    updated_vals, logaccLum, AccRate, scaleRef, lineFlux = CASPAR_derive_newMdot(row, unc=unc).rederiveRun()
    
    if unc:
        logaccLum, logaccLumerr = logaccLum
        AccRate, AccRateerr = AccRate
        row['Log Accretion Luminosity err'] = logaccLumerr
        row['Accretion Rate err'] = AccRateerr

    row['Log Accretion Luminosity'] = logaccLum
    row['Accretion Rate'] = AccRate
    row['Scaling Relation Reference'] = scaleRef

    ## ALL LINES
    for line in updated_vals.keys():
        if unc:
            row[f'{line} Log Accretion Luminosity'] = updated_vals[line][0][0]
            row[f'{line} Log Accretion Luminosity err'] = updated_vals[line][0][1]
            row[f'{line} Accretion Rate'] = updated_vals[line][1][0]
            row[f'{line} Accretion Rate err'] = updated_vals[line][1][1]
        else:
            row[f'{line} Log Accretion Luminosity'] = updated_vals[line][0]
            row[f'{line} Accretion Rate'] = updated_vals[line][1]

    if np.isfinite(lineFlux):
        line = row['Tracer']
        if unc:
            lineFlux, lineFluxerr = lineFlux
            row[f'{line} Line Flux err'] = lineFlux
            row[f'{line} Log Accretion Luminosity err'] = logaccLum
            row[f'{line} Accretion Rate err'] = AccRate
        row[f'{line} Line Flux'] = lineFlux
        row[f'{line} Log Accretion Luminosity'] = logaccLum
        row[f'{line} Accretion Rate'] = AccRate
    return row
        
    
def accretion_properties(df, unc=True):
    df = df.parallel_apply(derive_newMdot, args=(unc,), axis=1, result_type='expand')
    
        
    df.drop([col for col in df.columns if 'Original' in col], axis=1, inplace=True)
    
    return df


    
