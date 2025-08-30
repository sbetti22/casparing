import pandas as pd
import numpy as np
import .rederiveaccretion as cra


def derive_newMdot(row, unc):

    updated_vals, logaccLum, AccRate, scaleRef, lineFlux = cra.CASPAR_derive_newMdot(row, unc=unc).rederiveRun()
    
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


def put_in_caspar_order(df):
    cc = pd.read_csv('../data/accretion_data/caspar_columns.txt', header=None)
    caspar_cols = cc[0].values
    finaldf = pd.DataFrame(columns=caspar_cols)
    for column in list(df.columns):
        if column in caspar_cols:
            finaldf[column] = df[column]
        else:
            print(f'{column} not in caspar.')
    return finaldf

def save_df(df, savename):
    df.to_csv(savename, index=False)
    
    
    
    
    
    
    
