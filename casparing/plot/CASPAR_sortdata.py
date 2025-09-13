import os
import gdown
import pandas as pd
import numpy as np
from astropy import units as u

import casparing.plot.CASPAR_util as cutil
color = cutil.AccDiagColor()

def CASPAR_loaddata(lit_database=None, caspar=None):
    if caspar is None:
        url = 'https://drive.google.com/uc?id=1QbJHcrndhaP2JIrBazy4zcUULpgClCdt76nMilKfRNs'
        output = "caspar.xlsx"
        gdown.download(url, output)
        
        df_caspar = pd.read_excel('caspar.xlsx', sheet_name='CASPAR', skiprows=[1])

    else:
        if '.csv' in caspar:
            df_caspar = pd.read_csv(caspar, skiprows=[1])
        elif '.xlsx' in caspar:
            df_caspar = pd.read_excel(caspar, sheet_name='CASPAR', skiprows=[1])
        else:
            raise ValueError('caspar must be a .csv or .xlsx with a sheet named CASPAR to load.')
        

        
    if lit_database is None:
        url = 'https://drive.google.com/uc?id=1QbJHcrndhaP2JIrBazy4zcUULpgClCdt76nMilKfRNs'
        output = "caspar.xlsx"
        gdown.download(url, output)
        
        df_lit = pd.read_excel('caspar.xlsx', sheet_name='Literature Database', skiprows=[1])

    else:
        if '.csv' in lit_database:
            df_lit = pd.read_csv(lit_database, skiprows=[1])
        elif '.xlsx' in lit_database:
            df_lit = pd.read_excel(lit_database, sheet_name='Literature Database', skiprows=[1])
        else:
            raise ValueError('lit_database must be a .csv or .xlsx with a sheet named Literature Database to load.')
        
    df_lit['Ha EW'] = df_lit['Ha EW'].astype(str).str.replace(",", "").astype(float)
    df_caspar['Ha EW'] = df_caspar['Ha EW'].astype(str).str.replace(",", "").astype(float)
    
    df_lit = df_lit.sort_values(by=['Unique ID'], ignore_index=True)
    df_caspar = df_caspar.sort_values(by=['Unique ID'], ignore_index=True)
    df_lit = df_lit.loc[(df_caspar['Upper Limit']!='LOW')]
    df_caspar = df_caspar.loc[(df_caspar['Upper Limit']!='LOW')]
    
    df_lit['log Mdot'] = np.log10(df_lit['Accretion Rate'])
    df_lit['log Mass'] = np.log10(df_lit['Mass'])
    df_lit['log Lum'] = np.log10(df_lit['Luminosity'])
    df_lit['log Age'] = np.log10((df_lit['Age'].values*u.Myr).to(u.yr).value) 
    df_lit['Sp Type Num'] = cutil.to_SpTyNum(df_lit['Sp Type'].values)

    df_caspar['log Mdot'] = np.log10(df_caspar['Accretion Rate'])
    df_caspar['log Mass'] = np.log10(df_caspar['Mass'])
    df_caspar['log Lum'] = np.log10(df_caspar['Luminosity'])
    df_caspar['log Age'] = np.log10((df_caspar['Age'].values*u.Myr).to(u.yr).value) 
    df_caspar['Sp Type Num'] = cutil.to_SpTyNum(df_caspar['Sp Type'].values)
    
    df_lit['AD'] = AccDiag_AssD(df_lit)[0]
    df_caspar['AD'] = AccDiag_AssD(df_caspar)[0]

    df_lit['Main Association'] = AccDiag_AssD(df_lit)[2]
    df_caspar['Main Association'] = AccDiag_AssD(df_caspar)[2]

    df_lit['color'] = AccDiag_AssD(df_lit)[3]
    df_caspar['color'] = AccDiag_AssD(df_caspar)[3]
    
    df_caspar['log Mass err'] = 0.434 * (df_caspar['Mass err']/df_caspar['Mass'])
    df_caspar['log Mdot err'] = 0.434 * (df_caspar['Accretion Rate err']/df_caspar['Accretion Rate'])    
    
    ageyr = (df_caspar['Age'].values*u.Myr).to(u.yr).value 
    ageerryr = (df_caspar['Age err'].values*u.Myr).to(u.yr).value 
    df_caspar['log Age err'] = 0.434 * ageerryr/ageyr
    df_caspar.loc[df_caspar['Age err']>df_caspar['Age'], 'log Age err'] = np.nanmean(0.434 * ageerryr/ageyr)

    df_lit['log Mass err'] = 0.434 * (df_lit['Mass err']/df_lit['Mass'])
    df_lit['log Mdot err'] = 0.434 * (df_lit['Accretion Rate err']/df_lit['Accretion Rate'])
    ageyr = ((df_lit['Age'].values)*u.Myr).to(u.yr).value
    ageerryr = ((df_lit['Age err'].values)*u.Myr).to(u.yr).value
    df_lit['log Age err'] = 0.434 * ageerryr/ageyr
    df_lit.loc[df_lit['Age err']>df_lit['Age'], 'log Age err'] = np.nanmean(0.434 * ageerryr/ageyr) 
    
    realupdate = df_caspar.dropna(subset=["log Mass", "log Mdot err"])
    medmdoterrupdate = np.nanmean(realupdate['log Mdot err'])
    medmasserrupdate = np.nanmean(realupdate['log Mass err'])

    realorig = df_lit.dropna(subset=["log Mass", "log Mdot err"])
    medmdoterrorig = np.nanmean(realorig['log Mdot err'])
    medmasserrorig = np.nanmean(realorig['log Mass err'])
    
    df_caspar.loc[df_caspar['Upper Limit']=='UPP', 'log Mdot err'] = medmdoterrupdate
    df_caspar.loc[df_caspar['log Mdot err'].isnull(), 'log Mdot err'] = medmdoterrupdate
    df_caspar.loc[df_caspar['log Mass err'].isnull(), 'log Mass err'] = medmasserrupdate
    
    df_lit.loc[df_lit['Upper Limit']=='UPP', 'log Mdot err'] = 0
    df_lit.loc[df_lit['log Mdot err'].isnull(), 'log Mdot err'] = medmdoterrorig
    df_lit.loc[df_lit['log Mass err'].isnull(), 'log Mass err'] = medmasserrorig

    df_caspar['Upper Limit bool'] = df_caspar['Upper Limit'] != 'UPP' # true is detected, false is upper limit
    df_lit['Upper Limit bool'] = df_lit['Upper Limit'] != 'UPP' # true is detected, false is upper limit
    
    planets = df_lit.loc[(df_lit['Companion']=='COM') & (df_lit['Mass']<=0.075)]
    planets = planets.astype(df_caspar.dtypes.to_dict())
    df_caspar.loc[planets.index, :] = planets[:]
    
    return df_lit, df_caspar

    
def AccDiag_AssD(df):
    AD = []
    COL = []
    for index, row in df.iterrows():
        if 'Line Luminosity' in row['Accretion Diagnostic']:
            if 'Ha photometry' in row['Tracer']:
                AD.append('Hα Photometric Luminosity')
                COL.append(color['Hα Photometric Luminosity'])
            else:
                AD.append('Line Luminosity')
                COL.append(color['Line Luminosity'])
        elif 'Continuum Excess' in row['Accretion Diagnostic']:
            AD.append('Continuum Excess')
            COL.append(color['Continuum Excess'])
        else:
            AD.append('Line Profile')
            COL.append(color['Line Profile'])
            
    Ass_d = []
    Ass_name = []
    for index, row in df.iterrows():
        if isinstance(row['Association'], str):
            if 'Sh 2-284' in row['Association']:
                Ass_d.append(1)
                Ass_name.append('Sh 2-284')
            elif 'Lagoon Nebula' in row['Association']:
                Ass_d.append(2)
                Ass_name.append('Lagoon Nebula')
            elif 'Centaurus' in row['Association']:
                Ass_d.append(3)
                Ass_name.append('Upper Centaurus Lupus')
            elif 'Tucana-Horologium' in row['Association']:
                Ass_d.append(4)
                Ass_name.append('Tucana-Horologium')
            elif 'IC 348' in row['Association']:
                Ass_d.append(5)
                Ass_name.append('IC 348')
            elif 'Taurus' in row['Association']:
                Ass_d.append(6)
                Ass_name.append('Taurus')
            elif '118 Tau' in row['Association']:
                Ass_d.append(7)
                Ass_name.append('118 Tau')
            elif 'σ Ori' in row['Association']:
                Ass_d.append(8)
                Ass_name.append('σ Ori')
            elif 'NGC 2024' in row['Association']:
                Ass_d.append(9)
                Ass_name.append('NGC 2024')
            elif 'Chamaeleon I' in row['Association']:
                Ass_d.append(10)
                Ass_name.append('Chamaeleon I')
            elif 'TW Hya' in row['Association']:
                Ass_d.append(11)
                Ass_name.append('TW Hya')
            elif 'Upper Scorpius' in row['Association']:
                Ass_d.append(12)
                Ass_name.append('Upper Scorpius')
            elif 'Lupus' in row['Association']:
                Ass_d.append(13)  
                Ass_name.append('Lupus')
            elif 'Argus' in row['Association']:
                Ass_d.append(14)
                Ass_name.append('Argus')
            elif 'p Oph' in row['Association']:
                Ass_d.append(15) 
                Ass_name.append('p Oph')
            elif 'Corona-Australis' in row['Association']:
                Ass_d.append(16) 
                Ass_name.append('Corona-Australis')
            elif 'n Chamaeleontis' in row['Association']:
                Ass_d.append(17) 
                Ass_name.append('n Chamaeleontis')
            elif '25 Orionis' in row['Association']:
                Ass_d.append(18) 
                Ass_name.append('25 Orionis')
            elif 'B Pictoris' in row['Association']:
                Ass_d.append(19) 
                Ass_name.append('B Pictoris')
            else:
                Ass_d.append(np.nan)
                Ass_name.append('Field')
        else:
            Ass_d.append(np.nan)
            Ass_name.append('Field')
    return AD, Ass_d, Ass_name, COL

    
def CASPAR_separateByMass(df_caspar):
    
    stars = df_caspar.loc[(df_caspar['Mass'] > 0.075) & (df_caspar['Companion'] != 'COM')]
    bds = df_caspar.loc[(df_caspar['Mass'] <= 0.075) & (df_caspar['Companion'] != 'COM')]
    planets = df_caspar.loc[(df_caspar['Companion']=='COM') & (df_caspar['Mass']<=0.075)]

    Supp = stars.loc[stars['Upper Limit'] == 'UPP']
    Sreg = stars.loc[stars['Upper Limit'] != 'UPP']

    BDupp = bds.loc[bds['Upper Limit'] == 'UPP']
    BDreg = bds.loc[bds['Upper Limit'] != 'UPP']

    Pupp = planets.loc[planets['Upper Limit'] == 'UPP']
    Preg = planets.loc[planets['Upper Limit'] != 'UPP']

    return stars, bds, planets, [Sreg, Supp], [BDreg, BDupp], [Preg, Pupp] 


def CASPAR_separateByUppLimit(df):

    REG = df.loc[df['Upper Limit'] != 'UPP']
    UPP = df.loc[df['Upper Limit'] == 'UPP']
   
    return REG, UPP


def CASPAR_separateByAge(df_caspar, agelow, ageupp):
    df_age = df_caspar.loc[(df_caspar['Age']>agelow) &
               (df_caspar['Age']<= ageupp)]
    return df_age