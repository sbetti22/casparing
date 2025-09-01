import os
import numpy as np
import pandas as pd


def add_astrometry_objects(object_names, colname=None):
    if isinstance(object_names, str):
        if ('.csv' in object_names) or ('.txt' in object_names) or ('.dat' in object_names):
            df = pd.read_csv(object_names)
            if colname is not None:
                df = df[colname].rename({colname:'obj'})
        else:
            df = pd.DataFrame({'obj':[object_names]})
    elif isinstance(object_names, (list, np.ndarray)):
        object_names = np.array(object_names)
        if object_names.ndim > 1:
            raise ValueError(f'if obj is a list or nd.array, it must be 1-d.  obj is current {obj.ndim}-d')
        df = pd.DataFrame({'obj':object_names})
    elif isinstance(obj, pd.DataFrame):
        df = object_names
        if colname is not None:
            df = df[colname].rename({colname:'obj'})
    else:
        raise ValueError(f'obj must either be: 1D text file with extension .csv, .txt. or .dat, a 1D list or numpy.ndarray, or a 1D pandas dataframe')
    return df
        
def add_physparams_objects(object_names, age=None, age_err=None, value_id=None, value=None, value_err=None, verbose=True):
    if isinstance(object_names, str):
        if ('.csv' in object_names) or ('.txt' in object_names) or ('.dat' in object_names):
            if verbose: print(object_names + ' path provided.  Opening file as pd.DataFrame.')
            df = pd.read_csv(object_names)
            cols = list(df.columns)
            
            if (('Reference Name' in cols) and ('Age' in cols) and ('Age err' in cols) and ( (('Sp Type' in cols) and ('Sp Type err' in cols)) or (('Mass' in cols) and ('Mass err' in cols))  or (('Teff' in cols) and ('Teff err' in cols)) )):  
                if verbose: print('All required columns are provided.')
                return df
        
            if verbose: print('All required columns not in dataframe.  checking if they are provided as variables.')
            msg = ''
            if ((age is None) and ('Age' not in cols)):
                msg += '"Age" must be a column in dataframe or provided as a list as "age" variable. '
            if ((age_err is None) and ('Age err' not in cols)):
                msg += '"Age err" must be a column in dataframe or provided as a list as "age_err" variable. '
            if ((value_id is None) and (value_id not in cols)):
                msg += '"Sp Type", "Mass" or "Teff" must be a column in dataframe or provided as value_id.'
            if  (value is None):
                msg += '"Sp Type", "Mass" or "Teff" must be a column in dataframe or provided as a list as "value" variable. '
            if (value_err is None):
                msg += '"Sp Type err", "Mass err" or "Teff err" must be a column in dataframe or provided as a list as "value_err" variable. '
            if msg != '':    
                raise ValueError(msg)
        else:
            
            object_names = [object_names]
            
    if isinstance(object_names, pd.DataFrame):
        if verbose: print('DataFrame added.')
        df = object_names
        cols = list(df.columns)
        if (('Reference Name' in cols) and ('Age' in cols) and ('Age err' in cols) and ( (('Sp Type' in cols) and ('Sp Type err' in cols)) or (('Mass' in cols) and ('Mass err' in cols))  or (('Teff' in cols) and ('Teff err' in cols)) )):  
            if verbose: 
                print('All required columns are provided.')
                display(df)
            return df
        else:
            if verbose: print('All required columns not in dataframe.  checking if they are provided as variables.')
            msg = ''
            if ((age is None) and ('Age' not in cols)):
                msg += '"Age" must be a column in dataframe or values provided as a list as "age" variable. '
            if ((age_err is None) and ('Age err' not in cols)):
                msg += '"Age err" must be a column in dataframe or values provided as a list as "age_err" variable. '
            if (value_id is None) and (('Mass' not in cols) or ('Teff' not in cols) or ('Sp Type' not in cols)):
                msg += '"Sp Type", "Mass" or "Teff" must be a column in dataframe or provided as value_id. '
                
            if (value is None) and ((('Mass' not in cols) or ('Teff' not in cols) or ('Sp Type' not in cols))):
                msg += '"Sp Type", "Mass" or "Teff" must be a column in dataframe or values provided as a list as "value" variable. '
                
            if (value_err is None) and ((('Mass err' not in cols) or ('Teff err' not in cols) or ('Sp Type err' not in cols))):
                msg += '"Sp Type err", "Mass err" or "Teff err" must be a column in dataframe or values provided as a list as "value_err" variable. '
            if msg != '':    
                raise ValueError(msg)
     
    if isinstance(object_names, list):
        if verbose: print('Object Names provided.  New dataframe being created.')
        msg = ''
        if (age is None):
            msg += '"age" must be a provided. '
        if (age_err is None):
            msg += '"age_err" must be provided. '
        if (value_id is None):
            msg += 'value_id must be "Sp Type", "Mass" or "Teff".'
        if  (value is None):
            msg += 'value must be provided. '
        if (value_err is None):
            msg += 'value_err must be provided. '
        if msg != '':    
            raise ValueError(msg)
        df = pd.DataFrame()

    cols = list(df.columns)
    if 'Reference Name' not in cols:
        if isinstance(object_names, list):
            if verbose: print('Reference Name added to dataframe')
            df['Reference Name'] = object_names
    else:
        if verbose: print('Reference Name already in dataframe')
        
    if 'Age' not in cols:
        if verbose: print('Age added to dataframe')
        if isinstance(age, (float, int)):
            age = [age]
        df['Age'] = age
    else:
        if verbose: print('Age already in dataframe.')

    if 'Age err' not in cols:
        if verbose: print('Age err added to dataframe')
        if isinstance(age_err, (float, int)):
            age_err = [age_err]
        df['Age err'] = age_err
    else:
        if verbose: print('Age err already in dataframe')

    if (value_id == 'sptype') or (value_id == 'SpType') or (value_id == 'Sp Type'):
        value_id = 'Sp Type'
    elif (value_id == 'mass') or (value_id == 'Mass'):
        value_id = 'Mass'
    elif (value_id == 'teff') or (value_id == 'temp') or (value_id == 'temperature') or (value_id == 'Teff'):
        value_id = 'Teff'
    else:
        raise ValueError(f'{value_id} not a valid option. Choices are "sptype", "SpType" and "Sp Type", "mass" and "Mass", or "teff", "temp", "temperature" and "Teff"')
        
    if value_id not in cols:
        if verbose: print(f'{value_id} added to dataframe')
        if isinstance(value, (float, int, str)):
            value = [value]
        df[value_id] = value
    else:
        if verbose: print(f'{value_id} already in dataframe')
        
    if value_id + ' err' not in cols:
        if verbose: print(f'{value_id} err added to dataframe')
        if isinstance(value_err, (float, int, str)):
            value_err = [value_err]
        df[value_id + ' err'] = value_err
    else:
        if verbose: print(f'{value_id} err already in dataframe')
    return df
    

def save_astrometry_df(df, savename):
    columns = ['Simbad-Resolvable Name', 'Reference Name', 'RA (J2000.0)',
       'Dec (J2000.0)', 'RA (J2016.0)', 'Dec (J2016.0)', 'Association', 'Association Probability Banyan Sigma', 
               'Association Census Reference', 'Association Age', 'Association Age err', 'GAIA DR2 Source ID',
       'GAIA DR2 Parallax', 'GAIA DR2 Parallax err',
       'GAIA DR2 Reliable Parallax', 'GAIA DR2 Distance',
       'GAIA DR2 Distance lower limit', 'GAIA DR2 Distance upper limit',
       'GAIA DR2 RA proper motion', 'GAIA DR2 RA proper motion err',
       'GAIA DR2 Dec proper motion', 'GAIA DR2 Dec proper motion err',
       'GAIA EDR3 Source ID', 'GAIA EDR3 Parallax', 'GAIA EDR3 Parallax err',
       'GAIA EDR3 Reliable Parallax', 'GAIA EDR3 Geometric Distance',
       'GAIA EDR3 Geometric Distance lower limit', 'GAIA EDR3 Geometric Distance upper limit',
       'GAIA EDR3 RA proper motion', 'GAIA EDR3 RA proper motion err',
       'GAIA EDR3 Dec proper motion', 'GAIA EDR3 Dec proper motion err',
       'radial velocity', 'radial velocity err', 'Jmag', 'Jmag err',
       'Hmag', 'Hmag err', 'Kmag', 'Kmag err', 'Links']
    
    df[['GAIA DR2 Source ID', 'GAIA EDR3 Source ID']] = df[['GAIA DR2 Source ID', 'GAIA EDR3 Source ID']].astype(str)
    df.to_csv(savename, columns=columns, index=False)
    
def put_in_caspar_order(df):
    HERE = os.path.dirname(os.path.abspath(__file__))
    fil = os.path.join(HERE, '../data/accretion_data/caspar_columns.txt')

    cc = pd.read_csv(fil, header=None)
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
    
    
    
    
    
    
