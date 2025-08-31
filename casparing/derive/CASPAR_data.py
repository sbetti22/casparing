


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
        
def add_accretion_objects(object_names, age=None, age_err=None, value_id=None, value=None, value_err=None):
    if isinstance(object_names, str):
        if ('.csv' in object_names) or ('.txt' in object_names) or ('.dat' in object_names):
            df = pd.read_csv(object_names)
            cols = list(df.columns)
            if ('Reference Name' not in cols) or ('Age' not in cols) or ('Age err' not in cols) or ('Sp Type' not in cols) or ('Mass' not in cols) or ('Teff' not in cols) or ('Mass err' not in cols) or ('Sp Type err' not in cols) or ('Teff err' not in cols):
                raise ValueError('column names must be: Reference Name, Age, Age err, and one of [Sp Type, Mass, Teff], and [Sp Type err, Mass err, Teff err].')
            return df
            
        else:
            object_names = [object_names]
            
    if isinstance(object_names, pd.DataFrame):
        df = object_names
        cols = list(df.columns)
        if ('Reference Name' not in cols) or ('Age' not in cols) or ('Age err' not in cols) or ('Sp Type' not in cols) or ('Mass' not in cols) or ('Teff' not in cols) or ('Mass err' not in cols) or ('Sp Type err' not in cols) or ('Teff err' not in cols):
            raise ValueError('column names must be: Reference Name, Age, Age err, and one of [Sp Type, Mass, Teff], and [Sp Type err, Mass err, Teff err].')
        
    if (age is None) or (age_err is None) or (value_id is None) or (value is None) or (value_err is None):
        raise ValueError('if object_names is not filepath or dataframe, age, age_err, value_id, value, and value_err must be provided' )
            
    if isinstance(age, (float, int)):
        age = [age]

    if isinstance(age_err, (float, int)):
        age_err = [age_err]

    if isinstance(value, (float, int, str)):
        value = [value]

    if (value_id == 'sptype') or (value_id == 'SpType'):
        value_id = 'Sp Type'
    if (value_id == 'mass'):
        value_id = 'Mass'
    if (value_id == 'teff') or (value_id == 'temp') or (value_id == 'temperature'):
        value_id = 'Teff'

    if isinstance(value_err, (float, int, str)):
        value_err = [value_err]

    d = {'Reference Name':object_name, 'Age':age, 'Age err':age_err, value_id:value, value_id + ' err':value_err}
    df = pd.DataFrame(d)
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
    df.to_csv(savename, columns=columns)
    
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
    
    
    
    
    
    
