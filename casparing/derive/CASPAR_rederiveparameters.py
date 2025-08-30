import pandas as pd
import numpy as np
from astropy import units as u
import .rederiveparameters as cpm

from pandarallel import pandarallel
pandarallel.initialize(progress_bar=True)

def _phys_quants(row):
    params =  cpm.physical_params(row['FIN AGE'], value=row['VAL'])

    sptype, specnum, Teff, mass, radius, lum, logg, STwon, ICwon = params.get_params(which_models=True)
    masserr, tefferr, raderr, lumerr = params.uncertainty(row['FIN AGE ERR'], row['VAL ERR'], sptype=sptype, mass=mass, radius=radius, lum=lum, logg=logg)

    row['Mass'] = mass.value
    row['Mass err'] = masserr.value
    row['Radius'] = radius.value
    row['Radius err'] = raderr.value
    row['Teff'] = Teff.value
    row['Teff err'] = tefferr.value
    row['Luminosity'] = lum.value
    row['Luminosity Err'] = lumerr.value
    row['log g'] = logg.value
    row['SpTemp Conversion Reference'] = STwon
    row['Evolutionary Models Reference'] = ICwon
    return row
    
    
def create_table(object_name, age, age_err, value_id, value, value_err):
    if isinstance(object_name, str):
        object_name = [object_name]

    if isinstance(age, (float, int)):
        age = [age]

    if isinstance(age_err, (float, int)):
        age_err = [age_err]

    if isinstance(value, (float, int, str)):
        value = [value]

    if (value_id == 'sptype') or (value_id == 'SpType'):
        value_id = 'Sp Type'

    if (value_id == 'mass'):
        value_id = Mass
    if (value_id == 'teff') or (value_id == 'temp') or (value_id == 'temperature'):
        value_id = 'Teff'

    if isinstance(value_err, (float, int, str)):
        value_err = [value_err]

    d = {'Reference Name':object_name, 'Age':age, 'Age err':age_err, value_id:value, value_id + ' err':value_err}
    df = pd.DataFrame(d)
    return df

class get_physical_params:
    def __init__(self, df=None):
        if df is None:
            self.df = pd.DataFrame()
        else:
            self.df = df
        

    def _extractParams(self, p1, p2, p3, val):
        conditions1 = []
        choices1= []
        conditions2 = []
        choices2= []
        if (p1 in self.df.columns) and (p1 + ' err' in self.df.columns):
            conditions1.append(self.df[p1].notna())
            choices1.append(self.df[p1])
            conditions2.append(self.df[p1 + ' err'].notna())
            choices2.append(self.df[p1 + ' err'])
        if (p2 in self.df.columns) and (p2 + ' err' in self.df.columns):
            conditions1.append(self.df[p2].notna())
            choices1.append(self.df[p2])
            conditions2.append(self.df[p2 + ' err'].notna())
            choices2.append(self.df[p2 + ' err'])
        if (p3 in self.df.columns) and (p3 + ' err' in self.df.columns):
            conditions1.append(self.df[p3].notna())
            choices1.append(self.df[p3])
            conditions2.append(self.df[p3 + ' err'].notna())
            choices2.append(self.df[p3 + ' err'])
            
        if len(conditions1) == 0:
            raise ValueError(f'Must have either {p1}, {p2}, or {p3} as a value or table column. Also check that you have corresonding error value or column as well: {p1} err, {p2} err, or {p3} err.')
        
        self.df[val] = np.select(conditions1, choices1, default=-1000)
        self.df[val + ' ERR'] = np.select(conditions2, choices2, default=-1000)
        
        
    def extractSptypeTempOrMass(self):
        self._extractParams('Sp Type', 'Teff', 'Mass', 'VAL')

    def extractAge(self):   
        self._extractParams('Age', 'Individual Age', 'Association Age', 'FIN AGE')    


    def params(self):
        self.df = self.df.parallel_apply(_phys_quants, axis=1, result_type='expand')
        self.df.drop(['VAL', 'VAL ERR', 'FIN AGE', 'FIN AGE ERR'], axis=1, inplace=True)
    
    @property
    def database(self):
        return self.df
    
    

