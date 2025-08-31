import os
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import scipy.interpolate as sinterp
import scipy.stats as st
from scipy.optimize import minimize

from astropy import units as u

from uncertainties.umath import *
from isochrones.interp import DFInterpolator

HERE = os.path.dirname(os.path.abspath(__file__))

baraffe_df = pd.read_csv(os.path.join(HERE, '../data/parameters_data/models/Baraffe_models.csv'), index_col=[0,1])

mesa_df = pd.read_csv(os.path.join(HERE, '../data/parameters_data/models/mist_models.csv'), index_col=[0,1])


class interpIsochrones:
    def __init__(self, df, age, mass=None, temp=None):
        ''' THESE FUNCTIONS ARE MODIFIED FROM ISOCHRONES PACKAGE TO INTERPOLATE OVER MASS OR TEMP
        USING MORE MODELS BESIDES MIST
        
        Purpose 
        ----------
        Interpolates over isochrones and finds the closest isochrone
        
        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame of objects in which to find isochrones
        age : astropy.units.Quantity, 
            age of object in Myr
        mass : astropy.units.Quantity, optional
            mass of object in Msun.  if temp is None, mass must be supplied 
        temp : astropy.units.Quantity, optional
            temperature of object in K.  if mass is None, temp must be supplied
    
        '''
        
        self.age = age
        self.mass = mass
        self.temp = temp

        if (self.mass is not None ) & (self.temp is not None):
            self.interp_over = 'Teff'
            self.VAL = temp
        elif (self.mass is None ) & (self.temp is not None):
            self.interp_over = 'Teff'
            self.VAL = temp
        elif (self.mass is not None ) & (self.temp is None):
            self.interp_over = 'mass'
            self.VAL = mass
        else:
            raise ValueError('mass or temp must be given')
        self.df = df


    def interp_value2(self, pars, valu):
        interp = DFInterpolator(self.df)
        return interp(pars, valu)

    def temp_age_resid2(self, eep):
        '''modified from isochrones/models.py'''
        VAL_interp = self.interp_value2([self.age, eep], [self.interp_over])
#        print((self.VAL - VAL_interp) ** 2   )
        return (self.VAL - VAL_interp) ** 2    

    def get_eep2(self, **kwargs):
        '''modified from isochrones/models.py'''
        grid = self.df
        if (isinstance(self.VAL, float) or isinstance(self.VAL, int)) and (isinstance(self.age, float) or isinstance(self.age, int)):
            return self.get_eep_accurate2(**kwargs)


    def get_eep_accurate2(self, eep0=15,
        resid_tol=90,
        method="nelder-mead",
        return_object=False,
        return_nan=True,
        **kwargs):
        '''
        modfied from isochrones/models.py
        
        '''

        eeps_to_try = [1, 10, 30, 100, 200]
        while np.isnan(self.temp_age_resid2(eep0)):
            
            try:
                eep0 = eeps_to_try.pop()
                
            except IndexError:
                if return_nan:
                    return np.nan
                else:
                    raise ValueError("eep0 gives nan for all initial guesses! {}".format((self.VAL, self.age)))

        result = minimize(self.temp_age_resid2, eep0, method=method, options=kwargs)

        if return_object:
            return result

        if result.success and result.fun < resid_tol ** 2:
#            print(float(result.x))
            return float(result.x)
        else:
            if return_nan:
                return np.nan
            else:
                raise RuntimeError("EEP minimization not successful: {}".format((self.VAL, self.age)))
    
class physical_params:
    def __init__(self, age, **kwargs):
        ''' 
        Purpose 
        ----------
        Calculate physical params using nearest isochrone and spectral type-temperature conversion
        
        Parameters
        ----------
        age : astropy.units.Quantity, 
            age of object in Myr
        
        kwargs: 
        - !one of these must be provided! 
            sptype : string
                spectral type of object in letter number format
            temp : astropy.units.Quantity
                temperature in K
            mass : astropy.units.Quantity
                mass in Msun
            value : astropy.units.Quantity,
                either spectral type, mass, or temperature.  The type of object will be self identified.  
        '''
        if isinstance(age, float):
            age = age * u.Myr
        if age <= 0.1 *u.Myr: # SINCE AGE ISN'T ~THAT~ ACCURATE, ROUND ANYTHING LESS THAT 0.5 MYR UP TO 0.5 MYR
            self.age = 0.1 * u.Myr
        else:
            self.age = age
 
        if 'sptype' in kwargs:
            self.method = 'sptype'
            self._sptype = kwargs['sptype']
            # self.specnum = self.to_SpTyNum()
            # self.Teff, self.ST_who_won = self.SpTy_to_Teff(which_model=True)
            # self.mass, self.modelBestTeff, self.radius, self.lum, self.logg, self.IC_who_won = self.TeffAge_to_Params(temp=self.Teff.value, which_model=True)
        elif 'mass' in kwargs:
            self.method = 'mass'
            self._mass = kwargs['mass']
            # _, self.Teff, self.radius, self.lum, self.logg, self.IC_who_won = self.TeffAge_to_Params(mass=mass.value, which_model=True) 
            # self.specnum, self.ST_who_won = self.Teff_to_SpTy(Teff=Teff, which_model=True)
            # self.sptype = self.Num_to_SpTy(specnum=specnum)
        elif 'temp' in kwargs:
            self.method = 'temp'
            self._Teff = kwargs['temp']
            # self.specnum, self.ST_who_won = self.Teff_to_SpTy(which_model=True)    
            # self.sptype = self.Num_to_SpTy(specnum=specnum)
            # self.mass, self.modelBestTeff, self.radius, self.lum, self.logg, self.IC_who_won = self.TeffAge_to_Params(temp=self.Teff.value, which_model=True) 

        elif 'value' in kwargs:
            if isinstance(kwargs['value'], str):
                self.method = 'sptype'
                self._sptype = kwargs['value']
                # self.specnum = self.to_SpTyNum()
                # self.Teff, self.ST_who_won = self.SpTy_to_Teff(which_model=True)
                # self.mass, self.modelBestTeff, self.radius, self.lum, self.logg, self.IC_who_won = self.TeffAge_to_Params(temp=self.Teff.value, which_model=True)
            elif kwargs['value'].unit == u.Msun:
                self.method = 'mass'
                self._mass = kwargs['value']
                # _, self.Teff, self.radius, self.lum, self.logg, self.IC_who_won = self.TeffAge_to_Params(mass=self.mass.value, which_model=True) 
                # self.specnum, self.ST_who_won = self.Teff_to_SpTy(Teff=self.Teff, which_model=True)
                # self.sptype = self.Num_to_SpTy(specnum=self.specnum)
            elif kwargs['value'].unit == u.K:
                self.method = 'temp'
                self._Teff = kwargs['value']
                # self.specnum, self.ST_who_won = self.Teff_to_SpTy(which_model=True)    
                # self.sptype = self.Num_to_SpTy(specnum=self.specnum)
                # self.mass, self.modelBestTeff, self.radius, self.lum, self.logg, self.IC_who_won = self.TeffAge_to_Params(temp=self.Teff.value, which_model=True) 
            else:
                raise ValueError('value is not a spectral type, mass, or temp.')

        # elif 'nophysicalparams' in kwargs:
        #     self.method=None
            
        else:
            raise ValueError('must have either sptype, temp, mass, or value')
    
    def get_params(self, which_models=False):
        ''' calculates all physical params '''
        if self.method == 'sptype':
            sptype = self._sptype
            specnum = self.to_SpTyNum()
            if which_models:
                Teff,ST_who_won = self.SpTy_to_Teff(which_model=which_models)
                mass, modelBestTeff, radius, lum, logg,IC_who_won = self.TeffAge_to_Params(temp=Teff.value, which_model=which_models) # Msun, Rsun, Lsun, logg
            else:
                Teff = self.SpTy_to_Teff(which_model=which_models)
                mass, modelBestTeff, radius, lum, logg = self.TeffAge_to_Params(temp=Teff.value, which_model=which_models)
        elif self.method == 'temp':
            if which_models:
                specnum, ST_who_won = self.Teff_to_SpTy(which_model=which_models) 
            else:
                specnum = self.Teff_to_SpTy(which_model=which_models) 
            sptype = self.Num_to_SpTy(specnum=specnum)
            Teff = self._Teff
            if which_models:
                mass, modelBestTeff, radius, lum, logg,IC_who_won = self.TeffAge_to_Params(temp=Teff.value, which_model=which_models) 
            else:
                mass, modelBestTeff, radius, lum, logg = self.TeffAge_to_Params(temp=Teff.value, which_model=which_models)
        
        elif self.method == 'mass':
            mass = self._mass
            if which_models:
                _, Teff, radius, lum, logg,IC_who_won = self.TeffAge_to_Params(mass=mass.value, which_model=which_models) 
                specnum,ST_who_won = self.Teff_to_SpTy(Teff=Teff,which_model=which_models)
            else:
                _, Teff, radius, lum, logg = self.TeffAge_to_Params(mass=mass.value, which_model=which_models) 
                specnum = self.Teff_to_SpTy(Teff=Teff,which_model=which_models)
            sptype = self.Num_to_SpTy(specnum=specnum)
            modelBestTeff=Teff
            
        else:
            raise ValueError('this should not have happened. how did you get his error?')
        if which_models:
            return sptype, specnum, Teff, mass, radius, lum, logg,ST_who_won, IC_who_won
        else:
            return sptype, specnum, Teff, mass, radius, lum, logg

    def which_models(self):
        ''' calculates which model is used (Herczeg/Pecaut, Baraffe/MIST)'''
        if self.method == 'sptype':
            Teff, ST_who_won = self.SpTy_to_Teff(which_model=True)
            _, modelBestTeff, _, _, _, IC_who_won  = self.TeffAge_to_Params(temp=Teff.value, which_model=True) # Msun, Rsun, Lsun, logg
        elif self.method == 'temp':
            specnum, ST_who_won = self.Teff_to_SpTy(which_model=True)  
            _, modelBestTeff, _, _, _, IC_who_won  = self.TeffAge_to_Params(temp=self.Teff.value, which_model=True) 
        elif self.method == 'mass':
            mass = self._mass
            _, Teff, _, _, _ , IC_who_won = self.TeffAge_to_Params(mass=mass.value, which_model=True) 
            specnum, ST_who_won = self.Teff_to_SpTy(Teff=Teff, which_model=True)
        return ST_who_won, IC_who_won
             
    def compare_temperature_results(self):
        '''compare given temperature to the best match modeled temperature'''
        if self.method == 'sptype':
            Teff = self.SpTy_to_Teff()
            _, modelBestTeff, _, _, _  = self.TeffAge_to_Params(temp=Teff.value) # Msun, Rsun, Lsun, logg
        elif self.method == 'temp':
            Teff = self._Teff
            _, modelBestTeff, _, _, _  = self.TeffAge_to_Params(temp=Teff.value) 
        elif self.method == 'mass':
            mass = self._mass
            _, Teff, _, _, _  = self.TeffAge_to_Params(mass=mass.value) 
            modelBestTeff = Teff
        return Teff, modelBestTeff

    @property
    def sptype(self):
        if self.method == 'sptype':
            sptype= self._sptype
        elif self.method == 'temp':
            specnum = self.Teff_to_SpTy()    
            sptype = self.Num_to_SpTy(specnum=specnum)
        else:
            _, Teff, _, _, _ = self.TeffAge_to_Params(mass=self.mass.value) 
            specnum = self.Teff_to_SpTy(Teff=Teff)
            sptype = self.Num_to_SpTy(specnum=specnum)
        return sptype
    @property     
    def specnum(self):
        if self.method == 'sptype':
            specnum = self.to_SpTyNum()
        elif self.method == 'temp':
            specnum = self.Teff_to_SpTy()  
        else:
            _, Teff, _, _, _ = self.TeffAge_to_Params(mass=self.mass.value) 
            specnum = self.Teff_to_SpTy(Teff=Teff)
        return specnum
    @property
    def Teff(self):
        if self.method == 'sptype':
            Teff = self.SpTy_to_Teff()
        elif self.method == 'temp':
            Teff = self._Teff
        else:
            _, Teff, _, _, _ = self.TeffAge_to_Params(mass=self.mass.value) 
        return Teff
    @property       
    def mass(self):
        if self.method == 'sptype':
            Teff = self.SpTy_to_Teff()
            mass, _, _, _, _ = self.TeffAge_to_Params(temp=Teff.value) #
        elif self.method == 'temp':
            mass, _, _, _, _ = self.TeffAge_to_Params(temp=self.Teff.value) 
        else:
            mass = self._mass
        return mass
    @property         
    def radius(self):
        if self.method == 'sptype':
            Teff = self.SpTy_to_Teff()
            _, _, radius, _, _ = self.TeffAge_to_Params(temp=Teff.value) #
        elif self.method == 'temp':
            _, _, radius, _, _ = self.TeffAge_to_Params(temp=self.Teff.value) 
        else:
            _, _, radius, _, _ = self.TeffAge_to_Params(mass=self.mass.value) 
        return radius
    @property 
    def luminosity(self):
        if self.method == 'sptype':
            Teff = self.SpTy_to_Teff()
            _, _, _, lum, _ = self.TeffAge_to_Params(temp=Teff.value) #
        elif self.method == 'temp':
            _, _, _, lum, _ = self.TeffAge_to_Params(temp=self.Teff.value) 
        else:
            _, _, _, lum, _ = self.TeffAge_to_Params(mass=self.mass.value) 
        return lum
             
    @property
    def logg(self):
        if self.method == 'sptype':
            Teff = self.SpTy_to_Teff()
            _, _, _, _, logg = self.TeffAge_to_Params(temp=Teff.value) #
        elif self.method == 'temp':
            _, _, _, _, logg = self.TeffAge_to_Params(temp=self.Teff.value) 
        else:
            _, _, _, _, logg = self.TeffAge_to_Params(mass=self.mass.value) 
        return logg
            
            
    def to_SpTyNum(self, **kwargs):
        '''
        This function takes a spectral type in the form - 'letternumber', i.e. 'M5'.
        It then translates that spec type into a spectral type number identifier used for interpolation.
        '''
        spty_dict = {'O' : 1,
                     'B' : 2,
                     'A' : 3,
                     'F' : 4,
                     'G' : 5,
                     'K' : 6,
                     'M' : 7,
                     'L' : 8}
        if self.method != 'sptype':
            if 'sptype' in kwargs:
                SPT = kwargs['sptype']
            else:
                raise ValueError('must include sptype')
        else:
            SPT = self.sptype
            
        if 'V' in SPT:
            SPT = SPT.split('V')[0]
        if '.' in SPT:
            dec = float('0.' + SPT.split('.')[-1])
            SPT = SPT.split('.')[0]
        else:
            dec = 0
        letter = SPT[0]
        number = SPT[1:]

        sptynum = spty_dict[letter]*10 + int(number) + dec

        return sptynum

    def Num_to_SpTy(self, **kwargs):
        '''
        This function takes a spectral type in the form - 'letternumber', i.e. 'M5'.
        It then translates that spec type into a spectral type number identifier used for interpolation.
        '''
        spty_dict = {1 : 'O',
                     2 : 'B',
                     3 : 'A',
                     4 : 'F',
                     5 : 'G',
                     6 : 'K',
                     7 : 'M',
                     8 : 'L'}
        
        if self.method != 'sptype':
            if 'specnum' in kwargs:
                SPN = kwargs['specnum']
            else:
                raise ValueError('must include sptype')
        else:
            SPN = self.to_SpTyNum()
        
        if int(SPN/10) > 8:
            letter = 'L'
            number = int(SPN%10)+10
        else:
            letter = spty_dict[int(SPN/10)]
            number = int(SPN%10)
        
        spty = letter+str(number)
        
        return spty
    
    def SpTy_to_Teff(self, which_model=False, **kwargs):
        '''
        This function will take a numerical spectral type identifier, and using interpolation from the tables in
        Herczeg and Hillenbrand 2014 it will calculate an effective temperature.
        '''
        if self.method != 'sptype':
            if 'specnum' in kwargs:
                SPN = kwargs['specnum']
            else:
                raise ValueError('must include specnum')
        else:
            if 'specnum' in kwargs:
                SPN = kwargs['specnum']
            else:
                SPN = self.to_SpTyNum()
        
        if (SPN < 45) | (SPN >= 80): 
            P = pd.read_csv(os.path.join(HERE,'../data/parameters_data/sptype_temp_conversion/Pecaut_SpTyTeff_Numericized.txt'), delimiter='\t', comment='#')
            SpTyNum = P['SpTyIdentifier']
            Teff= P['Teff_K']
            who_won = 'Pecaut'
        else:
            HH = pd.read_csv(os.path.join(HERE,'../data/parameters_data/sptype_temp_conversion/HerczegHillenbrand_SpTyTeff_Numericized.txt'), delimiter='\t')
            SpTyNum = HH['SpTyIdentifier']+20
            Teff= HH['Teff_K']
            who_won = "Herczeg & Hillenbrand 2014"
        
        spl = sinterp.UnivariateSpline(SpTyNum, Teff)

        teff = spl(SPN)
       
        if which_model:
          
            return teff * u.K, who_won
        else:
           
            return teff * u.K

    def Teff_to_SpTy(self, which_model=False, **kwargs):
        '''
        This function will take an effective temperature, and using interpolation from the tables in
        Herczeg and Hillenbrand 2014 it will calculate a numerical spectral type identifier.
        '''
        if self.method != 'temp':
            if 'Teff' in kwargs:
                temp = kwargs['Teff']
            else:
                raise ValueError('must include Teff')
        else:
            temp = self.Teff
        
        if (temp.value > 6600) | (temp.value < 2500): 
            P = pd.read_csv(os.path.join(HERE,'../data/parameters_data/sptype_temp_conversion/Pecaut_SpTyTeff_Numericized.txt'), delimiter='\t', comment='#')
            SpTyNum = P['SpTyIdentifier'].to_list()
            Teff= P['Teff_K'].to_list()
            who_won = 'Pecaut'
        else:
            HH = pd.read_csv(os.path.join(HERE,'../data/parameters_data/sptype_temp_conversion/HerczegHillenbrand_SpTyTeff_Numericized.txt'), delimiter='\t')
            SpTyNum = (HH['SpTyIdentifier']+20).to_list()
            Teff= HH['Teff_K'].to_list()
            who_won = "Herczeg & Hillenbrand 2014"
            
        SpTyNum.reverse()
        Teff.reverse()
        
        spl = sinterp.UnivariateSpline(Teff, SpTyNum)
        
        spty = spl(temp)
        if which_model:
            return spty, who_won
        else:
            return spty
    
    
    def TeffAge_to_Params(self, mass=None, temp=None, age=None, which_model=False):  
        '''
        convert temperature and age to physical params.  
        Must have either mass or temp be an argument!  
        '''
        if (mass is None ) & (temp is None):
            raise ValueError('mass or temp must be an argument of TeffAge_to_Params()')
        
        if age is None:
            age = self.age
        
        
        II = interpIsochrones(baraffe_df, age.value, mass=mass, temp=temp)
        eep = II.get_eep2()
        if np.isnan(eep):
#            print('---------------- trying MESA')
            age_mesa = np.log10(age.to(u.yr).value)
            
            II_mesa = interpIsochrones(mesa_df, age_mesa, mass=mass, temp=temp)
            eep_mesa = II_mesa.get_eep2()
            
            if np.isnan(eep_mesa):
                
                print('cannot get physical information.  Outside model parameters')
                print(mass, temp, age_mesa)
                if which_model:
                    return np.nan*u.Msun, np.nan*u.K, np.nan*u.Rsun, np.nan*u.Lsun, np.nan*u.cm/u.s**2, 'no one'
                else:
                    return np.nan*u.Msun, np.nan*u.K, np.nan*u.Rsun, np.nan*u.Lsun, np.nan*u.cm/u.s**2
            else:
                interp = DFInterpolator(mesa_df)
                mass, teff, radius, loglum, logg = interp([age_mesa, eep_mesa], ['mass', 'Teff', 'radius', 'logL', 'logg'])
                
                who_won = 'MIST MESA'
        else:
            interp = DFInterpolator(baraffe_df)
            
            mass,teff, radius, loglum, logg = interp([age.value, eep], ['mass', 'Teff', 'radius', 'logL', 'logg'])
            
            who_won = 'Baraffe 2015'

        
        lum = 10**loglum * u.Lsun
        mass = mass*u.Msun
        Teff = teff * u.K
        radius = radius * u.Rsun
        logg = logg * u.cm / u.s**2
        if which_model:
            return mass, Teff, radius, lum, logg, who_won
        else:
            return mass, Teff, radius, lum, logg
        
        
    def uncertainty(self, ageerr, SpTerr, **kwargs):
        ''' calculates uncertainties for all parameters.'''
        if 'autoparams' in kwargs:
            sptype, specnum, Teff, mass, radius, lum, logg = self.get_params()
        else:
            sptype, mass, radius, lum, logg = kwargs['sptype'], kwargs['mass'], kwargs['radius'], kwargs['lum'], kwargs['logg']
            specnum = self.to_SpTyNum(sptype=sptype)

        if isinstance(ageerr, float):
            ageerr = ageerr * u.Myr
        
        SPN = specnum
        # SPECTRAL TYPE UNCERAINTY
        if (SpTerr is None) or (np.isnan(SpTerr)) or (isinstance(SpTerr, u.Quantity)):
            if SPN > 71:
                SpT_err = 0.2
            elif (SPN <= 71) | (SPN >= 68): 
                SpT_err = 0.5
            else:
                SpT_err = 1

        else:
            SpT_err = SpTerr
            
        # TEMPERATURE UNCERTAINTY
        # 1. get specnum upper and lower
       
        spnum_upp = SPN - SpT_err
        spnum_lower = SPN + SpT_err
        Teff_lower = self.SpTy_to_Teff(specnum=spnum_lower)
        Teff_upp = self.SpTy_to_Teff(specnum=spnum_upp)
        Teff_err_low = self.Teff-Teff_lower
        Teff_err_upp = Teff_upp -self.Teff
        Teff_err = np.nanmean([Teff_err_low.value, Teff_err_upp.value]) * u.K

        # MASS, RADIUS, LUMINOSITY Uncertainty 
        ## 1. low age/low temp, 2. low age/high temp, 3. high age/low temp, 4. high age/high temp
        
        
        if self.age-ageerr > 0:
            
            # 1. low age/low temp
            masserr1, _, radiuserr1, lumerr1, loggerr1 = self.TeffAge_to_Params(temp=Teff_lower.value, age=self.age-(ageerr))

            # 2. low age/high temp
            #print(Teff_upp, self.age-(ageerr*u.Myr))
            masserr2, _, radiuserr2, lumerr2, loggerr2 = self.TeffAge_to_Params(temp=Teff_upp.value, age=self.age-(ageerr))
        else:
            masserr1, _, radiuserr1, lumerr1, loggerr1 =np.nan, np.nan, np.nan, np.nan,np.nan
            masserr2, _, radiuserr2, lumerr2, loggerr2 = np.nan, np.nan, np.nan, np.nan,np.nan
        
        # 3. high age/low temp
        #print(Teff_lower, self.age+(ageerr*u.Myr))
        masserr3, _, radiuserr3, lumerr3, loggerr3 = self.TeffAge_to_Params(temp=Teff_lower.value, age=self.age+(ageerr))
        
        # 4. high age/high temp
        #print(Teff_upp, self.age+(ageerr*u.Myr))
        masserr4, _, radiuserr4, lumerr4, loggerr4= self.TeffAge_to_Params(temp=Teff_upp.value, age=self.age+(ageerr))
        
        
        mass_err = np.nanmean([abs(mass-masserr1).value, abs(mass-masserr2).value, abs(mass-masserr3).value, abs(mass-masserr4).value])
        
        rad_err = np.nanmean([abs(radius-radiuserr1).value, abs(radius-radiuserr2).value, abs(radius-radiuserr3).value, abs(radius-radiuserr4).value])
        
        lum_err = np.nanmean([abs(lum-lumerr1).value, abs(lum-lumerr2).value, abs(lum-lumerr3).value, abs(lum-lumerr4).value])
        
        logg_err = np.nanmean([abs(logg-loggerr1).value, abs(logg-loggerr2).value, abs(logg-loggerr3).value, abs(logg-loggerr4).value])
        
        #print(Teff_err)
        #print(mass_err*u.Msun)
        #print(rad_err*u.Rsun)
        #print(lum_err*u.Lsun)
        #print(logg_err)
 
        return mass_err*u.Msun, Teff_err, rad_err*u.Rsun, lum_err*u.Lsun
    
    