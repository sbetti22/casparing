### written by Sarah Betti 2022
### functions to derive physical and accretion properties for literature values 
### used primarily for CASPAR 
### use with FAST_get_physical_pararms

# must have 
# 1. alcala2017_linear_fits.csv
# 2. StellarParams/Pecaut_SpTyTeff_Numericized.txt
# 3. StellarParams/HerczegHillenbrand_SpTyTeff_Numericized.txt


#interpIsochrones class is heavily modified from isochrones package.  All parts that are not mine are taken from that package and property of Morton2015: 
#https://isochrones.readthedocs.io/en/latest/
#https://ui.adsabs.harvard.edu/abs/2015ascl.soft03010M/abstract (Morton2015)

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import scipy.stats as st

from astropy import constants as const
from astropy import units as u

from uncertainties import ufloat
from uncertainties.umath import *

import extinction as ex

class rederive_Mdot:
    def __init__(self, mass, radius, distance):
        '''

        age, sptype, specnum, Teff, mass, radius, lum, logg, distance, 

        Purpose 
        ----------
        Calculate accretion properties (line flux, line luminosity, accretino luminosity, accretion rate)
        
        Parameters
        ----------
        distance : astropy.units.Quantity, 
            distance of object in pc
        radius : astropy.units.Quantity
            radius in Rsun
        mass : astropy.units.Quantity
            mass in Msun

        '''

        self.mass = mass
        self.distance = distance
        self.radius = radius
        
        if isinstance(self.mass, (int, float)):
            self.mass = self.mass * u.Msun
        if isinstance(self.distance, (int, float)):
            self.distance = self.distance * u.pc
        if isinstance(self.radius, (int, float)):
            self.radius = self.radius * u.Rsun
    


    def emission_lines(self, line):
        '''return wavelenth in Angstrom of line'''
        emission_lines = {'Ha':6563, 'Hb':4861, 'Hg':4341, 'Hd':4102, 'He':3970 , 'H11':3371, 
          'PaB':1282, 'PaG':1094, 'PaD':1005, 'BrG':2165, 'PfB':4650, 
          'HeI587':5876, 'HeI667':6678, 'HeI706':7065, 'CaIIK':3934, 'CaIIH':3969, 
         'CaII866':8662, 'Li':6708, 'OI':8446}
        if line in emission_lines.keys():
            return emission_lines[line]
        else:
            raise ValueError(f'line does not exist.  must use: {list(emission_lines.keys())}')


    def lineflux_to_Lacc(self, flux, line, Av=None, A=None, B=None, wave=None, **kwargs):
        '''
        This function will turn a line flux into an accretion luminosity estimate using the Lacc-Lline fits derived
        by Alcala et al 2017. 
        
        Inputs:
        flux: astropy.Quantity
            line flux in erg/(s*cm^2)
        line: string
            name of emission line
        A : float
            If you want to input the parameters for your own line flux vs Lacc relationship
        B : float
            If you want to input the parameters for your own line flux vs Lacc relationship
        Av : float
            extinction 
        wave : float
            wavelength of line in A
        
        kwargs:
            uncertainty : boolean
                True to calculate uncertainty.  If True, below quantites must be provided
            fluxerr : astropy.Quantity
                error on line flux in erg/s/cm^2
            disterr : astropy.Quantity
                error on distance in pc 
                if not given, provide lower and upper bounds
            disterrlow : astropy.Quantity
                lower limit on distance in pc 
            disterrupp : astropy.Quantity
                upper limit on distance in pc 
        '''
        if isinstance(flux, float):
            flux = flux * u.erg/u.s/u.cm**2.
        
        alcala_lines = pd.read_csv('../data/accretion_data/alcala2017_linear_fits.csv', comment='#')

        if line is None:
            a = A
            b = B
            wave=None
            if a is None:
                raise ValueError('must provide A and B')
        elif line in alcala_lines.Diagnostic.values: 
            a = alcala_lines['a'].loc[alcala_lines['Diagnostic']==line].values[0]
            b = alcala_lines['b'].loc[alcala_lines['Diagnostic']==line].values[0]
            a_err = alcala_lines['a_err'].loc[alcala_lines['Diagnostic']==line].values[0]
            b_err = alcala_lines['b_err'].loc[alcala_lines['Diagnostic']==line].values[0]
            wave = alcala_lines['wavelength'].loc[alcala_lines['Diagnostic']==line].values[0]
        else:            
            raise ValueError(f'Line not found: please use either None or one of the following lines: {alcala_lines.Diagnostic.to_list()}')

        if Av is not None:
            #extinction correction
            if not np.isnan(Av):
                deredflux = ex.remove(ex.ccm89(np.array([wave*10]), Av, 3.1), flux.value)[0] * flux.unit
            else:
                deredflux = flux
            
        else:
            deredflux = flux
            
        ## find Lline in erg/s
        Lline = deredflux * (4*np.pi*self.distance**2)
        #convert to solar luminosity
        Lline = Lline.to(u.Lsun)
        #Find Lacc using Alcala relationships
        logLacc = a*np.log10(Lline.value)+b
        #solar luminosity
        Lacc = 10**(logLacc) * u.Lsun
        
        if kwargs.get('uncertainty'):
            if not np.isnan(kwargs['fluxerr']):
                if isinstance(kwargs['fluxerr'], float):
                    fluxerr = kwargs['fluxerr'] * u.erg/u.s/u.cm**2.
                else:
                    fluxerr = kwargs['fluxerr']
                
                if Av is not None:
                    #extinction correction
                    if not np.isnan(Av):
                        deredfluxerr = ex.remove(ex.ccm89(np.array([wave*10]), Av, 3.1), fluxerr.value)[0] * fluxerr.unit
                    else:
                        deredfluxerr = fluxerr
                else:
                    deredfluxerr = fluxerr
            else:
                fluxerr = 0 * u.erg/u.s/u.cm**2.
                deredfluxerr = fluxerr
                
            if 'disterr' in kwargs:
                disterr = kwargs['disterr']
            elif 'disterrlow' in kwargs:
                disterrlow = abs(self.distance - kwargs['disterrlow'])
                disterrupp = abs(self.distance - kwargs['disterrupp'])
            else:
                raise ValueError('must have either disterr, or disterrlow/disterrupp')

            Lacc_Lsun_err =  []
            for disterr in [disterrlow, disterrupp]:
                sigmaLline = Lline * np.sqrt((deredfluxerr/deredflux)**2. + ((2*disterr)/self.distance)**2.)  # Lsun

                sigmaLogLline = 0.434 * sigmaLline/Lline

                log10L = np.log10(Lline.value)
                sigma_log10L = 0.434 * sigmaLline.value/Lline.value

                alog10L = a * log10L
                alog10Lb = a * log10L + b 

                logsigma_Lacc = np.sqrt( b_err**2. + (alog10L * np.sqrt( (a_err/a)**2. + (sigma_log10L/log10L)**2. ))**2.) 

                sigma_Lacc = 2.303 * (10**(alog10Lb)) * logsigma_Lacc
                Lacc_Lsun_err.append(sigma_Lacc)
            
            return Lacc, np.nanmean(Lacc_Lsun_err)* u.Lsun
            
        else:
            return Lacc 
    
    

    
    def Lacc_to_lineflux(self, logL_acc_orig, dist_orig, A, B, **kwargs):
        ''' 
        if a paper only gives EW and logLacc, go backwards from logL to line luminosity using the paper"s scaling relation and then to line flux using paper's distance
        
        Inputs:
        logL_acc_orig: astropy.Quantity
            original log accretion luminosity in Lsun
        dist_orig: astropy.Quantity
            original distance in pc
        A : float
            input the parameters for your own line flux vs Lacc relationship
        B : float
            input the parameters for your own line flux vs Lacc relationship
            
        kwargs:
            uncertainty : boolean
                True to calculate uncertainty.  If True, below quantites must be provided
            Aerr: float
                error on line flux vs Lacc relationship slope
            Berr:
                error on line flux vs Lacc relationship intercept

        '''

        
        # step 1. go from log L_acc to line luminosity
        log_line_lum = ((logL_acc_orig.value - B) / A) * logL_acc_orig.unit
        line_lum = 10 ** (log_line_lum.value) * log_line_lum.unit
        
        # step 2: go from line_lum to line flux
        line_flux = line_lum / (4 * np.pi * dist_orig**2.)
        
        if kwargs.get('uncertainty'):
            if 'Aerr' in kwargs:
                Aerr=kwargs['Aerr']
            if 'Berr' in kwargs:
                Berr = kwargs['Berr']
                
            ## NEED TO FIGURE OUT HOW TO FINISH!!
            
        return line_flux.to(u.erg/u.s/u.cm**2.)
    
    def Mdot_to_lineflux(self, Mdot_orig, mass_orig, radius_orig, dist_orig, A, B, ignore_Mdot_constants=False, Rin = 5* u.Rsun, **kwargs):
        ''' 
        if a paper only gives EW and Mdot, go backwards from Mdot to logL to line luminosity using the paper"s scaling relation and then to line flux using paper's distance
        
        Inputs:
        Mdot_orig : astropy.Quantity
            original log accretion luminosity in Lsun
        mass_orig : astropy.Quantity
            original mass in Msun
        radius_orig : astropy.Quantity
            original radius in Rsun
        dist_orig : astropy.Quantity
            original distance in pc
        A : float
            input the parameters for your own line flux vs Lacc relationship
        B : float
            input the parameters for your own line flux vs Lacc relationship
        ignore_Mdot_constants : boolean
            if True, ignore the 1-R/Rin in calculating Lacc (default False)
        Rin : astropy.Quantity
            truncation radius (default = 5 Rsun)
            
        kwargs:
            uncertainty : boolean
                True to calculate uncertainty.  If True, below quantites must be provided
            Aerr: float
                error on line flux vs Lacc relationship slope
            Berr:
                error on line flux vs Lacc relationship intercept

        '''
        
        # step 1: go from Mdot to Lacc ()
        if ignore_Mdot_constants:
            Lacc = (Mdot_orig * const.G * mass_orig) / radius_orig
        else:
            Lacc = ((1-(self.radius/Rin)) * Mdot_orig * const.G * mass_orig) / radius_orig
        Lacc = Lacc.to(u.Lsun)
        log_Lacc = np.log10(Lacc.value) 
        
        # step 2: go from Lacc to line luminosity
        log_line_lum = (log_Lacc - B) / A
        line_lum = 10 ** (log_line_lum)* Lacc.unit
        
        # step 3: go from line_lum to line flux
        line_flux = line_lum / (4 * np.pi * dist_orig**2.)
        line_flux = line_flux.to(u.erg/u.s/u.cm**2)
        if kwargs.get('uncertainty'):
            if 'Aerr' in kwargs:
                Aerr=kwargs['Aerr']
            if 'Berr' in kwargs:
                Berr = kwargs['Berr']
                
                
        return log_Lacc, line_flux
            
    def rescale_Lacc(self, logL_acc_orig, dist_orig, Rin=5*u.Rsun):
        '''
        OLD --- This function will scale a log accretion luminosity in Lsun by a new distance
        and convert it into Mdot.  This is ONLY for 
        Excess Balmer accretion luminosities or if Lacc originally calculated from Alcala+2017 scaling relations.  We scale the accretion luminosity by the distance.  
        
        Inputs:
        logL_acc_orig:  astropy.Quantity 
            original accretion luminosity [Lsun]
        dist_orig : astropy.Quantity  
            original distance to object [pc]
        Rin: astropy.Quantity, optional
            magnetospheric radius (default 5 [Rsun])

        '''
        if isinstance(logL_acc_orig, float):
            logL_acc_orig = logL_acc_orig * u.Lsun
        if isinstance(dist_orig, float):

            dist_orig = dist_orig * u.pc
        elif np.isnan(dist_orig):
            dist_orig = self.distance

        L_acc_orig = ((10**(logL_acc_orig.value)) * logL_acc_orig.unit)#.to(u.erg/u.s)
        L_acc_new = L_acc_orig * (self.distance**2. / dist_orig**2.)
        
        return L_acc_new.to(u.Lsun)

    def photHaEW_to_Lacc(self, EW, Hamag,  returnFline=False, **kwargs):
        '''
        Calculate Lacc from Ha EW from Kalari 2015 paper ONLY
        
        Inputs:
        EW:  astropy.Quantity 
            original Ha equlivalent width [A]
        Hamag : float  
            Ha magnitude 
        returnFline: boolean, optional
            if True, return the line flux and luminosity as well as the accretion luminosity
            
        kwargs:
            uncertainty : boolean
                True to calculate uncertainty.  If True, below quantites must be provided
            disterr : astropy.Quantity
                error on distance in pc 
                if not given, provide lower and upper bounds
            disterrlow : astropy.Quantity
                lower limit on distance in pc 
            disterrupp : astropy.Quantity
                upper limit on distance in pc  
            Hamagerr : float
                error on Ha magnitude
            EWerr : float
                error on EW in A

        '''
        W = 103.76 * u.AA
        Ftotal = (1.84e-7 * u.erg/u.s/u.cm**2.) * (10**(-0.4*(Hamag + 0.03)))
        F_Ha = Ftotal * ( (-EW/W)/ (1-(EW/W)))
        
        L_Ha = F_Ha * 4*np.pi * self.distance**2.
        L_Ha = L_Ha.to(u.Lsun)
        
        # use Alcala 2017 conversion instead of Barentsen 2011 (A=1.13, B=1.93 ) like in paper (uses different intercept, 
        # but slope is the same) 
        
        # T = mass_to_Teff(self.mass.value, self.age.value)
        # m2, radius, l, logg= TeffAge_to_Params(T, age.value)
        # Lacc = self.lineflux_to_Lacc(F_Ha, line='Ha', Av=None)
        logLacc = 1.13*np.log10(L_Ha.value) + 1.74
        Lacc = 10**(logLacc) * u.Lsun
        
        # Mdot = Lacc_to_Mdot(Lacc, mass, radius)
        if kwargs.get('uncertainty'):
            if 'disterr' in kwargs:
                disterr = kwargs['disterr']
            elif 'disterrlow' in kwargs:
                disterrlow = abs(self.distance - kwargs['disterrlow'])
                disterrupp = abs(self.distance - kwargs['disterrupp'])
            else:
                raise ValueError('must have either disterr, or disterrlow/disterrupp')
                
            Hamagwerr = ufloat(Hamag, kwargs['Hamagerr'])
            Ftotalwerr = (1.84e-7) * (10**(-0.4*(Hamagwerr + 0.03)))
            Ftotalerr = Ftotalwerr.s * u.erg/u.s/u.cm**2.
      
            EWwerr = ufloat(EW.value, kwargs['EWerr'])
            F_Hawerr = Ftotalwerr * ( (-EWwerr/W.value)/ (1-(EWwerr/W.value)))
            F_Haerr = F_Hawerr.s * u.erg/u.s/u.cm**2.
            
            LaccLsunerr =  []
            for disterr in [disterrlow, disterrupp]:
                L_Ha_err = L_Ha * np.sqrt((F_Haerr/F_Ha)**2. + ((2*disterr)/self.distance)**2.)  # Lsun

                sigmaLogLline = 0.434 * L_Ha_err/L_Ha

                log10L = np.log10(L_Ha.value)
                sigma_log10L = 0.434 * L_Ha_err.value/L_Ha.value
                a,b=1.13, 1.74
                a_err, b_err = 0.05, 0.19
                alog10L = a * log10L
                alog10Lb = a * log10L + b 

                logsigma_Lacc = np.sqrt( b_err**2. + (alog10L * np.sqrt( (a_err/a)**2. + (sigma_log10L/log10L)**2. ))**2.) 

                sigma_Lacc = 2.303 * (10**(alog10Lb)) * logsigma_Lacc
                LaccLsunerr.append(sigma_Lacc)
          
            Lacc_Lsun_err = np.nanmean(LaccLsunerr)* u.Lsun

            if returnFline:
                return [F_Ha, F_Haerr], [L_Ha, L_Ha_err], [Lacc, Lacc_Lsun_err]
            else:
                return Lacc, Lacc_Lsun_err
        else:
            if returnFline:
                return F_Ha, L_Ha, Lacc
            else:
                return Lacc

    def Uband_to_Mdot(self, logLacc, Rin = 5*u.Rsun):
        '''
        rederive U band Mdot value 
        '''
        if isinstance(logLacc, float):
            logLacc = logLacc * u.Lsun
        Lacc = 10**(logLacc.value) * logLacc.unit
        
        Mdot = (1 - (self.radius/Rin))**(-1) * ((Lacc * self.radius) / (const.G * self.mass))
        Mdot = Mdot.to(u.Msun/u.yr)
        return Mdot
    
    def ExcessV_to_Mdot(self, logLacc, Rin = 5*u.Rsun):
        '''
        rederive V band Mdot value 
        '''
        if isinstance(logLacc, float):
            logLacc = logLacc * u.Lsun
        Lacc = 10**(logLacc.value) * logLacc.unit
        
        Mdot = (1 - (self.radius/Rin))**(-1) * ((Lacc * self.radius) / (const.G * self.mass))
        Mdot = Mdot.to(u.Msun/u.yr)
        return Mdot
    


    def logLacc(self, Lacc):
        '''
        convert accretion luminosity to log Lacc
        '''
        if isinstance(Lacc, list):
            L_acc, L_accerr = Lacc
            return np.log10(L_acc.value) * L_acc.unit, 0.434*L_accerr/L_acc
        else:
            return np.log10(Lacc.value) * Lacc.unit


    def Lacc_to_Mdot(self, Lacc, Rin=5*u.Rsun, ignore_Mdot_constants=False, **kwargs):
        '''
        This function will turn an accretion luminosity into a mass accretion rate estimate using the widely
        accepted relationship. The inputs of the function are:
        
        Lacc : astropy.Quantity 
            Accretion luminosity [Lsun]
        Rin : astropy.Quantity, optional
            magnetospheric radius (default 5 [Rsun])
        ignore_Mdot_constants : boolean, optional
            if True, ignore the 1-R/Rin in calculating Lacc (default False)
            
        kwargs:
            uncertainty : boolean
                True to calculate uncertainty.  If True, below quantites must be provided
            Radiuserr : astropy.Quantity 
                error in radius in Rsun
            Masserr : astropy.Quantity 
                error in mass in Msun
        '''
        if kwargs.get('uncertainty'):
            L_acc, L_accerr = Lacc
        else:
            L_acc = Lacc
        if self.radius > 5*u.Rsun:
            Mdot = (self.radius*L_acc)/(const.G*self.mass)

        else:
            Mdot = ((1-(self.radius/Rin))**-1) * (self.radius*L_acc)/(const.G*self.mass)
        Mdot = Mdot.to(u.Msun/u.yr)
        
        if kwargs.get('uncertainty'):
            if 'Radiuserr' in kwargs:
                Radius_Rsun_err = kwargs['Radiuserr'] 
            if 'Masserr' in kwargs:
                Mass_Msun_err = kwargs['Masserr']
            Mdot_err = Mdot * np.sqrt((L_accerr/L_acc)**2. + (Radius_Rsun_err/self.radius)**2. + (Mass_Msun_err/self.mass)**2.)
            return Mdot, Mdot_err
        else:
            return Mdot

class CASPAR_derive_newMdot:
    def __init__(self, row, unc=True, line_list = None):
        if line_list is None:
            self.line_list = ['Ha', 'Hb', 'Hg', 'Hd', 'Hep', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13', 'H14', 'H15', 
                              'PaB', 'PaG', 'PaD', 'Pa8', 'Pa9', 'Pa10', 'BrG', 'HeI 402.6','HeI 447.1','HeI 471.3','HeI 501.6',
                              'HeI 587.6','HeI 667.8','HeI 706.5','HeII','CaII K','CaII H','CaII 849.8','CaII 854.2','CaII 866.2',
                              'NaI 589.0','NaI 589.6','OI 777.3','OI 844.6', 'CIV 154.9']
            
        self.row = row
        self.unc = unc

        self.method = row['Accretion Diagnostic']
        self.tracer = row['Tracer']
        self.distance = row['Distance'] * u.pc
        self.distanceerrlow = row['Distance err lower limit']*u.pc
        self.distanceerrupp = row['Distance err upper limit']*u.pc
        self.mass = row['Mass'] *u.Msun
        self.masserr = row['Mass err'] * u.Msun
        self.radius = row['Radius'] * u.Rsun
        self.radiuserr = row['Radius err'] * u.Rsun

        self.rederive_mdot = rederive_Mdot(self.mass, self.radius, self.distance)

        self.logLacc_orig = row['Original Log Accretion Luminosity'] * u.Lsun
        self.logLacc_origerr = row['Original Log Accretion Luminosity err'] * u.Lsun
        self.Mdot_orig = row['Original Accretion Rate'] * u.Msun / u.yr
        self.mass_orig = row['Original Mass']  * u.Msun
        self.Teff_orig = row['Original Teff'] * u.K

        self.radius_orig = row['Original Radius'] * u.Rsun 
        self.dist_orig = row['Original Distance'] * u.pc
        if np.isnan(self.logLacc_orig):
            self.Lacc_orig = (Mdot_orig * const.G * mass_orig) / (radius_orig*(1-(radius_orig/(5*u.Rsun))))
            self.Lacc_orig = Lacc_orig.to(u.Lsun)
            self.logLacc_orig = np.log10(Lacc_orig.value)*u.Lsun

            
    def rederiveAllLines(self):
        # DO ALL LINE FLUXES 
        updated_vals = {}
        for line in self.line_list:
            if f'{line} Line Flux' in self.row.keys():
                if not np.isnan(self.row[f'{line} Line Flux']):

                    Lacc = self.rederive_mdot.lineflux_to_Lacc(self.row[f'{line} Line Flux'], 
                                                        line=line, 
                                                        uncertainty=self.unc, fluxerr=self.row[f'{line} Line Flux err'],
                                                        disterrlow=self.distanceerrlow, disterrupp = self.distanceerrupp)

                    Mdot = self.rederive_mdot.Lacc_to_Mdot(Lacc, uncertainty=self.unc, 
                                                            Radiuserr=self.radiuserr, Masserr=self.masserr)
                    if self.unc:
                        linelogLaccerr = (0.434 * Lacc[1]/ Lacc[0]).value
                        lineMdoterr = Mdot[1].value
                        linelogLacc = np.log10(Lacc[0].value)
                        lineMdot = Mdot[0].value
                        updated_vals[line] = [[linelogLacc, linelogLaccerr], [lineMdot, lineMdoterr]]
                    else:
                        linelogLacc = np.log10(Lacc.value)
                        lineMdot = Mdot.value
                        updated_vals[line] = [linelogLacc, lineMdot]
        return updated_vals

    
    def rederiveRescale(self):
        LaccFL = self.rederive_mdot.rescale_Lacc(self.logLacc_orig, self.dist_orig)
        AccRateLF = self.rederive_mdot.Lacc_to_Mdot(LaccFL)
        logLaccLF = self.rederive_mdot.logLacc(LaccFL)
        if self.unc:
            logLaccLFerr = self.logLacc_origerr
            LaccLFerr = logLaccLFerr * LaccFL.value / 0.434
            Mdot = self.rederive_mdot.Lacc_to_Mdot([LaccFL, LaccLFerr], uncertainty=self.unc, 
                                                    Radiuserr=self.radiuserr, Masserr=self.masserr)
            return [logLaccLF.value, logLaccLFerr.value], [Mdot[0].value, Mdot[1].value] 
        else:
            return logLaccLF.value, AccRateLF.value
    
    def rederiveComboLineFlux(self):
        formallines = self.tracer.split(', ')
        F2 = [sub.replace('λ', ' ') for sub in formallines]
        lines = [i.split(' ')[0] + ' ' + str(int(i.split(' ')[-1])/10) if ' ' in i  else i for i in F2 ]

        FIN_Lacc = []
        FIN_Lacc_err = []
        
        print(lines)
        for line in lines:
            if f'{line} Line Flux' in self.row.keys():
                if not np.isnan(self.row[f'{line} Line Flux']):
                    Lacc = self.rederive_mdot.lineflux_to_Lacc(self.row[f'{line} Line Flux'], 
                                                    line=line, uncertainty=unc,
                                                        fluxerr=self.row[f'{line} Line Flux err'],
                                                    disterrlow=self.distanceerrlow, disterrupp = self.distanceerrupp)

                    if self.unc:
                        FIN_Lacc.append(Lacc[0].value)
                        FIN_Lacc_err.append(Lacc[1].value)
                    else:
                        FIN_Lacc.append(Lacc.value)
        if len(FIN_Lacc) == 0:
            return 'no line fluxes in database', -1000
        else:
            LaccFL = np.nanmean(np.array(FIN_Lacc))* u.Lsun
            logLaccLF =  self.rederive_mdot.logLacc(LaccFL)
            AccRateLF = self.rederive_mdot.Lacc_to_Mdot(LaccFL)

            if self.unc:
                if np.nanstd(np.array(FIN_Lacc)) != 0:
                    LaccLFerr =  np.nanstd(np.array(FIN_Lacc))* u.Lsun
                else:
                    LaccLFerr =  np.nanmean(np.array(FIN_Lacc_err))* u.Lsun
                logLaccLFerr = (0.434* LaccLFerr/ LaccFL)

                Mdot = self.rederive_mdot.Lacc_to_Mdot([LaccFL, LaccLFerr], uncertainty=self.unc, 
                                                        Radiuserr=self.radiuserr, Masserr=self.masserr)
                return [logLaccLF.value, logLaccLFerr], [Mdot[0].value, Mdot[1].value] 
            else:
                return logLaccLF.value, AccRateLF.value

    def rederiveFromLaccMdot(self, origab=None, ignore_Mdot_constants=True):
        if origab is None:
            origab = {'Ha':[1, 1.72], 
                'BrG':[0.9, 2.9],
                'PaB':[1.36, 4.],
                'PfB':[0.91, 3.29],
                'PaG':[1.36, 4.1]} ##PaG=Gatti2008
                
        formallines = self.tracer.split(', ')
        F2 = [sub.replace('λ', ' ') for sub in formallines]
        lines = [i.split(' ')[0] + ' ' + str(int(i.split(' ')[-1])/10) if ' ' in i  else i for i in F2 ]
        if len(lines) > 0:
            raiseValueError('can only be 1 line')
        else:
            A,B = origab[line][0], origab[line][1]
            print(A,B)
            if not np.isnan(self.logLacc_orig.value):
                line_flux = self.rederive_mdot.Lacc_to_lineflux(self.logLacc_orig, self.dist_orig, A, B)
                print(line_flux, self.logLacc_orig)
            else: # ORIG MDOT
                newlogLacc, line_flux = self.rederive_mdot.Mdot_to_lineflux(self.Mdot_orig, self.mass_orig, 
                                                                           self.radius_orig, self.dist_orig, A,B, ignore_Mdot_constants=True)
                print(line_flux, newlogLacc)
            newLacc = self.rederive_mdot.lineflux_to_Lacc(line_flux, line=line)
            LogLaccLF = self.rederive_mdot.logLacc(newLacc)
            AccRateLF = self.rederive_mdot.Lacc_to_Mdot(newLacc)

            return line_flux.value, LogLaccLF.value, AccRateLF.value
    
    def rederiveHaphot(self):
        if self.row['Original Reference'] != 'Kalari 2015':
                raise ValueError('This attribute can only used with Kalari 2015 items')
        else:
            line = 'Ha'    
            EW = float(self.row['Ha EW']) * u.AA
            Hamag = self.row['Ha mag'] - 0.8
            F_Ha, L_Ha, Lacc = self.rederive_mdot.photHaEW_to_Lacc(EW, Hamag, returnFline=True,
                                                                                            uncertainty=self.unc, disterrlow=self.distanceerrlow, 
                                                                                            disterrupp = self.distanceerrupp,
                                                                                            Hamagerr=self.row['Ha mag err'], EWerr=self.row['Ha EW err'])

            logLaccphot = self.rederive_mdot.logLacc(Lacc)
            AccRatephot = self.rederive_mdot.Lacc_to_Mdot(Lacc, uncertainty=self.unc, 
                                                            Radiuserr=self.radiuserr, Masserr=self.masserr)
            if unc:     
                return [[F_Ha[0].value, F_Ha[1].value], [logLaccphot[0].value, 0.434* Lacc[1]/Lacc[0]], 
                        [AccRatephot[0].value, AccRatephot[1].value]]
            else:
                return F_Ha[0].value, logLaccphot.value, AccRatephot.value

    def rederiveRun(self, origab=None, ignore_Mdot_constants=True):
        # do all lines
        updated_vals = self.rederiveAllLines()
        lineFlux = np.nan
        if self.method == 'Line Luminosity':
            if self.tracer == 'Ha photometry':
                lineFlux, logaccLum, AccRate = self.rederiveHaphot()
            elif self.row['Scaling Relation Reference'] == 'Alcala 2017':
                logaccLum, AccRate = self.rederiveRescale()

            else:
                logaccLum, AccRate = self.rederiveComboLineFlux()
                if logaccLum == 'no line fluxes in database':
                    lineFlux, logaccLum, AccRate =  self.rederiveFromLaccMdot(origab=None, ignore_Mdot_constants=True)
            
            if self.tracer == 'PfB':
                scaleRef = 'Salyk 2009'
            else:
                scaleRef = 'Alcala 2017'
            
        elif self.method == 'Excess Continuum':
            logaccLum, AccRate = self.rederiveRescale()
            scaleRef = ' '

        elif self.method == 'Line Profile':
            logaccLum, AccRate, scaleRef = self.logLacc_orig, self.Mdot_orig, 'Natta 2004'
        
        return updated_vals, logaccLum, AccRate, scaleRef, lineFlux
        

