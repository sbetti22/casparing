import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import astroquery
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack

from astropy import units as u
from astropy import constants as const
from tqdm import tqdm

import logging
logger = logging.getLogger("astroquery")
logger.setLevel(logging.WARNING)

import warnings
warnings.filterwarnings('ignore')

from banyansigma_core import membership_probability

def query_simbad(obj, verbose=False):
    '''
    inputs:
        obj: str -- name of object to query
    '''
    if verbose:
        print('object: ' + obj)
        
    # initialize SIMBAD and find values of interest
    customSimbad = Simbad()
    customSimbad.add_votable_fields('flux')


    # search though SIMBAD within radius arcseconds around given coordinate
    result_table = customSimbad.query_object(obj)
    # if no object, return table with just name
    if result_table is None:
        if verbose:
            print('No results found with SIMBAD')
        fin_result_table = pd.DataFrame({'Reference Name':[obj]})
        return fin_result_table
    else:
        rt = result_table[['main_id', 'ra', 'dec', 'flux.filter', 'flux', 'flux_err']].to_pandas()
        collapsed_rt = pd.pivot_table(rt, values=['flux', 'flux_err'], index=['main_id', 'ra', 'dec'],
                       columns=['flux.filter'], aggfunc="min", fill_value=np.nan)
        fin_result_table = collapsed_rt.reset_index()
        fin_result_table.columns = [f'{col[1]} {col[0]}'.strip() for col in fin_result_table.columns]
        
        if 'J flux' in list(fin_result_table.columns):
            fin_result_table.rename(columns={'J flux':'Jmag'}, inplace=True)
        if 'J flux_err' in list(fin_result_table.columns):
            fin_result_table.rename(columns={'J flux_err':'Jmag err'}, inplace=True)
            
        if 'H flux' in list(fin_result_table.columns):
            fin_result_table.rename(columns={'H flux':'Hmag'}, inplace=True)
        if 'H flux_err' in list(fin_result_table.columns):
            fin_result_table.rename(columns={'H flux_err':'Hmag err'}, inplace=True)
            
        if 'K flux' in list(fin_result_table.columns):
            fin_result_table.rename(columns={'K flux':'Kmag'}, inplace=True)
        if 'K flux_err' in list(fin_result_table.columns):
            fin_result_table.rename(columns={'K flux_err':'Kmag err'}, inplace=True)


        # query object
        names = Simbad.query_objectids(obj)
        # get 2MASS name 
        mass2 = [N for N in names['id'] if '2MASS' in N ]
        link = 'https://simbad.u-strasbg.fr/simbad/sim-id?Ident='
        if len(mass2) > 0:
            fin_result_table['Simbad-Resolvable Name'] = [mass2[0]]
            sname = mass2[0].replace('+', '%2B').replace(' ', '+')
            fin_result_table['Links'] = [link+sname]
        elif len(names) > 0:
            fin_result_table['Simbad-Resolvable Name'] = [names['id'][0]]
            sname = names['id'][0].replace('+', '%2B').replace(' ', '+')
            fin_result_table['Links'] = [link+sname]
            
        else:
            fin_result_table['Simbad-Resolvable Name'] = []
            fin_result_table['Links'] = []
        # get Gaia DR2 and Gaia EDR3 source IDs
        gaia2 = [N for N in names['id'] if 'Gaia DR2' in N ]
        gaia3 = [N for N in names['id'] if 'Gaia EDR3' in N ]
        if len(gaia2) > 0:
            gaia2_source_id = gaia2[0].split(' ')[-1]
        else:
            gaia2_source_id = 'none'
        if len(gaia3) > 0:
            gaia3_source_id = gaia3[0].split(' ')[-1]
        else:
            gaia3_source_id = 'none'
        # get RA and Dec J2000 
        RA = str(fin_result_table['ra'].values[0])
        DEC = str(fin_result_table['dec'].values[0])
        RA, DEC
        c = SkyCoord(RA+' ' +DEC, unit=(u.deg, u.deg), frame='icrs')
        fin_result_table['RA (J2000.0)'] = [c.ra.deg]
        fin_result_table['Dec (J2000.0)'] = [c.dec.deg]
        
        ########################## GAIA DR 2 ###############################
        # search Gaia DR2 database
        from astroquery.gaia import Gaia
        radius = u.Quantity(5.0, u.arcsecond)
        Gaia.MAIN_GAIA_TABLE = "gaiadr2.gaia_source" 
        j = Gaia.cone_search_async(coordinate=c, radius=radius)
        r = j.get_results()
        # if it finds an object, confirm its id matches simbad
        if len(r) > 0:
            if verbose:
                if str(r['source_id'][0]) == gaia2_source_id:
                    print('  --> DR2 IDs match: gaia=', r['source_id'][0], 'simbad=', gaia2_source_id )
                else:
                    print('  --> DR2 no match: gaia=', r['source_id'][0], 'simbad=', gaia2_source_id )
            # SQL query to get the object properties
            qry = "select * from gaiadr2.gaia_source as gaia WHERE gaia.source_id = " + str(r['source_id'][0])
            job = Gaia.launch_job( qry, verbose=False )
            tblGaiaDR2 = job.get_results()      #Astropy table
            fin_result_table['GAIA DR2 Source ID'] = [str(r['source_id'][0])]
            # determine if the parallax is usable to get distance
            if tblGaiaDR2['parallax'][0]*0.25 > tblGaiaDR2['parallax_error'][0]:
                rel_par = 'True'
            else:
                rel_par = ' '
            # get all property information
            fin_result_table['GAIA DR2 Reliable Parallax'] = [rel_par]
            fin_result_table['GAIA DR2 Parallax'] = tblGaiaDR2['parallax']
            fin_result_table['GAIA DR2 Parallax err'] = tblGaiaDR2['parallax_error']
            fin_result_table['GAIA DR2 RA proper motion'] = tblGaiaDR2['pmra']
            fin_result_table['GAIA DR2 RA proper motion err'] = tblGaiaDR2['pmra_error']
            fin_result_table['GAIA DR2 Dec proper motion'] = tblGaiaDR2['pmdec']
            fin_result_table['GAIA DR2 Dec proper motion err'] = tblGaiaDR2['pmdec_error']
            fin_result_table['GAIA DR2 Dec proper motion err'] = tblGaiaDR2['pmdec_error']
            fin_result_table['radial velocity'] = tblGaiaDR2['radial_velocity']
            fin_result_table['radial velocity err'] = tblGaiaDR2['radial_velocity_error']
            
            # extract distances from Bailer-Jones distance catalogue
            if rel_par== 'True':
                qry = "select * from external.gaiadr2_geometric_distance as dist WHERE dist.source_id = " + str(tblGaiaDR2['source_id'][0])
                job = Gaia.launch_job( qry, verbose=False )
                tblGaiaDR2_dist = job.get_results()      #Astropy table
                fin_result_table['GAIA DR2 Distance'] = tblGaiaDR2_dist['r_est']
                fin_result_table['GAIA DR2 Distance lower limit'] = tblGaiaDR2_dist['r_lo']
                fin_result_table['GAIA DR2 Distance upper limit'] = tblGaiaDR2_dist['r_hi']
            else:
                fin_result_table['GAIA DR2 Distance'] = [np.ma.masked]
                fin_result_table['GAIA DR2 Distance lower limit'] = [np.ma.masked]
                fin_result_table['GAIA DR2 Distance upper limit'] = [np.ma.masked]
            
        else:
            # if nothing, return empty table
            fin_result_table['GAIA DR2 Source ID'] = [' ']
            fin_result_table['GAIA DR2 Reliable Parallax'] = [' ']
            fin_result_table['GAIA DR2 Parallax'] = [np.ma.masked]
            fin_result_table['GAIA DR2 Parallax err'] = [np.ma.masked]
            fin_result_table['GAIA DR2 RA proper motion'] = [np.ma.masked]
            fin_result_table['GAIA DR2 RA proper motion err'] = [np.ma.masked]
            fin_result_table['GAIA DR2 Dec proper motion'] = [np.ma.masked]
            fin_result_table['GAIA DR2 Dec proper motion err'] =[ np.ma.masked]
            fin_result_table['GAIA DR2 Dec proper motion err'] = [np.ma.masked]
            fin_result_table['radial velocity'] = [np.ma.masked]
            fin_result_table['radial velocity err'] = [np.ma.masked]
            fin_result_table['GAIA DR2 Distance'] =[ np.ma.masked]
            fin_result_table['GAIA DR2 Distance lower limit'] = [np.ma.masked]
            fin_result_table['GAIA DR2 Distance upper limit'] = [np.ma.masked]
            

    
        #################### GAIA EDR 3 #########################################
        # search Gaia EDR3 database
        radius = u.Quantity(5.0, u.arcsecond)
        Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
 
        jj = Gaia.cone_search_async(coordinate=c, radius=radius)
        rr = jj.get_results()
        # if it finds an object, confirm its id matches simbad
        if len(rr) > 0:
            if verbose:
                if str(rr['source_id'][0]) == gaia3_source_id:
                    print('  --> EDR3 IDs match: gaia=', rr['source_id'][0], 'simbad=', gaia3_source_id )
                else:
                    print('  --> EDR3 no match: gaia=', rr['source_id'][0], 'simbad=', gaia3_source_id )
            # SQL query to get the object properties
            qry = "select * from gaiaedr3.gaia_source as gaia WHERE gaia.source_id = " + str(rr['source_id'][0])
            job = Gaia.launch_job( qry )
            tblGaiaEDR3 = job.get_results()      #Astropy table
            fin_result_table['GAIA EDR3 Source ID'] = [str(rr['source_id'][0])]
            fin_result_table['RA (J2016.0)'] = tblGaiaEDR3['ra']
            fin_result_table['Dec (J2016.0)'] = tblGaiaEDR3['dec']
            # determine if the parallax is usable to get distance
            if tblGaiaEDR3['parallax'][0]*0.25 > tblGaiaEDR3['parallax_error'][0]:
                rel_par = 'True'
            else:
                rel_par = ' '
            # get all property information
            fin_result_table['GAIA EDR3 Reliable Parallax'] = [rel_par]
            fin_result_table['GAIA EDR3 Parallax'] = tblGaiaEDR3['parallax']
            fin_result_table['GAIA EDR3 Parallax err'] = tblGaiaEDR3['parallax_error']
            fin_result_table['GAIA EDR3 RA proper motion'] = tblGaiaEDR3['pmra']
            fin_result_table['GAIA EDR3 RA proper motion err'] = tblGaiaEDR3['pmra_error']
            fin_result_table['GAIA EDR3 Dec proper motion'] = tblGaiaEDR3['pmdec']
            fin_result_table['GAIA EDR3 Dec proper motion err'] = tblGaiaEDR3['pmdec_error']
            fin_result_table['GAIA EDR3 Dec proper motion err'] = tblGaiaEDR3['pmdec_error']
            # extract distances from Bailer-Jones distance catalogue
            if rel_par=='True':
                qry = "select * from external.gaiaedr3_distance as dist WHERE dist.source_id = " + str(tblGaiaEDR3['source_id'][0])
                job = Gaia.launch_job( qry )
                tblGaiaEDR3_dist = job.get_results()      #Astropy table
                fin_result_table['GAIA EDR3 Geometric Distance'] = tblGaiaEDR3_dist['r_med_geo']
                fin_result_table['GAIA EDR3 Geometric Distance lower limit'] = tblGaiaEDR3_dist['r_lo_geo']
                fin_result_table['GAIA EDR3 Geometric Distance upper limit'] = tblGaiaEDR3_dist['r_hi_geo']
            else:
                fin_result_table['GAIA EDR3 Geometric Distance'] = [np.ma.masked]
                fin_result_table['GAIA EDR3 Geometric Distance lower limit'] = [np.ma.masked]
                fin_result_table['GAIA EDR3 Geometric Distance upper limit'] = [np.ma.masked]

            fin_result_table['Reference Name'] = [str(obj)]
        else:
            # if nothing, return empty table
            fin_result_table['GAIA EDR3 Source ID'] = [' ']
            fin_result_table['GAIA EDR3 Reliable Parallax'] = [' ']
            fin_result_table['GAIA EDR3 Parallax'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 Parallax err'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 RA proper motion'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 RA proper motion err'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 Dec proper motion'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 Dec proper motion err'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 Dec proper motion err'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 Geometric Distance'] = [np.ma.masked]
            fin_result_table['GAIA EDR3 Geometric Distance lower limit'] =[ np.ma.masked]
            fin_result_table['GAIA EDR3 Geometric Distance upper limit'] = [np.ma.masked]

        return fin_result_table[['Reference Name','Simbad-Resolvable Name', 'RA (J2000.0)', 'Dec (J2000.0)', 'RA (J2016.0)', 'Dec (J2016.0)',
                            'GAIA DR2 Source ID','GAIA DR2 Parallax', 'GAIA DR2 Parallax err', 
                             'GAIA DR2 Reliable Parallax', 'GAIA DR2 Distance',
                             'GAIA DR2 Distance lower limit', 'GAIA DR2 Distance upper limit',
                             'GAIA DR2 RA proper motion', 'GAIA DR2 RA proper motion err', 
                             'GAIA DR2 Dec proper motion', 'GAIA DR2 Dec proper motion err',
                             'GAIA EDR3 Source ID','GAIA EDR3 Parallax', 'GAIA EDR3 Parallax err', 
                             'GAIA EDR3 Reliable Parallax', 'GAIA EDR3 Geometric Distance',
                             'GAIA EDR3 Geometric Distance lower limit', 'GAIA EDR3 Geometric Distance upper limit',
                             'GAIA EDR3 RA proper motion', 'GAIA EDR3 RA proper motion err', 
                             'GAIA EDR3 Dec proper motion', 'GAIA EDR3 Dec proper motion err',
                             'radial velocity', 'radial velocity err', 
                             'Jmag', 'Jmag err', 
                             'Hmag', 'Hmag err', 'Kmag', 'Kmag err', 'Links']]
    

def run_banyan(df):
    '''
    A function that takes an input and output dataframe and runs the banyan_sigma code from jgagneastro
    Keyword inputs:
    df - CASPAR input for query

    inputs - list of CASPAR columns with EDR3 data and RV and ERV data
    DR2_inputs - list of CASPAR columns with DR2 versions of EDR3 data
    '''
    inputs = ['RA (J2016.0)','Dec (J2016.0)','GAIA EDR3 RA proper motion','GAIA EDR3 Dec proper motion',
          'GAIA EDR3 RA proper motion err','GAIA EDR3 Dec proper motion err',
          'radial velocity','radial velocity err',
          'GAIA EDR3 Parallax','GAIA EDR3 Parallax err']
    
    DR2_inputs = ['RA (J2000.0)','Dec (J2000.0)','GAIA DR2 RA proper motion','GAIA DR2 Dec proper motion',
                  'GAIA DR2 RA proper motion err','GAIA DR2 Dec proper motion err']

    BSnames = {'118TAU':'118 Tau', 'ABDMG':'AB Doradus', 'BPMG':'B Pictoris', 
           'CAR':'Carina', 'CARN':'Carina-Near', 'CBER':'Coma Berenices',
           'COL':'Columba', 'CRA':'Corona Australis', 'EPSC':'e Chamaeleontis',
           'ETAC':'n Chamaeleontis', 'HYA':'Hyades', 'LCC':'Lower Centaurus Crux',
           'OCT':'Octans', 'PL8':'Platais 8', 'PLE':'Pleides Cluster',
           'ROPH':'p Ophiuchus', 'THA':'Tucana-Horologium', 'THOR':'32 Orionis', 
           'TWA':'TW Hya', 'UCL':'Upper Centaurus Lupus', 'UCRA':'Upper CrA',
           'UMA':'Ursa Major cluster', 'USCO':'Upper Scorpius', 'TAU':'Taurus', 
           'XFOR':'X For', 'FIELD':'FIELD'}

    # loop through each object and run it through BANYAN Σ.
    for ind in tqdm(df.index):
        warnings.filterwarnings(action="ignore")        
        dic = {}
        #fill dictionary with key as EDR3 columns and value as EDR3 data
        for col in inputs:
            dic[col] = float(df.loc[ind,col])
        #make standard input GAIA ERD3 data
        input0 = dic[inputs[0]] # RA
        input1 = dic[inputs[1]] # Dec
        input2 = dic[inputs[2]] # pm ra
        input3 = dic[inputs[3]] # pm dec
        input4 = dic[inputs[4]] # pm ra err
        input5 = dic[inputs[5]] # pm de cerr
        #these won't change, we only have 1 version of each data here
        input6 = dic[inputs[6]] # rad vel
        input7 = dic[inputs[7]] # rad vel err
        input8 = dic[inputs[8]] # parallax
        input9 = dic[inputs[9]] # parallax err
        #Check for required inputs using EDR3 data, if nan check DR2 and use, if nan move to next object
        #RA check
        if np.isnan(input0):
            DR2_RA = float(df.loc[ind,DR2_inputs[0]])
            if not np.isnan(DR2_RA):
                input0 = DR2_RA
            else:
                input0 = None
        #Dec check
        if np.isnan(input1):
            DR2_Dec = float(df.loc[ind,DR2_inputs[1]])
            if not np.isnan(DR2_Dec):
                input1 = DR2_Dec
            else:
                input1 = None
        #PMRA check
        if np.isnan(input2):
            DR2_PMRA = float(df.loc[ind,DR2_inputs[2]])
            if not np.isnan(DR2_PMRA):
                input2 = DR2_PMRA
            else:
                input2 = None
        #PMDEC check
        if np.isnan(input3):
            DR2_PMDec = float(df.loc[ind,DR2_inputs[3]])
            if not np.isnan(DR2_PMDec):
                input3 = DR2_PMDec
            else:
                input3 = None
        #EPMRA check
        if np.isnan(input4):
            DR2_PMRA_err = float(df.loc[ind,DR2_inputs[4]])
            if not np.isnan(DR2_PMRA_err):
                input4 = DR2_PMRA_err
            else:
                input4 = None
        #EPMDEC check
        if np.isnan(input5):
            DR2_PMDec_err = float(df.loc[ind,DR2_inputs[5]])
            if not np.isnan(DR2_PMDec_err):
                input5 = DR2_PMDec_err
            else:
                input5 = None
        #run banyan ∑
        #if radial velocity or radial velocity err is nan, make both nan, can't have 1/2
        if np.isnan(input6) or np.isnan(input7):
            input6 = np.nan
            input7 = np.nan
        # #if parallax or parallax err is nan, make both nan, can't have 1/2
        if np.isnan(input8) or np.isnan(input9):
            input8 = np.nan
            input9 = np.nan
        
        # if any value (RA, DEC, pms) are neg, can't run the code, so just return nans
        if input3 is None:
            df.loc[ind,'Association'] = np.nan
            df.loc[ind,'Association Probability Banyan Sigma'] = np.nan
        #run banyan sigma with valid inputs, collect output in df
        else:
            result = membership_probability(ra=input0,dec=input1,pmra=input2,pmdec=input3,epmra=input4,epmdec=input5,plx=input8,eplx=input9)
            
            ## OLD 
#             result = bs.banyan_sigma(ra=input0,dec=input1,pmra=input2,pmdec=input3,epmra=input4,epmdec=input5, plx=input8,eplx=input9)
            # if the FIELD has a prob > 70% use it, otherwise, use the association
            if result['ALL']['FIELD'][0] > 0.7:
                Assoc = 'FIELD'
                PROB = result['ALL']['FIELD'][0]
            else:
                Assoc = result['BEST_YA']['Global'][0]
                PROB = result['ALL'][result['BEST_YA']['Global'][0]][0]
            
            df.loc[ind,'Association'] = BSnames[Assoc]
            df.loc[ind,'Association Probability Banyan Sigma'] = PROB
            df.loc[ind,'Association Census Reference'] = 'Gange2018'



def get_ages(tab, verbose=False):
    ages = pd.read_csv('../data/parameters_data/Association Ages and Distances.csv')
    df2 = {'Association': ['FIELD'], 'Age (Myr)': [np.nan], 'Age err (Myr)': [np.nan]}
    df2 = pd.DataFrame(df2)
    ages = pd.concat([ages,df2], ignore_index = True)
   
    for index, row in tqdm(tab.iterrows()):
        ass = row['Association']
        if row['Simbad-Resolvable Name']==u' ':
            name = row['Reference Name']
        else:
            name = row['Simbad-Resolvable Name']    
        if verbose: print(name, ass)
        
        if ass=='Taurus':
            Taurus_groups = {'C1':'Taurus/L1551', 'C2':'Taurus/L1495', 'C3':'Taurus/L1517-Halo', 'C4':'Taurus/L1517-Center', 
                         'C5':"Taurus/L1546", 'C6':'Taurus/L1524', 'C7':'Taurus/L1527', 'C8':'Taurus/B213', 'C9':'Taurus/118TauE',
                         'C10':'Taurus/118TauW', 'D1':'Taurus/L1544', 'D2':'Taurus/L1558', 'D3':'Taurus/South',
                         'D4':'Taurus/North', 'D5':'Taurus/Central', 'D6':'Taurus/North', 'D7':'Taurus/East'}

            v = Vizier(columns=["Name", "GMM"], catalog="J/AJ/162/110/table1")
            result = v.query_region(name, radius="1s")
            if len(result)>0:
                group = result[0].to_pandas()['GMM'].map(Taurus_groups)   
                if verbose: print('  --> Taurus:', group.values[0])
                tab.at[index, 'Association'] = group.values[0]
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
            else:
                if verbose: print(' -->', row['Association'])
                if verbose: print(ages['Age (Myr)'].loc['Taurus '==ages['Association']])
        
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc['Taurus '==ages['Association']].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc['Taurus ' == ages['Association']].values[0]                           


        elif ass == 'p Ophiuchus':
            pOph_groups = {'off':'p Oph', 'L1688':'p Oph/L1688', 'RhoOph':'p Oph', 'L1689':'p Oph/L1689', 
              'L1709':'p Oph/L1709'}

            v = Vizier(columns=["2MASS", "Region"], catalog="J/AJ/159/282/table1.dat")
            result = v.query_region(name, radius="1s")
            if len(result)>0:
                
                group = result[0].to_pandas()['Region'].map(pOph_groups)
                if verbose: print(' --> p Oph', group.values[0])
                tab.at[index, 'Association'] = group.values[0]
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
            else:
                if verbose: print(' --> ', row['Association'])
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==row['Association']].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==row['Association']].values[0]                           
        
        elif (ass == 'FIELD') or (ass==u''):
            good = 0
            
            # LUPUS
            v = Vizier(columns=["GaiaDR2", "Cloud"], catalog="J/A+A/643/A148/tablea1")
            result = v.query_region(name, radius="1s")
            if len(result)>0:
                group = result[0].to_pandas()['Cloud']
                if verbose: print('  --> lupus', group.values[0])
                tab.at[index, 'Association'] = group.values[0]
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
                tab.at[index, 'Association Census Reference'] = 'Galli2020'
                tab.at[index, 'Association Probability Banyan Sigma'] = np.nan


                good = 1
                
            # σ Ori, NGC 2024, Orion Nebular Cluster, 25 Ori
            if good == 0:
                v = Vizier(columns=["2MASS", "Group"], catalog="J/AJ/156/84/table1")
                result = v.query_region(name, radius="5s")
                if len(result)>0:
                    group = result[0].to_pandas()['Group'].str.split('-').str[0]
                    if (group.values[0] == 'NGC 2024') or (group.values[0] == 'Orion Nebular Cluster') or (group.values[0]=='25Ori'):
                        group = result[0].to_pandas()['Group']#.str.split('-').str[0]
                    else:
                        group = ages['Association'].loc[ages['Association']=='σ Ori']
    
                    if verbose: print('  --> sig Ori:', group.values[0])
                    tab.at[index, 'Association'] = group.values[0]
                    tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                    tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
                    tab.at[index, 'Association Census Reference'] = 'Kounkel2018'
                    tab.at[index, 'Association Probability Banyan Sigma'] = np.nan
                    good = 1
            # IC 348 & sigma Ori
            if good == 0:
                ic348_groups = {32:'IC 348', 93:'σ Ori', 90:'σ Ori'}
                v = Vizier(columns=["*"], catalog="J/MNRAS/434/806/members")
                result = v.query_region(name, radius="5s")
                
                if len(result)>0:
                    group = result[0].to_pandas()['Field'].astype(int).map(ic348_groups)
                    if verbose: 
                        if group.values[0] == 'IC 348':
                            print('  --> ic348:', group.values[0])
                        elif group.values[0] == 'σ Ori':
                            print('  --> sig Ori:', group.values[0])
                        else:
                            print(result[0].to_pandas()['Field'])
                    tab.at[index, 'Association'] = group.values[0]
                    tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                    tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
                    tab.at[index, 'Association Census Reference'] = 'Bell2013'
                    tab.at[index, 'Association Probability Banyan Sigma'] = np.nan
                    good = 1
            # Cha I
            if good ==0:
                v = Vizier(columns=["Name", "Cloud"], catalog="J/A+A/646/A46/tablea7")
                result = v.query_region(name, radius="1s")
                if len(result)>0:
                    group = ages['Association'].loc[ages['Association']=='Chameleon I']
                    if verbose: print('  --> chai', group.values[0])
                    tab.at[index, 'Association'] = group.values[0]
                    tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                    tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
                    tab.at[index, 'Association Census Reference'] = 'Galli 2021'
                    tab.at[index, 'Association Probability Banyan Sigma'] = np.nan
                    good = 1
                else:
                    v = Vizier(columns=["Name", "Cloud"], catalog="J/AJ/154/46")
                    result = v.query_region(name, radius="1s")
                    if len(result)>0:
                        group = ages['Association'].loc[ages['Association']=='Chameleon I']
                        if verbose: print('  --> chai 2', group.values[0])
                        tab.at[index, 'Association'] = group.values[0]
                        tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                        tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                           
                        tab.at[index, 'Association Census Reference'] = 'Galli 2021'
                        tab.at[index, 'Association Probability Banyan Sigma'] = np.nan
                        good = 1
            # FIELD
            if good ==0:
                group = ages['Association'].loc[ages['Association']=='FIELD']
                if verbose: print('  --> FIELD')
                good = 1
                tab.at[index, 'Association'] = group.values[0]
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                                   
                tab.at[index, 'Association Census Reference'] = ' '
                tab.at[index, 'Association Probability Banyan Sigma'] = np.nan
        # FIELD
        elif isinstance(ass, float):
            if verbose: print('  --> FIELD')
            group = ages['Association'].loc[ages['Association']=='FIELD']
            tab.at[index, 'Association'] = group.values[0]
            tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==group.values[0]].values[0]
            tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==group.values[0]].values[0]                                      
            tab.at[index, 'Association Census Reference'] = ' '
            tab.at[index, 'Association Probability Banyan Sigma'] = np.nan
        else: 
            if verbose: print('other  --> ', row['Association'])
            if row['Association'] in ages['Association'].values:
                if verbose: print('age',ages['Age (Myr)'].loc[ages['Association']==row['Association']].values[0])
                tab.at[index, 'Association Age'] = ages['Age (Myr)'].loc[ages['Association']==row['Association']].values[0]
                tab.at[index, 'Association Age err'] = ages['Age err (Myr)'].loc[ages['Association']==row['Association']].values[0]                           
            else:
                if verbose: print('nan ages')
                tab.at[index, 'Association Age'] = np.nan
                tab.at[index, 'Association Age err'] = np.nan                           


        if verbose: print('---')
    return tab


def save_df(df, savename):
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
    df.to_csv(savename,
          columns=columns)
    

            