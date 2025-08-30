# casparing
---
casparing is the toolkit used to rederive and extract parameters for the Comprehensive Archive of SubStellar and Planetary Accretion Rates (CASPAR). 

There are tools available to extract SIMBAD ra, dec, magnitudes, and GAIA distances and parallaxes.  Additionally, the toolkit uses Baraffe+2015 and MIST MESA evolutionary models to extract physical parameters such as stellar Mass and Radius. Finally, with a known accretion line flux or luminosity accretion rates can be rederived to be self consistent with the rest of CASPAR. 

---
## Installation 
To use the module, download or clone the repository above and install it with the following 
```
$ git clone https://github.com/sbetti22/casparing.git
$ cd casparing
$ pip install -e .
````

After installation, install ```linmix``` [https://linmix.readthedocs.io/en/latest/install.html](https://linmix.readthedocs.io/en/latest/install.html) following the specific installation instructions.


## Requirements
```casparing``` primarily uses the ```pandas``` package as well as several astronomical and data specific packages.
  - astroquery
  - isochrones
  - extinction
  - uncertainties
  - pandarallel
  - [astropy](https://www.astropy.org)
  - numpy 
  - pandas 
  - matplotlib 
  - astroquery 
  - tqdm
  - pandarallel
  - scipy
  - uncertainties 
  - extinction 
  - isochrones
  - gdown

```casparing``` additionally uses Banyan Σ [https://github.com/jgagneastro/banyan_sigma](https://github.com/jgagneastro/banyan_sigma).  While easy to install on Intel chips, we found issues when installing on apple M1 chips.  Therefore, we have copied the banyan_sigma.core file and banyan_sigma.data folder to ```casparing``` in order to use it. Those file and data are from commit 84391a8. 


## To Run
Example notebooks are in the works.  Please check back from tutorials.  

To run any of the plotting notebooks, you need to download CASPAR and the Literature Database as a .csv file which can be found at [https://doi.org/10.5281/zenodo.8393053](https://doi.org/10.5281/zenodo.8393053) or through [https://ui.adsabs.harvard.edu/abs/2023AJ....166..262B/abstract](https://ui.adsabs.harvard.edu/abs/2023AJ....166..262B/abstract). 

Then, run 
```
from casparing.plot import CASPAR_sortdata as csort
lit_database = 'Comprehensive Archive of Substellar and Planetary Accretion Rates (CASPAR) Betti+2023 - Literature Database.csv'
caspar = 'Comprehensive Archive of Substellar and Planetary Accretion Rates (CASPAR) Betti+2023 - CASPAR.csv'
df_lit, df_caspar = csort.CASPAR_loaddata(lit_database=lit_database, caspar=caspar)
```
With this, you can then use any of the other plotting modules.
```
from casparing.plot import CASPAR_plot_variability as cvar
cvar.plot_allvariability(df_caspar)
```

To run the rederivation, you need to start with known quantities, such as the object name, age, and spectral type.   
```
from casparing.derive import CASPAR_astrometry as cmetry
from casparing.derive import CASPAR_rederiveparameters as cderpar
df = pd.DataFrame({'obj':['CI Tau', 'CQ Tau']})
table = df['obj'].progress_apply(cmetry.query_simbad)
tab = pd.concat([t for t in table], ignore_index=True)
cmetry.run_banyan(tab)
cmetry.get_ages(tab)

newdf = cderpar.create_table(object_name='CI Tau', age=2, age_err=0.32, 
                             value_id='sptype', value='M7', value_err=0.5)
cdp = cderpar.get_physical_params(df=tab)
cdp.extractSptypeTempOrMass()
cdp.extractAge()
cdp.params()
```


## Attribution

Please cite Betti et al. (2023) whenever results from CASPAR or ```casparing``` are used in a publication as well as the many relevant papers and catalogues that the derivations use, including 2MASS, Gaia, Alcala et al., (2017), and Banyan Σ.    


## License

Copyright 2025 Sarah Betti.

```casparing``` is distributed under the MIT License. See LICENSE for the terms and conditions.
