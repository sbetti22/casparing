### fitting codes
from scipy import optimize
import emcee, corner
from scipy.optimize import curve_fit
import numpy as np
import statsmodels.api as sm
import linmix
import corner
import pandas as pd
import matplotlib.pyplot as plt
from scipy import odr as odr
import os

HERE = os.path.dirname(os.path.abspath(__file__))
fil = os.path.join(HERE, 'Linear_fitting_ALLValues_LINMAX.csv')
linfit = pd.read_csv(fil, delimiter=',', comment='#', skiprows=[0])

linfit = linfit.set_index('obj')
linfit = linfit.replace(-1000, np.nan)

## MCMC
def log_likelihood(params, x, obs):
    """ Calculate log likelihood assuming iid Gaussian errors.
    """    
    # Extract parameter values
    a_est, b_est, sigma_est = params
    
    # Calculate deterministic results with these parameters
    sim = a_est*x + b_est
    
    # Calculate log likelihood
    ll = np.sum(norm(sim, sigma_est).logpdf(obs))
        
    return ll

def log_prior(params):
    """ Calculate log prior.
    """   
    # Extract parameter values
    a_est, b_est, sigma_est = params
    a_min, a_max = -3, 3
    b_min, b_max = -20, -5
    sigma_min, sigma_max = 0, 5

    # If all parameters are within allowed ranges, return a constant 
    # (anything will do - I've used 0 here)
    if ((a_min <= a_est < a_max) and 
        (b_min <= b_est < b_max) and 
        (sigma_min <= sigma_est < sigma_max)):
        return 0
    # Else the parameter set is invalid (probability = 0; log prob = -inf)
    else:
        return -np.inf

def log_posterior(params, x, obs):
    """ Calculate log posterior.
    """  
    # Get log prior prob
    log_pri = log_prior(params)

    # Evaluate log likelihood if necessary
    if np.isfinite(log_pri):
        log_like = log_likelihood(params, x, obs)
        
        # Calculate log posterior
        return log_pri + log_like
    
    else:
        # Log prior is -inf, so log posterior is -inf too
        return -np.inf  
def neg_log_posterior(params, x, obs):
    """ Negative of log posterior.
    """
    return -log_posterior(params, x, obs)

def find_map(init_guess, x, obs):
    """ Find max of posterior.
        init_guess [a, b, sigma]
    """
    # Run optimiser
    param_est = optimize.fmin(neg_log_posterior, 
                              init_guess, 
                              args=(x, obs))  
    return param_est

def run_mcmc(XX, YY, param_guess):
    '''mcmc for linear fitting code'''

    # Run optimiser
    param_est = find_map(param_guess, XX, YY)

    # emcee parameters
    n_dim = 3      # Number of parameters being calibrated
    n_walk = 10    # Number of "walkers"/chains
    n_steps = 3000  # Number of steps per chain
    n_burn = 50   # Length of burn-in to discard

    # Generate starting locations for the chains by adding a small
    # amount of Gaussian noise to optimised MAP
    starting_guesses = [param_est + 1e-4*np.random.randn(n_dim) 
                        for i in range(n_walk)]

    # Prepare to sample. The params are automatically passed to log_posterior
    # as part of n_dim. "args" lists the other params that are also necessary
    sampler = emcee.EnsembleSampler(n_walk, n_dim, log_posterior, 
                                    args=[XX, YY])

    # Run sampler
    pos, prob, state = sampler.run_mcmc(starting_guesses, n_steps)

    # Print some stats. based on run properties
    print( '\n')
    print( 'Average acceptance fraction: ', np.mean(sampler.acceptance_fraction))
    print( 'Autocorrelation time: ', sampler.acor)

    # Discard burn-in
    samples = sampler.chain[:, n_burn:, :].reshape((-1, n_dim))

    # Print estimates and confidence intervals
    mcmc_df = pd.DataFrame(data=samples, columns=['a_mcmc', 'b_mcmc', 'sigma_mcmc'])
    mcmc_des = mcmc_df.describe(percentiles=[0.16, 0.5, 0.84])
    print(mcmc_des)

    # Store output data in lists
    conf = []
    x = np.linspace(XX.min(), XX.max(), 100)
    # Pick parameter sets at random from the converged chains
    for a, b, sigma in samples[np.random.randint(len(samples), size=1000)]:
        # Simulate values
        sim = a*x + b + norm.rvs(loc=0, scale=sigma, size=len(x))
        df = pd.DataFrame(data={'Sim':sim})
        # Add to conf
        conf.append(df)

    # Concatenate results
    conf = pd.concat(conf, axis=1)

    # Get 2.5 and 97.5 percentiles for plotting
    conf = conf.T.describe(percentiles=[0.16, 0.5, 0.84]).T[['16%', '50%', '84%']]

    # Plot predicted
    ax.fill_between(x, np.median(conf['16%'].values), np.median(conf['84%'].values), color='g', alpha=0.3)
    ax.plot(x, np.ones_like(x)*np.median(conf['50%'].values), 'g-', label='Estimated')
    


# STATSMODELS
def run_statsmodels(XX, YY, plot=False, ax=None, extend=True, listitems=False, color='k'):
    x = sm.add_constant(XX)
    model = sm.OLS(YY,x).fit()
    N = len(XX)
    print(model.summary())
    intercept = model.params[0]
    slope = model.params[1]
    intercepterr = model.bse[0]
    slopeerr = model.bse[1]
    stdfit =np.sqrt(model.scale)
    r2= model.rsquared
    if plot:
        pred_ols = model.get_prediction()
        iv_l = pred_ols.summary_frame(alpha=0.318)["obs_ci_lower"]
        iv_u = pred_ols.summary_frame(alpha=0.318)["obs_ci_upper"]
        ax.plot(XX, model.fittedvalues, color=color, alpha=1, zorder=15)
        ind = np.argsort(XX)
        XXs =XX[ind]
        iv_ls, iv_us = iv_l.values[ind], iv_u.values[ind]

        ax.fill_between(XXs, iv_ls, iv_us, 
                       color=color, alpha=0.2, edgecolor='none', zorder=10)
        if extend:
            x2 = np.array([-2.5, XX.min()])
            ax.plot(x2, model.params[1]*x2+model.params[0], color=color, alpha=0.4, zorder=5, linestyle='--')

    if listitems:
        return [N, slope, slopeerr, intercept, intercepterr, stdfit, r2]

    else:
        return N, slope, slopeerr, intercept, intercepterr, stdfit, r2


##### curve_fit no slope ####



# FIT INTERCEPT
def fstar(x, b):
    slope = linfit.at['CASPAR stars', 'slope'] #=2.1715174
    return slope*x+b
def fbd(x,b):
    slope = linfit.at['CASPAR BDs', 'slope'] #=3.1941512
    return slope*x+b
def f(x,b):
    slope = linfit.at['CASPAR all', 'slope'] #=2.0191426
    return slope*x+b
def fall(x,b):
    slope = linfit.at['CASPAR all', 'slope'] #=2.0191426
    return slope*x+b
        
def run_interceptfit(XX, YY, obj, yerr=None):
    if obj=='star':
        B, pcov = curve_fit(fstar, XX, YY)
        M = 1.84
        Berr= np.sqrt(pcov[0][0])
    elif obj=='bd':
        B, pcov = curve_fit(fbd, XX, YY)
        M = 3.08
        Berr= np.sqrt(pcov[0][0])
    elif obj=='all':
        B, pcov = curve_fit(fall, XX, YY)
        M = 1.75
        Berr = np.sqrt(pcov[0][0])
    return M, B, Berr


# LINMIX          
def linmix_fitting(df2, fitname='doLinearFit', xval='logMass', plot=False, ax=None, extend=True, as_list=True, color='k', sig3 = False, avgMdoterr=0.3658954994942055):
   
    if fitname=='doLinearFit':
        
        df = df2.copy()
        y = df['log Mdot'].values
        # tiered upper limit error limits
        # mass > 0.075 => 0.3
        # mass < 0.075 => 0.6

        df.loc[(df['log Mdot err']==0)&(df['log Mass']>np.log10(0.075)), 'log Mdot err'] = 0.3 
        df.loc[(df['log Mdot err']==0)&(df['log Mass']<=np.log10(0.075)), 'log Mdot err'] = 0.6
        
        ysig = df['log Mdot err'].values
        upper_limit_boolarray = df['Upper Limit bool'].values
        upplims = np.logical_not(upper_limit_boolarray) # upper limits

        if xval=='Age':
            x = df['log Age'].values
            df.loc[df['log Age err'].isnull(), 'log Age err'] = 0.16 # median log age err
            xsig = df['log Age err'].values

            idx = np.isfinite(x)
            x, xsig, y, ysig, upper_limit_boolarray = x[idx], xsig[idx], y[idx], ysig[idx], upper_limit_boolarray[idx]
        else:
            x = df['log Mass'].values
            df.loc[df['log Mass err']==0, 'log Mass err'] = 0.2
            xsig = df['log Mass err'].values
        plt.figure()
        plt.errorbar(x,y,xerr=xsig, yerr=ysig,fmt='o')
        plt.show()

        lmcens  = linmix.LinMix(x, y, xsig, ysig, delta=upper_limit_boolarray, K=2, seed=8)
        lmcens.run_mcmc(silent=False)
        print('finished!!')
        if plot:
            xs = np.linspace(x.min()-1, x.max()+1, 100)
            for i in range(0, len(lmcens.chain), 25):
                ys = lmcens.chain[i]['alpha'] + xs * lmcens.chain[i]['beta']
                ax.plot(xs, ys, color=color, alpha=0.02, zorder=15)
            ax.plot(xs, np.nanmean(lmcens.chain['alpha']) + xs*np.nanmean(lmcens.chain['beta']), 
                    color=color, alpha=1, zorder=15, linewidth=2)

            if extend:
                x2 = np.array([-2.5, x.min()])
                ax.plot(x2, np.nanmean(lmcens.chain['alpha']) + x2*np.nanmean(lmcens.chain['beta']), color=color, alpha=0.4, zorder=5, linestyle='--')


        data = lmcens.chain

        slopevals = np.percentile(data['beta'], [16, 50, 84])
        slopeerrs = np.diff(slopevals)
        interceptvals = np.percentile(data['alpha'], [16, 50, 84])
        intercepterrs = np.diff(interceptvals)

        N = len(x)
        slope = slopevals[1]
        slopeerrlow = slopeerrs[0]
        slopeerrupp = slopeerrs[1]

        intercept = interceptvals[1]
        intercepterrlow = intercepterrs[0]
        intercepterrupp = intercepterrs[1]

        stdfit = np.median(data['sigsqr'])
        r2 = np.mean(data['corr'])
#
    else: #method=<table row name>
        df = df2.copy()

        N, slope, slopeerrlow, slopeerrupp, intercept, intercepterrlow, intercepterrupp, stdfit, r2 = linfit.loc[fitname].values
        
#        slope = linfit.at[fitname, 'slope']  
#        intercept = linfit.at[fitname, 'intercept']  
#        stdfit = linfit.at[fitname, 'stdfit'] 

        if plot:
            if xval=='Age':
                xvals = df['log Age'].values

                xx = np.linspace(np.nanmin(xvals)-1, np.nanmax(xvals)+1, 100)
            else:
                xvals = df['log Mass'].values
                xx = np.linspace(np.nanmin(xvals), np.nanmax(xvals), 100)

            yy = slope*xx + intercept
            ax.plot(xx, yy, color=color, alpha=1, zorder=150)
            ax.fill_between(xx, yy-stdfit, yy+stdfit, 
                           color=color, alpha=0.2, edgecolor='none', zorder=100)
            if extend:
                x2 = np.array([-2.5, xx.min()])
                ax.plot(x2, slope*x2+intercept, color=color, alpha=0.4, zorder=50, linestyle='--')
                if xx.max() < -1.2:
                    x3 = np.array([xx.max(), -1.12])
                    ax.plot(x3, slope*x3+intercept, color=color, alpha=0.4, zorder=50, linestyle='--')            
    if as_list:
        return [N, slope, slopeerrlow, slopeerrupp, intercept, intercepterrlow, intercepterrupp, stdfit, r2]
    else:
        return N, slope, slopeerrlow, slopeerrupp, intercept, intercepterrlow, intercepterrupp, stdfit, r2
    
    
# ODR  
def fODR_all(B, x):
    '''Linear function y = m*x + b'''
    slope = linfit.at['CASPAR all', 'slope'] #=2.0191426 
    return slope*x + B[0]
def fODR_star(B, x):
    '''Linear function y = m*x + b'''
    slope = linfit.at['CASPAR stars', 'slope'] #=2.1715174
    return slope*x + B[0]
def fODR_bd(B, x):
    '''Linear function y = m*x + b'''
    slope = linfit.at['CASPAR BDs', 'slope'] #=3.1941512
    return slope*x + B[0] 
def fODR_xy(B, x):
    '''Linear function y = m*x + b'''
    return B[1]*x + B[0]

def run_odrI(x_var, y_var, x_err, y_err, obj, job=0, avgMdoterr= 0.3658954994942055):

    if not isinstance(x_var, list):
        x_var, y_var = x_var.values, y_var.values
        x_err, y_err = x_err.values, y_err.values
    IND = np.where(y_err==0)
    
    y_err[IND]=avgMdoterr

    data = odr.RealData(x_var,y_var,sx = x_err, sy = y_err)
    if obj == 'star':
        quad_model = odr.Model(fODR_star)
    elif obj == 'bd':
        quad_model = odr.Model(fODR_bd)
    elif obj == 'xy':
        quad_model=odr.Model(fODR_xy)
    else:
        quad_model = odr.Model(fODR_all)
    
    if obj != 'xy':
        Odr = odr.ODR(data, quad_model,beta0 =[-8])
    else:
        print('xy')
        Odr = odr.ODR(data, quad_model,beta0 =[-8,2])
    Odr.set_job(fit_type=job)
    out =Odr.run()

    #out.sd_beta is a list which contains the values of the errors of the parameters 
    return out.beta, out.sd_beta
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            