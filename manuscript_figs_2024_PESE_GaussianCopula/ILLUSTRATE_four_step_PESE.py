'''
    SCRIPT TO DEMONSTRATE THE PROBIT SPACE RESAMPLING STRATEGY
    ---------------------------------------------------------- 
'''



''' 
    LOAD USEFUL LIBRARIES
'''
import numpy as np
import scipy.stats as stats
import random

from matplotlib import use 
use('agg')
import matplotlib.pyplot as plt





'''
    FUNCTION TO GENERATE GAUSSIAN RESAMPLING COEFFICIENTS 

    Based on M.-Y. Chan et al 2020 Monthly Weather Review paper on Bi-Gaussian 
    EnKFs

    Inputs:
    1) N        -- Original ensemble size
    2) M        -- Number of additional members

    Note that M must be greater than N!
'''
def compute_gaussian_resampling_coefficients( N, M ):


    # Generating matrix of noise (Appendix B Step 1)
    W = np.random.normal( size=[N,M] )

    # Appendix B Step 2
    W = np.matrix( (W.T - np.mean(W, axis=1)).T )

    

    # Appendix B Step 3 and 4
    C_W = W * W.T
    inv_L_W = np.linalg.inv(np.linalg.cholesky( C_W ))
    inv_L_W = np.matrix( inv_L_W )


    # Appendix B Step 5
    k = np.sqrt( 
                ( M+N-1. )/(N-1.)
                )
    C_E = np.eye(N)* M/(N-1.)
    C_E -= (k-1)**2/M
    L_E = np.matrix(np.linalg.cholesky( C_E ))
    

    # Appendix B Step 6
    E_prime = L_E * inv_L_W * W


    # Appendix B Step 7
    E = E_prime + (k-1)/M

    return np.array( E )


















# FRAMES FOR ANIMATIONS
nFrames = 2
x1lims = [-3,0.5]
x2lims = [-0.5,4]





'''
    USEFUL PARAMETERS CONTROLLING THE SYSTEMATIC TESTS
'''


np.random.seed(0)
random.seed(0)

# Basic settings
n_ens_src = 50               # Size of small ens
n_ens_big = 110             # Size of expanded ens



# Probit-space covariance for the model ensemble
probit_prior_cov = np.matrix(
                    [
                        [ 1.0, 0.8 ],
                        [ 0.8, 1.0 ]
                    ]
                )
sqrt_probit_prior_cov = np.linalg.cholesky( probit_prior_cov ) 


# Set up true prior pdf (used to generate ensemble members)
true_prior_pdf = {}
true_prior_pdf[0] = stats.skewnorm( -7, scale=1. )
#stats.norm(  loc=-1, scale=1. )
true_prior_pdf[1] = stats.gamma( 2, scale=0.5 )


# Proposal pdf for prior
proposal_prior_pdf = {}
proposal_prior_pdf[0] = stats.skewnorm
proposal_prior_pdf[1] = stats.gamma




































'''
    BINS USED TO PLOT JOINT PDFS
'''
# Bins to generate joint pdf
x2e = np.linspace( -1,4, 101)
x1e = np.linspace( -3,1, 102)
p1e = np.linspace(-4,4,101 )
p2e = np.linspace(-4,4,102 )










































'''
    INITIALIZE DICTIONARY TO HOLD PDFS USED TO TRANSFORM MODEL VARIABLES
'''
sample_prior_pdf = {}
sample_prior_pdf['expanded_ens'] = {}
sample_prior_pdf['small_ens'] = {}







































'''
    GENERATE SMALL ENSEMBLE
'''

# Generate correlated Gaussian noise
probit_noise = np.random.normal( size= [n_ens_src,2] )
probit_noise = (probit_noise - np.mean( probit_noise, axis=0)).T
noise_cov = np.cov( probit_noise, ddof=1 )
cov_transform = (
    np.matrix( sqrt_probit_prior_cov )
    * np.matrix( np.linalg.inv( np.linalg.cholesky(noise_cov) ) )
)
probit_noise = np.array( cov_transform * probit_noise )

# Map correlated Gaussian noise to targetted marginal distributions
quantile_noise = stats.norm.cdf( probit_noise )
prior_ens_ORIGINAL = probit_noise * 0.
prior_ens_ORIGINAL[0,:] = true_prior_pdf[0].ppf( quantile_noise[0,:] )
prior_ens_ORIGINAL[1,:] = true_prior_pdf[1].ppf( quantile_noise[1,:] )


            

            


        































'''
    EXPAND ENSEMBLE VIA CAC2020 RESAMPLING IN PROBIT SPACE
'''
small_ens_probits = prior_ens_ORIGINAL*0.



# First, fit the proposed pdfs to the small ens and map to probit space
for i in range(2):

    # Fit data
    fit_params = proposal_prior_pdf[i].fit( prior_ens_ORIGINAL[i,:] * 1. )
    
    # Generate fitted pdf
    if ( proposal_prior_pdf[i] == stats.norm ):
        prior_std = np.std( prior_ens_ORIGINAL[i,:], ddof=1 )
        prior_avg = np.mean( prior_ens_ORIGINAL[i,:])
        sample_prior_pdf['expanded_ens'][i] = stats.norm( loc=prior_avg, scale=prior_std )

    elif( len( fit_params ) == 2):
        sample_prior_pdf['expanded_ens'][i] = proposal_prior_pdf[i]( 
            loc = fit_params[0], scale = fit_params[1]
        )
    elif( len( fit_params) == 3 ):
        sample_prior_pdf['expanded_ens'][i] = proposal_prior_pdf[i]( 
            fit_params[0], loc = fit_params[1], scale = fit_params[2]
        )

    # Now convert to probits
    small_ens_probits[i,:] = (
        stats.norm.ppf(
            sample_prior_pdf['expanded_ens'][i].cdf( prior_ens_ORIGINAL[i,:]*1. )
        )
    )
    small_ens_probits[i,:] -= np.mean(small_ens_probits[i,:])
    small_ens_probits[i,:] /= np.std( small_ens_probits[i,:], ddof=1 )

    
# --- Finished fitting pdfs and mapping to probit space

    


# Use CAC20 trick to generate new ens members
F_CAC20 = compute_gaussian_resampling_coefficients( n_ens_src, n_ens_big - n_ens_src )
small_ens_perts = (small_ens_probits.T - np.mean(small_ens_probits.T, axis=0)).T
ens_CAC20_probits = np.zeros([2, n_ens_big] )
ens_CAC20_probits[:,:n_ens_src] = small_ens_perts * 1.
ens_CAC20_probits[:,n_ens_src:] = np.array(np.matrix( small_ens_perts ) * np.matrix( F_CAC20 ))
ens_CAC20_probits[:,:] = (ens_CAC20_probits.T + np.mean(small_ens_probits.T, axis=0)).T



# Map expanded ens from probit space to native space
ens_CAC20 = ens_CAC20_probits*0.
for i in range(2):
    ens_CAC20[i,:] = (
        sample_prior_pdf['expanded_ens'][i].ppf(
            stats.norm.cdf( ens_CAC20_probits[i,:])
        )
    )





























































'''
    FIT MARGINAL PDFS TO DATA
'''

# Init figure
fig,axs = plt.subplots(nrows=2, ncols=2, figsize=(4,4))
axs[1,0].remove()

axs[0,1].set_xlim( x1lims)
axs[0,1].set_ylim( x2lims)
axs[0,1].set_position([0.4, 0.4, 0.5, 0.5] )
axs[0,1].set_title('2D Model Space')

axs[1,1].set_xlim( x1lims)
axs[1,1].set_position([0.4,0.15,0.5,0.15] )
axs[1,1].set_xlabel('Model variable 1')
axs[1,1].set_ylabel('PDF')

axs[0,0].set_ylim( x2lims)
axs[0,0].set_position([0.15, 0.4, 0.15, 0.5] )
axs[0,0].set_ylabel('Model variable 2')
axs[0,0].set_xlabel('PDF')



# Plot out 2d ens members
axs[0,1].scatter( prior_ens_ORIGINAL[0,:], prior_ens_ORIGINAL[1,:],
             marker='x', color='r', s=20, zorder=100)


# Plot out x1 marginals
x1_range = np.linspace(x1lims[0],x1lims[1],101)
axs[1,1].scatter( prior_ens_ORIGINAL[0,:], prior_ens_ORIGINAL[0,:]*0., 
                  marker='x', color='r', s=20 )
axs[1,1].plot( x1_range, sample_prior_pdf['expanded_ens'][0].pdf( x1_range ), '-r', alpha=0 )



# Plot out x2 marginals
x2_range = np.linspace(x2lims[0], x2lims[1],101)
axs[0,0].scatter( prior_ens_ORIGINAL[1,:]*0., prior_ens_ORIGINAL[1,:], 
                  marker='x', color='r', s=20 )
#axs[0,0].set_xlim( [-0.05,0.6])
axs[0,0].plot( sample_prior_pdf['expanded_ens'][1].pdf( x2_range ), x2_range, '-r',alpha=0)


plt.savefig('ILLUSTRATE_four_step_PESE/ILLUSTRATE_PESE_stage0.svg')





# Plot fitted PDF for x1
axs[1,1].plot( x1_range, sample_prior_pdf['expanded_ens'][0].pdf( x1_range ), '-r',)


# Plot fitted PDF for x2
axs[0,0].plot( sample_prior_pdf['expanded_ens'][1].pdf( x2_range ), x2_range, '-r')
plt.savefig('ILLUSTRATE_four_step_PESE/ILLUSTRATE_PESE_stage1.svg')

plt.close()











'''
    PLOT PPI-TRANSFORMED ENSEMBLE
'''

# Init figure
fig,axs = plt.subplots(nrows=2, ncols=2, figsize=(4,4))
axs[1,0].remove()

axs[0,1].set_xlim( x1lims)
axs[0,1].set_ylim( x2lims)
axs[0,1].set_position([0.4, 0.4, 0.5, 0.5] )
axs[0,1].set_title('2D Probit Space')

axs[1,1].set_xlim( x1lims)
axs[1,1].set_position([0.4,0.15,0.5,0.15] )
axs[1,1].set_xlabel('Probit variable 1')
axs[1,1].set_ylabel('PDF')

axs[0,0].set_ylim( x2lims)
axs[0,0].set_position([0.15, 0.4, 0.15, 0.5] )
axs[0,0].set_ylabel('Probit variable 2')
axs[0,0].set_xlabel('PDF')


# Generate plot
axs[0,1].scatter( 
    small_ens_probits[0],
    small_ens_probits[1],
    marker='x', color='r', s=20, zorder=100)
axs[0,1].set_ylim( [-3, 3] )
axs[0,1].set_xlim( [-3,3] )



# Plot out x2 data points
axs[0,0].scatter( 
    prior_ens_ORIGINAL[0,:]*0., 
    small_ens_probits[1],
    marker='x', color='r', s=20 )
axs[0,0].set_ylim( [-3,3] )

# Plot out PDF
xc = np.linspace(-3,3,1000)
axs[0,0].plot( stats.norm.pdf(xc), xc ,'-r')



# Plot out p1 marginals
axs[1,1].scatter( small_ens_probits[0], small_ens_probits[0]*0., 
                    marker='x', color='r', s=20 )
axs[1,1].plot( xc, stats.norm.pdf(xc), '-r' )
axs[1,1].set_xlim( [-3,3])


plt.savefig('ILLUSTRATE_four_step_PESE/ILLUSTRATE_PESE_stage2.svg')








'''
    ADD NEW PROBIT MEMEBRS
'''



# Stick on new ensemble members
# -----------------------------
# Plot out 2d ens members
axs[0,1].scatter( ens_CAC20_probits[0,n_ens_src:], ens_CAC20_probits[1,n_ens_src:],
                marker='o', color='dodgerblue', s=5, zorder=1000)
axs[1,1].scatter( ens_CAC20_probits[0,n_ens_src:], ens_CAC20_probits[1,n_ens_src:]*0,
                marker='o', color='dodgerblue', s=5, zorder=1000)
axs[0,0].scatter( ens_CAC20_probits[0,n_ens_src:]*0., ens_CAC20_probits[1,n_ens_src:],
                marker='o', color='dodgerblue', s=5, zorder=1000)
plt.savefig('ILLUSTRATE_four_step_PESE/ILLUSTRATE_PESE_stage3.svg')
plt.close()





























'''
    INVERT PPI TRANSFORM
'''
iframe=nFrames



# Init figure
fig,axs = plt.subplots(nrows=2, ncols=2, figsize=(4,4))
axs[1,0].remove()

axs[0,1].set_xlim( x1lims)
axs[0,1].set_ylim( x2lims)
axs[0,1].set_position([0.4, 0.4, 0.5, 0.5] )
axs[0,1].set_title('2D Model Space')

axs[1,1].set_xlim( x1lims)
axs[1,1].set_position([0.4,0.15,0.5,0.15] )
axs[1,1].set_xlabel('Model variable 1')
axs[1,1].set_ylabel('PDF')

axs[0,0].set_ylim( x2lims)
axs[0,0].set_position([0.15, 0.4, 0.15, 0.5] )
axs[0,0].set_ylabel('Model variable 2')
axs[0,0].set_xlabel('PDF')

# Generate 2d plot
axs[0,1].scatter( 
    prior_ens_ORIGINAL[0,:],
    prior_ens_ORIGINAL[1,:n_ens_src],
    marker='x', color='r', s=20, zorder=100)
axs[0,1].scatter( 
    ens_CAC20[0,n_ens_src:],
    +ens_CAC20[1,n_ens_src:],
    marker='o', color='dodgerblue', s=5, zorder=1000)

axs[0,1].set_ylim( x2lims )
axs[0,1].set_xlim( x1lims)


# Plot out x2 data points & marginal
axs[0,0].scatter( 
    small_ens_probits[1,:]*0,
    prior_ens_ORIGINAL[1,:n_ens_src],        
    marker='x', color='r', s=20, zorder=100)
axs[0,0].scatter( 
    ens_CAC20_probits[1,n_ens_src:]*0,
    ens_CAC20[1,n_ens_src:],
    marker='o', color='dodgerblue', s=5, zorder=1000)

x2c = np.linspace(x2lims[0], x2lims[1], 1000)
axs[0,0].plot( sample_prior_pdf['expanded_ens'][1].pdf(x2c),x2c,'-r')
axs[0,0].set_ylim( x2lims )



# Plot out x1 data & marginal
x1c = np.linspace(x1lims[0], x1lims[1], 1000)
axs[1,1].scatter( ens_CAC20[0,n_ens_src:], ens_CAC20[0,n_ens_src:]*0., 
                marker='o', color='dodgerblue', s=5, zorder=1000 )
axs[1,1].scatter( prior_ens_ORIGINAL[0,:], prior_ens_ORIGINAL[0,:]*0., 
                marker='x', color='r', s=20 )
axs[1,1].plot( x1c, sample_prior_pdf['expanded_ens'][0].pdf(x1c), '-r' )
axs[1,1].set_xlim(x1lims)


plt.savefig('ILLUSTRATE_four_step_PESE/ILLUSTRATE_PESE_stage4.svg')
plt.close()
    
