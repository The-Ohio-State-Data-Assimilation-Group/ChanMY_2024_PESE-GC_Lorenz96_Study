'''
    DEMO: PESE WITH DIFFERENT KINDS OF BIVARIATE RELATIONSHIPS
    ----------------------------------------------------------
'''

import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt

from scipy.stats import norm, gamma
import pylib_rank_histogram_filter as rhf
from copy import deepcopy
from os.path import exists
from math import sqrt, erf, pi, isinf

from sklearn.mixture import GaussianMixture
from numba import njit


'''
    Generate colorbars for the figure
'''
fig, axs = plt.subplots( nrows=1, ncols=2, figsize=(7,3))

xmesh, ymesh = np.meshgrid( np.linspace(0,1,100), np.linspace(0,2,101) )
zmesh = np.sqrt(xmesh**2 + ymesh**2)

cnf1 = axs[0].contourf(xmesh, ymesh, zmesh, np.linspace(0.1,0.9,5), extend='max', cmap='Reds')
cnf2 = axs[1].contourf(xmesh, ymesh, zmesh, np.linspace(0.1,0.9,5), extend='max', cmap='Blues')

fig.subplots_adjust(bottom=0.3)

cbar1 = fig.add_axes( [0.05,0.1,0.4,0.05])
fig.colorbar( cnf1, cax = cbar1, orientation='horizontal')
cbar2 = fig.add_axes( [0.55,0.1,0.4,0.05])
fig.colorbar( cnf2, cax = cbar2, orientation='horizontal')
plt.savefig('DEMO_PESE_GMM_with_diff_nonlinearities/DEMO_PESE_GMM_with_diff_nonlinearities_COLORBARS.svg')
plt.close()





'''
    FUNCTION TO GENERATE RESAMPLING COEFFICIENTS WITH DIFFERENT DRAWING PROCEDURES
'''

# Function to generate resampling coefficients based on pre-calculated random draws (W)
def compute_resampling_coefficients( N, M ):

    # Draw samples
    W = np.random.normal(size=[N,M])

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
    C_E -= (k-1)**2/N
    L_E = np.matrix(np.linalg.cholesky( C_E ))
    

    # Appendix B Step 6
    E_prime = L_E * inv_L_W * W


    # Appendix B Step 7
    E = E_prime + (k-1)/M

    # print( np.matrix(E_prime) * np.matrix(E_prime).T)
    # print( np.sum(E, axis=1), k)

    return np.array( E )








'''
    APPROXIMATE QUANTILE FUNCTION OF GAUSSIAN MIXTURE MODEL
'''


# Cumulative density function of 1D gaussian mixture model 
@njit
def cum_density_func_GMM1d( kernel_weights, kernel_means, kernel_sigmas, data1d ):

    # Determine number of kernels
    nKernels = kernel_weights.shape[0]

    # Determine number of data points
    nData = data1d.shape[0]

    # Init array to hold CDF values
    cdf_vals = data1d *0.

    # Compute CDF by looping all kernels
    for kk in range(nKernels):
        for dd in range(nData):
            cdf_vals[dd] += (
                kernel_weights[kk] * 0.5
                * (1.
                   + erf(
                        (data1d[dd]-kernel_means[kk])
                        / (kernel_sigmas[kk]*sqrt(2))
                    )
                )
            )
    # --- End of loop over kernels

    return cdf_vals





# Quantile function of 1D gaussian mixture model 
def quantile_func_GMM1d( weights, means, std_devs, input_cdf_vals ):

    ref_xvals = np.linspace(-10,10,10000)
    ref_cdf_vals = cum_density_func_GMM1d( 
        np.array(weights), np.array(means), np.array(std_devs), np.array(ref_xvals) 
    )
    
    # Approximate evaluation of quantile function
    return np.interp( input_cdf_vals, ref_cdf_vals, ref_xvals )





'''
    GENERAL INVERSION ALGORITHM
'''
# Approximate inversion of numerically-evaluated function
def approx_function_inv( ref_xvals, ref_yvals, target_yvals ):
    return np.interp( target_yvals, ref_yvals, ref_xvals )










'''
    APPLY PESE WITH GAUSSIAN MIXTURE MODEL MARGINALS AND PLOT OUT BEHAVIOR
'''
def visualize_pese_in_native_and_probit_space( truth_samples, N, E, xrange, yrange, n_gmm_components=1 ):

    # Generate small initial ensemble
    small_samples = deepcopy( truth_samples[:,:N] )

    # NUmber of variables considered
    nVariables = truth_samples.shape[0]

    # Init array to hold all probits and virtual members
    truth_probits = np.zeros( truth_samples.shape )
    small_probits = np.zeros( small_samples.shape)
    virtual_samples = np.zeros( [nVariables,  E.shape[1]])
    virtual_probits = np.zeros( [nVariables,  E.shape[1]])

    # Init dictionary to hold GMM parameters
    gmm_params = {}
    for sampletype in ['truth','small']:
        gmm_params[sampletype] = {}
        for param in ['weights','sigmas','ctrs']:
            gmm_params[sampletype][param] = np.zeros( [nVariables, n_gmm_components])


    # Probit space Ensemble Size Expansion
    # ------------------------------------
    for i in range( nVariables ):

        # For each variable, fit a GMM distribution to small samples
        # ----------------------------------------------------------

        # Run scikit's GMM fitting code on small samples
        tmp_gmm_func = GaussianMixture(n_components=n_gmm_components)
        tmp_gmm_func.fit( (small_samples[i:i+1,:]).T )
        gmm_params['small']['weights'][i]  = tmp_gmm_func.weights_[:]
        gmm_params['small']['sigmas'][i]   = np.sqrt(tmp_gmm_func.covariances_[:,0,0])
        gmm_params['small']['ctrs'][i]     = tmp_gmm_func.means_[:,0]


        # Map to probit space
        # --------------------
        small_probits[i] = norm.ppf(
            cum_density_func_GMM1d( 
                gmm_params['small']['weights'][i], gmm_params['small']['ctrs'][i], gmm_params['small']['sigmas'][i],
                small_samples[i,:]
            )
        )


        # Generate virtual probits
        # -------------------------
        small_probits[i] -= np.mean(small_probits[i])
        small_probits[i] /= np.std(small_probits[i], ddof=1)
        virtual_probits[i] = np.array(np.matrix(small_probits[i:(i+1),:]) * np.matrix(E))


        # Map virtual probits to real space
        # ----------------------------------
        virtual_samples[i] = quantile_func_GMM1d( 
            gmm_params['small']['weights'][i], gmm_params['small']['ctrs'][i], gmm_params['small']['sigmas'][i],
            norm.cdf( virtual_probits[i]) 
        )

    # ---- End of loop over variables



    # Construct probits for truth samples via RH approach
    truth_probits = truth_samples *0.
    for i in range(truth_samples.shape[0]):
        truth_probits[i] =  norm.ppf(
            ( np.argsort( np.argsort(truth_samples[i]) ) + 1.) 
            / ( truth_samples.shape[1] + 1. )
        )




    # FIGURE GENERATION SECTION
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(4,5))
    axs = [val for sub in axs for val in sub]

    # Plot histogram of truth samples
    hist, xe, ye = np.histogram2d( truth_samples[0], truth_samples[1],
                                   bins = (xrange, yrange), density=True )
    xc = (xe[1:]+xe[:-1])/2
    yc = (ye[1:]+ye[:-1])/2
    hist /= hist.max()
    cnf1 = axs[0].contourf( 
        xc, yc, hist.T, np.linspace(0.1,0.9,5), cmap='Reds', extend='max' 
    )
    axs[0].set_xlabel( 'Model variable 1')
    axs[0].set_ylabel( 'Model variable 2')
    axs[0].set_title('True PDF')

    # Plot histogram of virtual samples
    hist, xe, ye = np.histogram2d( virtual_samples[0], virtual_samples[1],
                                   bins = (xrange, yrange), density=True )
    xc = (xe[1:]+xe[:-1])/2
    yc = (ye[1:]+ye[:-1])/2
    hist /= hist.max()
    axs[1].contourf( 
        xc, yc, hist.T, np.linspace(0.1,0.9,5), cmap='Reds', extend='max' 
    )
    axs[1].yaxis.set_ticklabels([])
    axs[1].set_xlabel( 'Model variable 1')
    axs[1].set_title('Virt Mems\' PDF')
    
        

    # Plot histogram of truth probits
    p1range = np.linspace(-3.3,3.3,30); p2range= np.linspace(-3.3,3.3,31)
    hist, xe, ye = np.histogram2d( truth_probits[0], truth_probits[1],
                                   bins = (p1range, p2range), density=True )
    xc = (xe[1:]+xe[:-1])/2
    yc = (ye[1:]+ye[:-1])/2
    hist /= hist.max()
    cnf2 = axs[2].contourf( 
        xc, yc, hist.T, np.linspace(0.1,0.9,5), cmap='Blues', extend='max' 
    )
    axs[2].set_xlabel( 'Probit variable 1')
    axs[2].set_ylabel( 'Probit variable 2')
    axs[2].set_title('True PDF')

    # Plot histogram of virtual probits
    hist, xe, ye = np.histogram2d( virtual_probits[0], virtual_probits[1],
                                   bins = (p1range, p2range), density=True )
    hist /= hist.max()
    axs[3].contourf( 
        xc, yc, hist.T, np.linspace(0.1,0.9,5), cmap='Blues', extend='max' 
    )
    axs[3].yaxis.set_ticklabels([])
    axs[3].set_xlabel( 'Probit variable 1')
    axs[3].set_title('Virt Mems\' PDF')    



    fig.subplots_adjust(left=0.2, right=0.98, bottom=0.2, wspace=0.1, hspace=0.6, top=0.93)

    # Return figure
    return fig, axs















'''
    CONTROLLING CONSTANTS
'''
np.random.seed(0)
N= 100
M= 1000000 #10000000
nTruth = M*10

# Generate E matrix
if exists( 'E_matrix.npy' ):
    E_matrix = np.load( 'E_matrix.npy')
    if ( (E_matrix.shape[0] != N) or  (E_matrix.shape[1] != M) ):
        E_matrix = compute_resampling_coefficients(N, M)
        np.save('E_matrix.npy', E_matrix)
else: 
    E_matrix = compute_resampling_coefficients(N, M)
    np.save('E_matrix.npy', E_matrix)



# Warm up a function
cum_density_func_GMM1d( np.ones(2)*0.5, np.linspace(0,1,2), np.arange(2)+1. , np.random.normal(size=30) )
print('warmup done')






'''
    GAUSSIAN SCENARIO
'''

# Generate truth samples
np.random.seed(0)
cov = np.matrix( [ [2.0, -1.4],  [-1.4, 2.0] ])
chol_cov = np.linalg.cholesky( cov )
truth_samples = np.array(
    np.matrix(chol_cov) * np.matrix(np.random.normal(size=[2,nTruth]))
)
xrange = np.linspace(-5.5,5.5,30)
yrange = np.linspace(-5.5,5.5,31)


# Plot results
gauss_fig, gauss_axs = visualize_pese_in_native_and_probit_space( truth_samples, N, E_matrix, xrange, yrange, n_gmm_components=1 )
plt.savefig('DEMO_PESE_GMM_with_diff_nonlinearities/DEMO_PESE_GMM_with_diff_nonlinearities_GaussianLinear.svg')
plt.close()




'''
    GAUSSIAN COPULA SCENARIO
'''

# Generate truth samples
np.random.seed(0)
cov = np.matrix( [ [1.0, -0.7],  [-0.7, 1.0] ])
chol_cov = np.linalg.cholesky( cov )
truth_samples = np.array(
    np.matrix(chol_cov) * np.matrix(np.random.normal(size=[2,nTruth]))
)
xrange = np.linspace(-5.5,5.5,30)
yrange = np.linspace(-5.5,5.5,31)


# PDF of x1 is gamma, PDF of x2 is bi-modal
truth_samples[0] = quantile_func_GMM1d(  [0.5,0.5], [-1.,3.], [2,1.], norm.cdf(truth_samples[0]) )
truth_samples[1] = quantile_func_GMM1d(  [0.5,0.5], [2,-2], [1.,2], norm.cdf(truth_samples[1]) )

# truth_line_x = gamma(2, loc=1, scale=2.).ppf( norm.cdf( probit_line_x ))
# truth_line_y = quantile_func_GMM1d(  [0.5,0.5], [-1.5,1.5], [1,0.5], norm.cdf(probit_line_y) )



# Plot results
gauss_fig, gauss_axs = visualize_pese_in_native_and_probit_space( truth_samples, N, E_matrix, xrange, yrange, n_gmm_components=2 )
plt.savefig('DEMO_PESE_GMM_with_diff_nonlinearities/DEMO_PESE_GMM_with_diff_nonlinearities_GaussianCopula.svg')
plt.close()








'''
    BIJECTIVE BI-GAUSSIAN SCENARIO
'''

# Defining parameters of reference PDF
cov1 = np.matrix( [ [2.0, -0.5],  [-0.5, 0.5] ])
chol_cov1 = np.linalg.cholesky( cov1 )
ctr1 = np.array( [-1, 2])
weight1 = 0.5

cov2 = np.matrix( [ [0.5, -0.5],  [-0.5, 2.] ])
chol_cov2 = np.linalg.cholesky( cov2 )
ctr2 = np.array( [3, -2])
weight2 = 1.-weight1


# Generate truth samples
cluster1_size = int( np.round(weight1*nTruth) )
cluster2_size = int( nTruth - cluster1_size )
truth_samples[:,:cluster1_size] = (np.array(
    np.matrix(chol_cov1) * np.matrix( (np.random.normal(size=[cluster1_size,2])).T )
).T + ctr1).T
truth_samples[:,cluster1_size:] = ( np.array(
    np.matrix(chol_cov2) * np.matrix( (np.random.normal(size=[cluster2_size,2])).T )
).T + ctr2).T
reorder = np.arange(nTruth)
np.random.shuffle( reorder )
truth_samples[:,:] = truth_samples[:,reorder]


# Ranges for this problem
xrange = np.linspace(-5.5,5.5,30)
yrange = np.linspace(-5.5,5.5,31)


# Plot results
bigauss_fig, bigauss_axs = visualize_pese_in_native_and_probit_space( truth_samples, N, E_matrix, xrange, yrange, n_gmm_components=2  )
plt.savefig('DEMO_PESE_GMM_with_diff_nonlinearities/DEMO_PESE_GMM_with_diff_nonlinearities_BiGauss_BiJective.svg')
plt.close()







'''
    NON-BIJECTIVE BI-GAUSSIAN SCENARIO
'''

# Defining parameters of reference PDF
cov1 = np.matrix( [ [2.0, -0.5],  [-0.5, 0.5] ])
chol_cov1 = np.linalg.cholesky( cov1 )
ctr1 = np.array( [0, -1])
weight1 = 0.5

cov2 = np.matrix( [ [0.5, -0.5],  [-0.5, 2.] ])
chol_cov2 = np.linalg.cholesky( cov2 )
ctr2 = np.array( [2.5, 0])
weight2 = 1.-weight1


# Generate truth samples
cluster1_size = int( np.round(weight1*nTruth) )
cluster2_size = int( nTruth - cluster1_size )
truth_samples[:,:cluster1_size] = (np.array(
    np.matrix(chol_cov1) * np.matrix( (np.random.normal(size=[cluster1_size,2])).T )
).T + ctr1).T
truth_samples[:,cluster1_size:] = ( np.array(
    np.matrix(chol_cov2) * np.matrix( (np.random.normal(size=[cluster2_size,2])).T )
).T + ctr2).T
reorder = np.arange(nTruth)
np.random.shuffle( reorder )
truth_samples[:,:] = truth_samples[:,reorder]


# Ranges for this problem
xrange = np.linspace(-5.5,5.5,30)
yrange = np.linspace(-5.5,5.5,31)


# Plot results
bigauss_fig, bigauss_axs = visualize_pese_in_native_and_probit_space( truth_samples, N, E_matrix, xrange, yrange, n_gmm_components=2  )
plt.savefig('DEMO_PESE_GMM_with_diff_nonlinearities/DEMO_PESE_GMM_with_diff_nonlinearities_BiGauss_NON_BiJective.svg')
plt.close()
