'''
    DEMO: USER-INPUTTED INFORMATION CAN IMPROVE PRIOR STATISTICS
    -------------------------------------------------------------

    Bi-variate case: 1 model variable and 1 observable quantity (i.e., x and h(x))
    
    Model variable prior population (i.e., reference) pdf is skewed normal (mean=-1, sigma=2, shape=5)
    Observation operator: h(x) = sign(x) sqrt( abs(x) )
    Observable prior reference pdf is estimated using 1,000,000 Monte Carlo samples drawn from reference x pdf.

    Forecast ensemble size is a command-line input

    Will use empirical distribution functions to measure statistics of prior ensemble and expanded ensemble.
    Will use 1-Wasserstein metric to measure statistical distance between prior/expanded ensemble and the reference distribution

'''

import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt

from scipy.stats import skewnorm, norm, gamma

from copy import deepcopy

from time import time

from sys import argv



'''
    USER INPUTS PESE MARGINALS
'''
pese_pdf_name = argv[1]

small_ens_size= int( argv[2] )

# Reference distribution: skewed normal
ref_dist = skewnorm( 5, scale=2, loc=-1 )











'''
    PRIMARILY USED FUNCTION: RUNNING PESE AND EXAMINING WHAT HAPPENS TO EMPIRICAL DISTRIBUTION FUNCS. AND X-Y RELATIONSHIPS
'''
def single_test_pese_gaussian_copula( ref_x_coord, ref_x_cdf, ref_y_coord, ref_y_cdf,
                                      sorted_small_ens_x, obs_operator,
                                      expanded_ens_size=1000, pese_marginal_name='rank histogram with flat tails'):
    

    # APPLY PESE-GC TO GENERATE EXPANDED ENSEMBLE 
    # -------------------------------------------
    small_ens_size = sorted_small_ens_x.shape[0]
    
    # Apply PPI transform
    if (pese_marginal_name == 'rank histogram with flat tails'):
        sorted_small_ens_x_probits = PPI_rank_histogram_with_flat_tails( sorted_small_ens_x )
    elif (pese_marginal_name == 'normal'):
        small_avg = np.mean( sorted_small_ens_x )
        small_std = np.std(sorted_small_ens_x, ddof=1)
        sorted_small_ens_x_probits = (sorted_small_ens_x-small_avg)/small_std
    elif (pese_marginal_name == 'gamma'):
        fitted_dist_params = gamma.fit( sorted_small_ens_x )
        fitted_dist = gamma( fitted_dist_params[0], loc=fitted_dist_params[1], scale=fitted_dist_params[2])
        sorted_small_ens_x_probits = norm.ppf( fitted_dist.cdf(sorted_small_ens_x) )

    

    # Adjust sorted probit mean and variance
    sorted_small_ens_x_probits -= np.mean( sorted_small_ens_x_probits)
    sorted_small_ens_x_probits /= np.std( sorted_small_ens_x_probits, ddof=1 )

    # Generate virtual probit ensemble in x
    expanded_ens_x = np.zeros([expanded_ens_size], dtype=float)
    expanded_ens_x[:small_ens_size] = deepcopy( sorted_small_ens_x[:] )
    resampling_coeffs = compute_gaussian_resampling_coefficients(small_ens_size, expanded_ens_size-small_ens_size)
    expanded_ens_x[small_ens_size:] = np.array(
        np.matrix( [sorted_small_ens_x_probits] ) * np.matrix( resampling_coeffs )
    )[0,:]

    # Invert PPI transform on the virtual probits to generate virtual members
    if (pese_marginal_name == 'rank histogram with flat tails'):
        expanded_ens_x[small_ens_size:] = invPPI_rank_histogram_with_flat_tails( 
            sorted_small_ens_x, expanded_ens_x[small_ens_size:], bounds=[-1,3]
        )
    elif (pese_marginal_name == 'normal'):
        expanded_ens_x[small_ens_size:] = (
            expanded_ens_x[small_ens_size:] * small_std
            + small_avg
        ) 
    elif (pese_marginal_name == 'gamma'):
        expanded_ens_x[small_ens_size:] = fitted_dist.ppf(
            norm.cdf(expanded_ens_x[small_ens_size:])
        )




    # COMPUTE OBSERVABLE QUANTITIES FOR BOTH SMALL AND EXPANDED ENSEMBLES
    # --------------------------------------------------------------------
    expanded_ens_y = obs_operator( expanded_ens_x )
    sorted_small_ens_y = obs_operator( sorted_small_ens_x )


    # SORT EXPANDED ENSEMBLE
    # -----------------------
    sorted_expanded_ens_x = np.sort( expanded_ens_x )
    sorted_expanded_ens_y = np.sort( expanded_ens_y )


    # GENERATE DICTIONARY TO HOLD EMPIRICAL CDF INFO
    # -----------------------------------------------
    ecdf_dict = {}
    ecdf_dict['small'] = {}
    ecdf_dict['expanded'] = {}


    # COMPUTE EMPIRICAL X CDFS FOR EXPANDED AND SMALL ENS
    # ---------------------------------------------------
    ecdf_dict['small']['xcoord'] , ecdf_dict['small']['xcdf'] = \
        generate_empirical_cdf( sorted_small_ens_x )

    ecdf_dict['expanded']['xcoord'] , ecdf_dict['expanded']['xcdf'] = \
        generate_empirical_cdf( sorted_expanded_ens_x )
    

    # COMPUTE EMPIRICAL Y CDFS FOR EXPANDED AND SMALL ENS
    # ---------------------------------------------------
    ecdf_dict['small']['ycoord'] , ecdf_dict['small']['ycdf'] = \
        generate_empirical_cdf( sorted_small_ens_y )

    ecdf_dict['expanded']['ycoord'] , ecdf_dict['expanded']['ycdf'] = \
        generate_empirical_cdf( sorted_expanded_ens_y )
    


    # COMPUTE WASSERSTEIN DISTANCES FOR EMPIRICAL X CDFS
    # --------------------------------------------------
    ecdf_dict['small']['x cdf dist'] = compute_areal_dist_of_two_curves( 
        ref_x_coord, ref_x_cdf, ecdf_dict['small']['xcoord'] , ecdf_dict['small']['xcdf']
    )
    ecdf_dict['expanded']['x cdf dist'] = compute_areal_dist_of_two_curves( 
        ref_x_coord, ref_x_cdf, ecdf_dict['expanded']['xcoord'] , ecdf_dict['expanded']['xcdf']
    )

    # COMPUTE WASSERSTEIN DISTANCES FOR EMPIRICAL Y CDFS
    # --------------------------------------------------
    ecdf_dict['small']['y cdf dist'] = compute_areal_dist_of_two_curves( 
        ref_y_coord, ref_y_cdf, ecdf_dict['small']['ycoord'] , ecdf_dict['small']['ycdf']
    )
    ecdf_dict['expanded']['y cdf dist'] = compute_areal_dist_of_two_curves( 
        ref_y_coord, ref_y_cdf, ecdf_dict['expanded']['ycoord'] , ecdf_dict['expanded']['ycdf']
    )


    # COMPUTE DISTANCE BETWEEN REFERENCE X-Y CURVE AND ENSEMBLE X-Y CURVE
    # -------------------------------------------------------------------
    ref_xy_coord = np.linspace(-10,10, 10000)
    ecdf_dict['small']['xy dist'] = compute_areal_dist_of_two_curves( 
        ref_xy_coord, obs_operator(ref_xy_coord), sorted_small_ens_x, sorted_small_ens_y
    )
    ecdf_dict['expanded']['xy dist'] = compute_areal_dist_of_two_curves( 
        ref_xy_coord, obs_operator(ref_xy_coord), sorted_expanded_ens_x, sorted_expanded_ens_y
    )


    return ecdf_dict, sorted_expanded_ens_x, sorted_expanded_ens_y


















'''
    OBSERVATION OPERATORS
'''

def hx_signed_power( x, pow=0.5 ):
    xp = x #*0.7+ ref_dist.rvs(size=1)[0]*0.2
    return np.sign(xp) * np.power( np.abs(xp), pow ) 

hx = hx_signed_power













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











'''
    PROBIT PROBABILITY INTEGRAL TRANSFORM & INVERSE-TRANSFORM FUNCTIONS FOR RANK HISTOGRAM WITH FLAT TAILS
'''
# PPI transform for Rank Histogram with Flat Tails
# Arguments: 
# 1) Vector of values defining the PDF (sorted_samples1d)
def PPI_rank_histogram_with_flat_tails( sorted_samples1d ):

    # Assign quantiles
    nEns = sorted_samples1d.shape[0]
    quantiles = (np.arange(nEns)+1.)/(nEns+1.)

    # Map quantiles to probit space
    return norm.ppf(quantiles) 


# Inverse PPI transform for Rank Histogram with Flat Tails
# Arguments: 
# 1) Vector of values defining the PDF (sorted_samples1d)
# 2) Vector of probits to apply inverse PPI on (probits)
# 3) Boundaries indicating the support of the PDF (bounds)
def invPPI_rank_histogram_with_flat_tails( sorted_samples1d, probits, bounds=[-10,10] ):

    # Assign quantiles
    nEns = sorted_samples1d.shape[0]
    bounded_pts = np.insert( sorted_samples1d, [0,nEns], bounds )
    bounded_quantiles = (np.arange(nEns+2))/(nEns+1.)
    bounded_quantiles[-1] =  1.0

    # Move probits to quantile space
    quantiles = norm.cdf( probits )


    # Map quantiles to native space
    return np.interp( quantiles, bounded_quantiles, bounded_pts)









'''
    FUNCTION TO GENERATE EMPIRICAL CDFS

    ArgumentS:
    1) 1D array of sorted sample values
    2) Boundaries of the approximate support of the distribution
'''
def generate_empirical_cdf( sorted_samples1d, bounds = [-10,10] ):
    
    nSamples = sorted_samples1d.shape[0]

    coords = np.insert(
        np.ravel( [sorted_samples1d-1e-5, sorted_samples1d], 'F' ),
        [0,nSamples*2], [-10,10]
    )
    cdf = np.insert(
        np.ravel( [ np.arange(nSamples)*1./nSamples, 
                    (np.arange(nSamples)+1.)/nSamples ], 
                'F'),
        [0,nSamples*2], [0,1.]
    )

    return coords, cdf







'''
    FUNCTION TO COMPUTE AREAL DISTANCE BETWEEN TWO CURVES

    Arguments:
    1) Coordinates of first curve (x1)
    2) First curve (y1)
    3) Coordinates of second  curve (x2)
    4) Second curve (y2)

    Assumes both curves have the same starting and ending coordinate values
    In other words, assume, x1[0] = x2[0] and x1[-1] = x2[-1]
'''
def compute_areal_dist_of_two_curves( x1, y1, x2, y2 ):
    
    # Determine the reference coordinate
    if (len(x1) > len(x2)):
        ref_x = x1; ref_y = y1
        arg_x = x2; arg_y = y2

    else:
        ref_x = x2; ref_y = y2
        arg_x = x1; arg_y = y1

    # Interpolate to reference coordinate and compute abs difference at every ref point
    absdiff = np.abs( np.interp( ref_x, arg_x, arg_y ) - ref_y )

    # Compute areal distance
    dx = ref_x[1:]-ref_x[:-1]
    dist = np.sum( (absdiff[:-1]+absdiff[1:])*dx*0.5 )

    return dist







'''
    FUNCTION: Generate RHF-approximated likelihood function

    Will approximate likelihood function in-btwn members with constant values
    Outside of members, approx likelihood functino has gaussian tails
'''
def rhf_approx_likelihood( y_ref_values, ref_likelihood_values, ens_y_values ):

    # Sort ens y values
    sorted_ens_y_values = np.sort( ens_y_values )

    # Evaluate likelihood function at the ens values
    ens_likelihood = np.interp( sorted_ens_y_values, y_ref_values, ref_likelihood_values )

    # Determine the interior likelihood density
    interior_likelihood_density = (ens_likelihood[1:] + ens_likelihood[:-1])/2.
    left_likelihood_density = ens_likelihood[0]/(sorted_ens_y_values[0]-y_ref_values[0])
    right_likelihood_density = ens_likelihood[-1]/(y_ref_values[-1]-sorted_ens_y_values[-1])

    # Prepare interpolation points
    ens_y_values_left_stag  = sorted_ens_y_values-1e-6
    approx_likelihood_left_stag = sorted_ens_y_values*0.
    approx_likelihood_left_stag[0] = left_likelihood_density
    approx_likelihood_left_stag[1:] = interior_likelihood_density

    ens_y_values_right_stag = sorted_ens_y_values+1e-6
    approx_likelihood_right_stag = sorted_ens_y_values*0.
    approx_likelihood_right_stag[:-1] = interior_likelihood_density
    approx_likelihood_right_stag[-1] = right_likelihood_density

    # Interpolate!
    ens_y_stag_coords = np.ravel([ens_y_values_left_stag, ens_y_values_right_stag],'F')
    approx_likelihood_no_interp = np.ravel([approx_likelihood_left_stag, approx_likelihood_right_stag],'F')
    interp_likelihood = np.interp( y_ref_values, ens_y_stag_coords, approx_likelihood_no_interp)
    

    interp_likelihood /= (np.sum(interp_likelihood)*(y_ref_values[1]-y_ref_values[0]))

    return interp_likelihood










































'''
    MAIN PROGRAM
'''
if __name__ == '__main__':

    # Generate reference samples 
    np.random.seed(0)

    # Generate reference obs distribution via Monte Carlo
    nRefs = 1000000
    ref_x_samples = ref_dist.rvs( size=nRefs )
    ref_y_samples = hx( ref_x_samples )
    ref_y_samples_sorted = np.sort( ref_y_samples )
    ref_x_samples_sorted = np.sort( ref_x_samples )

    expanded_ens_size = 1000

    # Generate reference CDFs
    ref_x_coord = np.insert(
        np.ravel( [ref_x_samples_sorted-1e-5, ref_x_samples_sorted], 'F' ),
        [0,nRefs*2], [-10,10]
    )
    ref_x_cdf = np.insert(
        np.ravel( [ np.arange(nRefs)*1./nRefs, 
                    (np.arange(nRefs)+1.)/nRefs ], 
                'F'),
        [0,nRefs*2], [0,1.]
    )
    ref_y_coord = np.insert(
        np.ravel( [ref_y_samples_sorted-1e-5, ref_y_samples_sorted], 'F' ),
        [0,nRefs*2], [-10,10]
    )
    ref_y_cdf = ref_x_cdf*1.

    # Generate small ensemble
    sorted_small_ens_x = np.sort( ref_dist.rvs( size = small_ens_size ) )
    sorted_small_ens_y = hx( sorted_small_ens_x )
    

    # Generate expanded ensemble and useful stats
    ecdf_dict, sorted_expanded_ens_x, sorted_expanded_ens_y = (
        single_test_pese_gaussian_copula( 
            ref_x_coord, ref_x_cdf, ref_y_coord, ref_y_cdf, 
            sorted_small_ens_x, hx,
            expanded_ens_size=expanded_ens_size, 
            pese_marginal_name=pese_pdf_name
        )
    )





    # Plot reference distribution & obs likelihood function and drawn members overlaid
    xvals = np.linspace(-3,5,1000)
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(4,5) )
    axs[0].plot( xvals, ref_dist.pdf(xvals), ':k', linewidth=2, label='True Curve')

    if (small_ens_size <= 20):
        axs[0].scatter( sorted_small_ens_x, sorted_small_ens_x*0., c='r', s=100, marker='x', zorder=100, label = 'Initial Ensemble Members')

    axs[0].set_ylabel('PDF')
    axs[0].set_xlabel('x')
    axs[0].set_title( 'True Fcst PDF (skewed normal)')
    axs[0].set_xlim([-2,4])



    # Plot observation likelihood function, and overlay the RHF approximations
    yvals = np.linspace(-100,100,10000)
    likelihood_dist = norm(loc=1,scale=0.5)


    # Plot fcst member approximation
    axs[1].plot( yvals, 
            rhf_approx_likelihood( yvals, likelihood_dist.pdf(yvals), sorted_small_ens_y*1. ),
            '-r', linewidth=2, label='Fcst. Ens. Curve'
    )

    # Plot virtual member approximation
    axs[1].plot( yvals, 
            rhf_approx_likelihood( yvals, likelihood_dist.pdf(yvals), sorted_expanded_ens_x ),
            color='dodgerblue', linewidth=2, label='Virt. Mems. Curve'
    )

    # Plot actual likelihood function
    axs[1].plot( yvals, likelihood_dist.pdf(yvals), ':k', linewidth=2, label='True Curve')


    if (small_ens_size <= 20):
        axs[1].scatter( sorted_small_ens_y, sorted_small_ens_y*0., c='r', s=100, marker='x', zorder=100, label = 'Fcst. Ens. Mems.')

    fig.subplots_adjust( left=0.2, bottom=0.35, right=0.95, top=0.9, hspace=0.75)
    axs[1].set_ylabel('PDF')
    axs[1].set_xlabel('y')
    axs[1].set_title( 'Obs. Likelihood Function [p(y|x)]')
    axs[1].set_xlim([-2,4])

    plt.legend( ncol=2, fontsize="10", bbox_to_anchor = (0.32,-0.7,0.7,0.1))

    plt.savefig('DEMO_PESE_impacts_on_prior_stats/DEMO_ref_xPDF_and_obs_likelihood.svg')
    plt.close()














   
    # Generate figure to hold single illustative test
    fig, axs0 = plt.subplots(nrows=2, ncols=2, figsize=(6,4))
    fig.delaxes(axs0[0,1])
    axs = [axs0[0,0], axs0[1,1], axs0[1,0]]

    # Plot Reference EMPIRICAL CDF in model space
    axs[0].plot( ref_x_coord, ref_x_cdf, 
                 color='k',  zorder=10, linewidth=3, linestyle=':' )

    # Plot small ens EMPIRICAL CDF in model space   
    axs[0].plot( ecdf_dict['small']['xcoord'] , ecdf_dict['small']['xcdf'], 
                 color='r',  zorder=2, linewidth=3, linestyle='-' )
    
    # Plot expanded ens EMPIRICAL CDF in model space
    axs[0].plot(  ecdf_dict['expanded']['xcoord'] , ecdf_dict['expanded']['xcdf'], 
                 color='dodgerblue',  zorder=1, linewidth=3, linestyle='-' )

    axs[0].set_xlim([-3,5])

    if (small_ens_size <= 20):
        axs[0].scatter( sorted_small_ens_x, sorted_small_ens_x*0-0.05, c='r', s=100,  zorder=100,marker='x')
    axs[0].set_title('CDF of x')
    axs[0].set_xlabel('x')
    axs[0].set_ylabel('CDF')



    # Plot Reference EMPIRICAL CDF in obs space
    axs[1].plot( ref_y_coord, ref_y_cdf, 
                 color='k',  zorder=10, linewidth=3, linestyle=':' )

    # Plot small ens EMPIRICAL CDF in obs space   
    axs[1].plot( ecdf_dict['small']['xcoord'] , ecdf_dict['small']['xcdf'], 
                 color='r',  zorder=2, linewidth=3, linestyle='-' )
    
    # Plot expanded ens EMPIRICAL CDF in obs space
    axs[1].plot(  ecdf_dict['expanded']['ycoord'] , ecdf_dict['expanded']['ycdf'], 
                 color='dodgerblue',  zorder=1, linewidth=3, linestyle='-' )


    axs[1].set_xlim([-3,4])

    if (small_ens_size <= 20):
        axs[1].scatter( sorted_small_ens_x, sorted_small_ens_y*0-0.05,  c='r', s=100,  zorder=100,marker='x')
    axs[1].set_title('CDF of y [=h(x)]')
    axs[1].set_xlabel('y [=h(x)]')





    # Generate probit regression curves in native space assuming bi-jective function
    ref_x_curve = np.linspace(-4,5,1000)
    ref_y_curve = hx( ref_x_curve )
    axs[2].plot( ref_x_curve, ref_y_curve, color='k',  zorder=10,linestyle=':', linewidth=3,  label='True Curve' )

    small_y_curve = np.interp( 
        ref_x_curve, sorted_small_ens_x, sorted_small_ens_y
    )
    axs[2].plot( ref_x_curve, small_y_curve, color='r', zorder=10, linewidth=3, label='Fcst. Ens. Curve' )

    if (small_ens_size <= 20):
        axs[2].scatter( sorted_small_ens_x, sorted_small_ens_y, 
                        c='r', s=100, marker='x', zorder=100, label='Fcst. Ens. Mems.')


    expanded_y_curve = np.interp( 
        ref_x_curve, sorted_expanded_ens_x, sorted_expanded_ens_y
    )
    axs[2].plot( ref_x_curve, expanded_y_curve, color='dodgerblue', zorder=5, linewidth=3, label='Virt. Mems. Curve'  )

    axs[2].set_xlim([-3,4])
    axs[2].set_ylim([-2,2.5])

    axs[2].set_title('y [=h(x)]   VS     x')
    axs[2].set_xlabel('x')
    axs[2].set_ylabel('y [=h(x)]')





    fig.subplots_adjust(bottom=0.15, top=0.9, right=0.95, hspace=0.6, wspace=0.3)


    handles, labels = axs[2].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower left', bbox_to_anchor=(0.57,0.6), ncol=1, fontsize="12")



    # Save figure
    plt.savefig('DEMO_PESE_impacts_on_prior_stats/DEMO_user_marginals_ILLUSTRATION_%s_EnsSize%04d.svg' % (pese_pdf_name.replace(' ','_'), small_ens_size) )
    plt.close()





 


    '''
        REPEAT TESTS AN ABSURD NUMBER OF TIMES
    '''
    nTrials = 1000

    # Dictionary to hold trials
    performance_metric_dict = {}
    performance_metric_dict['small x cdf'] = np.zeros(nTrials)
    performance_metric_dict['small y cdf'] = np.zeros(nTrials)
    performance_metric_dict['small xy curve'] = np.zeros(nTrials)
    performance_metric_dict['expanded x cdf'] = np.zeros(nTrials)
    performance_metric_dict['expanded y cdf'] = np.zeros(nTrials)
    performance_metric_dict['expanded xy curve'] = np.zeros(nTrials)


    # Loop over trials
    print_interval = int(nTrials/10)
    for tt in range(nTrials):

        # Run test
        sorted_small_ens_x = np.sort( ref_dist.rvs( size = small_ens_size ) )
        ecdf_dict, sorted_expanded_ens_x, sorted_expanded_ens_y = (
            single_test_pese_gaussian_copula( 
                ref_x_coord, ref_x_cdf, ref_y_coord, ref_y_cdf, 
                sorted_small_ens_x, hx,
                expanded_ens_size=expanded_ens_size, 
                pese_marginal_name=pese_pdf_name
            )
        )
        if ( (tt+1) % print_interval == 0 ):
            print('%5d out of %5d trials done' % (tt+1, nTrials))

        # Store Wasserstein metrics
        performance_metric_dict['small x cdf'][tt] = ecdf_dict['small']['x cdf dist']
        performance_metric_dict['small y cdf'][tt] = ecdf_dict['small']['y cdf dist']
        performance_metric_dict['small xy curve'][tt] = ecdf_dict['small']['xy dist']
        performance_metric_dict['expanded x cdf'][tt] = ecdf_dict['expanded']['x cdf dist']
        performance_metric_dict['expanded y cdf'][tt] = ecdf_dict['expanded']['y cdf dist']
        performance_metric_dict['expanded xy curve'][tt] = ecdf_dict['expanded']['xy dist']

    # ---- End of loop over trials

    # Print out useful info
    for stype in ['x cdf','y cdf','xy curve']:
        for etype in ['small','expanded']:
            
            print('Performance of %10s for %10s ensemble: %10f (std err: %10f)' 
                   % ( stype, etype, 
                       np.mean( performance_metric_dict['%s %s' % (etype, stype)]),
                       np.std( performance_metric_dict['%s %s' % (etype, stype)], ddof=1)/np.sqrt(nTrials)
                      )
            )