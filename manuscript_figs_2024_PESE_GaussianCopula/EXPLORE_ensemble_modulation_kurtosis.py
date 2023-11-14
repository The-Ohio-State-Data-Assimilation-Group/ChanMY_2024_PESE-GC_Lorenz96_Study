'''
    EXPLORATION OF ENSEMBLE MODULATION

    Suppose 1000 element state (i.e., 1000 locations).
    
    Will draw a 100-member Gaussian ensemble from 1000-dimensional Gaussian distribution
    with zero mean and identity covariance.
    
    Localization matrix created by applying Gaussian convolution on identity matrix
'''

import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from scipy.linalg import sqrtm as matrix_sqrt
from scipy.stats import kurtosis, norm
from scipy.ndimage import gaussian_filter1d



'''
    USER SETTINGS
'''
# Number of state elements (i.e., locations)
nX = 1000

# Number of memebrs
nEns = 100

# Initial "white noise draw"
white_noise_ens = norm.ppf( (np.arange(nEns)+0.18)/(nEns) )
print( kurtosis(white_noise_ens, bias=False, fisher=False) )



# Range of localization radii to test
loc_radius_list = []
loc_radius_list.append(0.0001)
loc_len = 1.
while loc_len < 2300:
    loc_radius_list.append(loc_len*1.)
    loc_len *= 1.1

loc_radius_list.append(9999)
loc_radius_list = np.array(loc_radius_list)


loc_matrix = np.eye(nX)*1.


'''
   TEST VARIOUS LOCALIZATION MATRICES AND STORE THE RESULTING KURTOSIS
'''
kurtosis_list = loc_radius_list * 0.

for ll in range(len(loc_radius_list)):

    loc_len = loc_radius_list[ll]


    '''
        Generate localization matrix
    '''

    # Construct gaussian convolution kernel matrix
    loc_matrix[:,:] = np.matrix( gaussian_filter1d(
        np.eye( nX ),
        loc_len/2,                 # Note that sigma here is half of the prior ens's
        mode='wrap',
        axis=0
    ) )[:,:]

    # Construct localization matrix
    loc_matrix[:,:] = loc_matrix * loc_matrix.T

    # Normalize localization matrix (diagonal elements must be one!)
    loc_diag = np.diagflat( np.diag( loc_matrix ) )
    loc_diag_inv = np.matrix( np.linalg.inv(np.sqrt(loc_diag)) )
    loc_matrix[:,:] = loc_diag_inv * loc_matrix * loc_diag_inv


    '''
        Generate square-root of localization matrix
    '''
    if loc_len < 2300:
        loc_square_root = np.matrix( matrix_sqrt( loc_matrix) )
    else:
        loc_square_root = np.ones( [nX,nX])/np.sqrt(nX)

    

    '''
        Generate the modulated ensemble at location site 0
    '''
    prior_ens_site0 = white_noise_ens*1. #/ np.sqrt(nEns-1)
    modulated_ens = np.zeros( nX*nEns )
    for ivec in range(nX):
        for imem in range(nEns):
            modulated_ens[ ivec*nEns + imem -1 ] = loc_square_root[0, ivec ] * prior_ens_site0[imem] * np.sqrt( nX*nEns) / np.sqrt(nEns-1)


    '''
        Store modulated ensemble's kurtosis and move on
    '''
    kurtosis_list[ll] = kurtosis( modulated_ens, bias=False, fisher=False )

    print( 'Processed localization radius %5f , kurtosis is %5d' % (loc_len, np.round(kurtosis_list[ll])))


# ----- End of loop over localization radii





'''
    PLOT KURTOSIS OF MODULATED ENSEMBLE
'''
fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(4,3))
original_kurtosis = kurtosis(prior_ens_site0, bias=False, fisher=False)
ax.plot( loc_radius_list[:-1]/nX, kurtosis_list[:-1], '-r')
ax.axhline( kurtosis_list[-1], color='r', linestyle='--')

#plt.xscale('log')
ax.set_yscale('log')
ax.set_ylim([1,4000])
ax.set_ylabel('Kurtosis')
ax.set_xlabel('( Localization Length Scale )/$N_x$')
ax.set_title('Modulated Ensemble\'s Kurtosis')
fig.subplots_adjust(left=0.2,bottom=0.2, top=0.8)

plt.savefig('EXPLORE_ensemble_modulation/kurtosis_of_modulated_ens.svg')
