'''
    SCRIPT TO EXPLORE THE OVERLAP BETWEEN GAUSSIAN FORECAST PDF AND THE STRONGLY-CURVED REGION OF THE OBS LIKELIHOOD FUNCTION
'''

import numpy as np
from scipy.stats import norm
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt


'''
    USER CONTROLS
'''
cr_list = np.linspace(0.5,1.5,51)
obs_err_var = 1.
fcst_err_var_list =  np.linspace(0.01,5.01,201)


'''
    COMPUTE PROBABILITY OF FORECAST MEMBERS SAMPLING THE STRONG-CURVED REGION OF THE LIKELIHOOD FUNCTION
'''
# Determine observation values (depends on CR value, fcst err var and obs err var)
fcst_err_mesh, cr_mesh = np.meshgrid( fcst_err_var_list, cr_list )
fcst_avg = 0.
obs_val_mesh = (np.sqrt( (fcst_err_mesh + obs_err_var) ) / cr_mesh )
obs_region_left_bound = obs_val_mesh - np.sqrt(obs_err_var)*1.5
obs_region_right_bound = obs_val_mesh + np.sqrt(obs_err_var)*1.5


# Compute sampling probabilities
sampling_probabilities = fcst_err_mesh * 0.
for ierr in range( len(fcst_err_var_list) ):
    for icr in range( len(cr_list) ):
        sampling_probabilities[icr, ierr] = (
            norm.cdf(
                obs_region_right_bound[icr, ierr],
                loc=0, scale = np.sqrt( fcst_err_var_list[ierr] )
            )
            -
            norm.cdf(
                obs_region_left_bound[icr, ierr],
                loc=0, scale = np.sqrt( fcst_err_var_list[ierr] )
            )
        )


'''
    PLOT SAMPLING PROBABILITIES
'''

fig = plt.figure(figsize=(4,3))
sampling_probabilities = np.ma.masked_where( sampling_probabilities < 0.01, sampling_probabilities )
cnf = plt.contourf( fcst_err_mesh, cr_mesh, sampling_probabilities, np.linspace(0.01,0.91,11), cmap = 'inferno_r',
                     extend='max' )
cbar = plt.colorbar(cnf)
plt.xlabel('Fcst Variance')
plt.ylabel('CR (spread / innovation)')
fig.subplots_adjust(top=0.9, bottom=0.3, left=0.2)
plt.savefig('EXPLORE_Gaussian_overlaps/EXPLORE_Gaussian_overlaps.svg')
