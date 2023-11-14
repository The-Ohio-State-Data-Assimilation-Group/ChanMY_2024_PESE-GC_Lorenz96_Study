'''
    SCRIPT TO GENERATE SPUN-UP STATES FOR LORENZ 1996 MODEL
'''

import numpy as np
from netCDF4 import Dataset as ncopen
import pylib_Lorenz96_with_tracer as pyL96
import time
import shutil


# FIXED CONSTANTS
num_spinup_steps = 100  # Number of time steps for spinup
ens_size = 80           # Number of ensemble members for each set (number of members in TEMPLATE_filter_input.nc)
num_grid_pts = 40       # Number of Lorenz 96 grid points

# Number of init conditions to create
id_list = np.arange(72)+1



# Timer variable
t0 = time.time()


# Loop over all ids
for set_id in id_list:

    # Generate and spin-up massive ensemble
    np.random.seed( set_id )
    state_ensemble = np.random.normal( size=(num_grid_pts, ens_size+1) )
    pyL96.multistep( state_ensemble, num_spinup_steps )

    print('Spin-up completed for %03d / %03d init conditions (%5d sec elapsed)'
            % (set_id+1, len(id_list), time.time() - t0 ) )
    



    # Generate truth file 
    # -------------------

    # File name of truth file
    truth_fname = 'init_conditions/id%03d_perfect_input.nc' % set_id

    # Copy truth file
    shutil.copy( 'TEMPLATE_perfect_input.nc', truth_fname )

    # Overwrite state in truth file
    f_true = ncopen( truth_fname, 'a' )
    f_true.variables['state'][0] = state_ensemble[:,0]

    # Flush updated truth file to disk
    f_true.close()

    print('Wrote perfect_input for %03d / %03d init conditions (%5d sec elapsed)'
            % (set_id, len(id_list), time.time() - t0 ) )


    # Generate filter_input file 
    # --------------------------

    # File name of truth file
    filter_fname = 'init_conditions/id%03d_filter_input.nc' % set_id

    # Copy truth file
    shutil.copy( 'TEMPLATE_filter_input.nc', filter_fname )

    # Overwrite state in filter_input file
    f_filt = ncopen( filter_fname, 'a' )
    f_filt.variables['state'][0,:,:] = (state_ensemble[:,1:]).T

    # Flush updated truth file to disk
    f_filt.close()
    print('Wrote filter_input for %03d / %03d init conditions (%5d sec elapsed)'
            % (set_id, len(id_list), time.time() - t0 ) )
