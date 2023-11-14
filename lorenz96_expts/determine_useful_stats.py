#!/bin/python

# =============================================
# SCRIPT TO CONSTRUCT USEFUL SUMMARY STATISTICS
# =============================================

# Load useful modules
import numpy as np
from netCDF4 import Dataset as ncopen
from copy import deepcopy

# Lorenz model forcing
F = 8.0




# Function to compute RMS dxdt of an ensemble
# -------------------------------------------
def rms_dxdt( ncfile ):

    # Load ensemble of space-time states
    x_ens = deepcopy( ncfile.variables['state'][500:,:,:] )

    # Initialize memory to hold rate of change 
    rate = np.zeros( x_ens.shape )

    # Compute advection term in interior points
    rate[:, :, 2:-1] = (x_ens[:, :, 3:] - x_ens[:, :, :-3]) * x_ens[:, :, 1:-2]
    # Compute advection term in boundary points
    rate[:, :,    0] = (x_ens[:, :,  1] - x_ens[:, :,  -2]) * x_ens[:, :,   -1]
    rate[:, :,    1] = (x_ens[:, :,  2] - x_ens[:, :,  -1]) * x_ens[:, :,    0]
    rate[:, :,   -1] = (x_ens[:, :,  0] - x_ens[:, :,  -3]) * x_ens[:, :,   -2]
    # Add in self-decay term
    rate -= x_ens
    # Add in forcing term
    rate += F

    # Compute time-ensemble averaged RMS rate
    rms_rate = np.mean(
        np.sqrt(
            np.mean(
                rate**2, axis=-1
            )
        )
    )
    
    return rms_rate 
# -------------- End of function to compute time-ensemble averaged RMS rate






# Function to compute RMSE and RMSS
# ----------------------------------
def compute_performance_stats( f0, f_truth ): 

    # Compute t-series of RMSE
    rmse_tseries = (
        np.sqrt(
            np.mean(
                (
                    np.mean( f0.variables['state'][500:,:,:], axis=1 )
                    - f_truth.variables['state'][500:,0,:]
                )**2
                , axis=-1
            )
        )
    )

    # Compute t-series of RMSS
    rmss_tseries = (
        np.sqrt(
            np.mean(
                np.var( f0.variables['state'][500:,:,:], 
                        axis=1, ddof=1
                ), axis=-1
            )
        )
    )

    # Compute time-averaged CR, RMSE & bias
    rmse_avg = np.mean( rmse_tseries )
    cr_avg = np.mean( rmss_tseries/rmse_tseries )
    bias_avg = (
        np.mean( 
            np.mean( f0.variables['state'][500:,:,:], axis=1 )
            - f_truth.variables['state'][500:,0,:]
        )
    )

    return rmse_avg, bias_avg, cr_avg

# --------------- End of function to compute performance stats






# Compute RMSE & RMSS of prior and posterior
# -------------------------------------------
f_prior = ncopen('forecast.nc', 'r')
f_truth = ncopen('true_state.nc', 'r')
f_poste = ncopen('analysis.nc','r')

rmse_prior, bias_prior, cr_prior = compute_performance_stats( f_prior, f_truth )
rmse_poste, bias_poste, cr_poste = compute_performance_stats( f_poste, f_truth )


# Compute RMS of rate for posterior ensembles
rms_dxdt_poste = rms_dxdt( f_poste )




# Save useful statistics
# -----------------------
f = open('useful_stats.txt','w')
string = ""
string += (
    "xf_RMSE: %f    xf_CR(=RMSS/RMSE): %f    xf_BIAS: %f\n" 
    % ( rmse_prior, cr_prior, bias_prior )
)
string += (
    "xa_RMSE: %f    xa_CR(=RMSS/RMSE): %f    xa_BIAS: %f\n" 
    % ( rmse_poste, cr_poste, bias_poste )
)
string += "RMS_dxdt: %f" % rms_dxdt_poste
f.write(string)
f.close()