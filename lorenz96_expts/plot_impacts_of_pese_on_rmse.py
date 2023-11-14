'''
    SCRIPT TO COMPARE RMSE PERFORMANCES OF VARIOUS SET-UPS BY READING IN A TEXT FILE CONTAINING ASYMPTOTIC RMSE VALUES
    
    This script plots out the <RMSE> quantity described in Chan's manuscript on "Improving EnsDA with PESE-GC"

    RMSE array is 7-dimensional. Dims: [ filter type, num of obs, ens size, obs err, cycling interval, PESE usage, trials ]
'''

import numpy as np
from matplotlib import use
use('agg')
from matplotlib.colors import LogNorm as mcolor_LogNorm
import matplotlib.pyplot as plt
from sys import argv
from copy import deepcopy
import time 

# User-inputted path to text file containing asymptotic RMSEs
textfile_path=argv[1]

# Create figure title
# setup_title = (textfile_path.split('/'))[1]
# setup_title = (setup_title.split('.'))[0]
setup_title = ""




# Pretty names for various setups
pretty_names = {
    '1_xGauss_yGauss': 'EAKF',
    '2_xGauss_yGauss': 'EnKF',
    '8_xGauss_yGauss': 'RHF',
    '8_xBNRH_yBNRH': 'PR',
}












# Load and interpret textfile
# ----------------------------
t0 = time.time()
print('\n%-100s (%10.4f seconds elapsed)' % ("Starting to load text file data", time.time()-t0))



# Read all lines
f = open( textfile_path)
all_lines = f.readlines()
f.close()
num_lines = len(all_lines)


# Determine all the setups contained in the text file
all_info = {}
all_info['filter_type'] = [] ; all_info['n_obs'] = []
all_info['ens_size'] = []     ; all_info['obs_err_var'] = []
all_info['cyc_int'] = []      ; all_info['PESE'] = []

for line in all_lines[::2]:
    entry = []                  # Entry for setup_info_lines
    multiitems = line.split()        # Split line by whitespaces
    for item in multiitems:
        elements = item.split(":")
        if ( not elements[1] in all_info[elements[0]] ):
            all_info[elements[0]].append( elements[1] )
        # --- End of if statement
    # --- End of loop over multiitems
# --- End of loop over lines 



# Sort the filter names in order of filter type and PESE factor
all_info['filter_type'].sort()
buffer_list = []
for name in all_info['filter_type']:
    if ( name[-3:] == '20X' ):
        buffer_list.append(name)
for name in all_info['filter_type']:
    if ( name[-3:] == '10X' ):
        buffer_list.append(name)
for name in all_info['filter_type']:
    if ( name[-3:] == '05X' ):
        buffer_list.append(name)
all_info['filter_type'] = buffer_list


# Swap the 8_xBNRH and 8_xGauss order
for i in range(3):
    buffer_list[ i*4+3 -1 ], buffer_list[ (i+1)*4 -1] = buffer_list[ (i+1)*4 -1], buffer_list[ i*4+3 -1 ]
# ----- End of swappin loop





# Determine number of trials
num_trials = len( all_lines[1].split())



# Set up array to hold all statistics
# various stats arrays are 7-dimensional. Dims: [ filter type, num of obs, ens size, obs err, cycling interval, PESE usage, trials ]
all_info['trial rmses'] = np.zeros( 
    [ len(all_info['filter_type']), len(all_info['n_obs']), len(all_info['ens_size']), 
      len(all_info['obs_err_var']), len( all_info['cyc_int']), len( all_info['PESE']),
      num_trials ]
)

all_info['trial cr'] = np.zeros( 
    [ len(all_info['filter_type']), len(all_info['n_obs']), len(all_info['ens_size']), 
      len(all_info['obs_err_var']), len( all_info['cyc_int']), len( all_info['PESE']),
      num_trials ]
)

all_info['trial bias'] = np.zeros( 
    [ len(all_info['filter_type']), len(all_info['n_obs']), len(all_info['ens_size']), 
      len(all_info['obs_err_var']), len( all_info['cyc_int']), len( all_info['PESE']),
      num_trials ]
)

all_info['trial dxdt'] = np.zeros( 
    [ len(all_info['filter_type']), len(all_info['n_obs']), len(all_info['ens_size']), 
      len(all_info['obs_err_var']), len( all_info['cyc_int']), len( all_info['PESE']),
      num_trials ]
)

all_info['trial roi'] = np.zeros( 
    [ len(all_info['filter_type']), len(all_info['n_obs']), len(all_info['ens_size']), 
      len(all_info['obs_err_var']), len( all_info['cyc_int']), len( all_info['PESE']),
      num_trials ]
)





# Read in statistics of trial 
for iline in range( int(num_lines/2) ):

    # Determine coordinate of the RMSE array
    infoline = all_lines[2*iline].split()
    ifilter = all_info['filter_type'].index( (infoline[0].split(':'))[1] )
    iobs = all_info['n_obs'].index( (infoline[1].split(':'))[1] )
    isize = all_info['ens_size'].index( infoline[2].split(':')[1] )
    ierr = all_info['obs_err_var'].index( infoline[3].split(':')[1] )
    icyc = all_info['cyc_int'].index( infoline[4].split(':')[1] )
    ipese = all_info['PESE'].index( infoline[5].split(':')[1] )

    # Read RMSE
    stats_entries = all_lines[2*iline+1].split(" ")
    if (len( stats_entries ) == num_trials):

        for itrial in range( num_trials ):

            # String of trial stats
            trial_string=stats_entries[itrial]

            # Parse string (see evaluate_setup_stats.py)
            trial_string = trial_string.split(",")
            for ss in range(5):
                sname = "trial " + ["roi","rmses","cr","bias","dxdt"][ss]
                all_info[sname][ifilter, iobs, isize, ierr, icyc, ipese, itrial] = (
                    float( trial_string[ss].split(":")[1])
                )                
    else:
        print('RMSE data for %s is missing' % all_lines[2*iline])
        print('RMSE vals: %s' % all_lines[2*iline+1] )
        for ss in range(5):
            sname = "trial " + ["roi","rmses","cr","bias","dxdt"][ss]
            all_info[sname][ifilter, iobs, isize, ierr, icyc, ipese, itrial] = np.nan
        
    # --- End of loop over trials
# --- End of loop over lines in the text file

print('\n%-100s (%10.4f seconds elapsed)' % ("Finished loading text file data", time.time()-t0))





# Useful mesh for plotting later
ens_mesh, cyc_mesh = np.meshgrid( np.arange(len(all_info['ens_size'])), np.arange(len(all_info['cyc_int'])) )





# Generating absolute biases
all_info['trial abs bias'] = np.abs( all_info['trial bias'] )
































# Now examine impacts of PESE on RMSE
# -----------------------------------

# Determine indices for when PESE is or is not in use.
iPESE_true = all_info['PESE'].index('true')
iPESE_false = all_info['PESE'].index('false')

# Now plot out impacts of PESE
ncols=4; nrows = int(np.ceil(len( all_info['filter_type'])/4))
fig, axs = plt.subplots( nrows=nrows, ncols=ncols, figsize = (1.5*ncols,2*nrows) )
if len( axs.shape ) == 2:
    axs = [ val for sub in axs for val in sub ]


# Compute diffs in RMSE
sname = "trial rmses"
stats_pese_T = all_info[sname][:,:,:,:,:,iPESE_true,:]
stats_pese_F = all_info[sname][:,:,:,:,:,iPESE_false,:]
stats_diff = stats_pese_T - stats_pese_F

# Compute statistics of statistics
stats_diff_avg = np.nanmean( stats_diff, axis=-1 )
stats_diff_std = np.nanstd(  stats_diff, axis=-1, ddof=1)
flags_valid = np.invert( np.isnan( stats_diff ) )
num_valid_trials = np.sum( flags_valid, axis=-1)
stats_diff_ste = stats_diff_std / np.sqrt( num_valid_trials )

# Compute relative difference in statistics 
stats_diff_avg_relative = stats_diff_avg / np.nanmean( stats_pese_F, axis=-1 )

# Determine which differences are significant and appreciable
flags_display = np.abs( stats_diff_avg / stats_diff_ste ) > 2.58

# Color ranges
max_color=20
min_color=0.5

# Go through each filter and output results
for ifilter in range( len( all_info['filter_type']) ):

    # Select panel
    ax = axs[ifilter]

    # Subset relative stats diff and display flags 
    sub_rel_stats_diff = stats_diff_avg_relative[ ifilter, iobs, :, ierr, : ]
    sub_flags_display = flags_display[ifilter, iobs, :, ierr, :]

    # Plotting out places where no valid trials were fonud
    flag_zero_valid = num_valid_trials[ifilter, iobs, :, ierr, : ] < 1
    ax.scatter( ens_mesh.T[flag_zero_valid], cyc_mesh.T[flag_zero_valid], marker='x', c='k', s=200 )

    # Plotting out places where less than 30 trials were fonud
    flag_30 = (num_valid_trials[ifilter, iobs, :, ierr, : ] < 30) * (num_valid_trials[ifilter, iobs, :, ierr, : ] > 0)
    ax.scatter( ens_mesh.T[flag_30], cyc_mesh.T[flag_30], marker='o', facecolors='none', edgecolor='k',s=300 )
        

    # Plotting out positive values (for RMSE: PESE made things worse)
    flags = (sub_rel_stats_diff > 0) * (sub_flags_display)
    sc_bad = ax.scatter( ens_mesh.T[flags], cyc_mesh.T[flags], c=sub_rel_stats_diff[flags]*100, marker='^', cmap='Reds',
                        edgecolor='k', linewidth=1., s=100, norm = mcolor_LogNorm(vmin=min_color, vmax=max_color) )

    # Plotting out negative values (for rmse PESE made things better)
    flags = (sub_rel_stats_diff < 0) * (sub_flags_display)
    sc_good = ax.scatter( ens_mesh.T[flags], cyc_mesh.T[flags], c=sub_rel_stats_diff[flags]*100*-1, marker='v', cmap='Blues',
                        edgecolor='k', linewidth=1., s=100, norm = mcolor_LogNorm(vmin =min_color, vmax=max_color) )
    
    # Label subplot
    for name in pretty_names.keys():
        if ( name in  all_info['filter_type'][ifilter]  ):
            title_string = all_info['filter_type'][ifilter].replace( 
                name, pretty_names[name]
            )
    # End of special handling
    ax.set_title( "%s" %
                    (title_string.replace("pese","") ) )


    # Label y-axis
    ax.set_ylim([-0.5, len(all_info['cyc_int'])-0.5])
    if ( ifilter % 4 == 0 ):
        ax.set_ylabel('Assim. period (sec)' )
        ax.set_yticks(np.arange(len(all_info['cyc_int'])))
        ax.set_yticklabels(all_info['cyc_int'])
    else:
        ax.set_yticks(np.arange(len(all_info['cyc_int'])))
        ax.set_yticklabels(["" for i in np.arange(len(all_info['cyc_int']))])
        
        
    # Label x-axis
    ax.set_xlim([-0.5, len(all_info['ens_size'])-0.5])
    if ( ifilter >= (nrows-1)*(ncols) ):
        ax.set_xlabel('Ens. size')
        ax.set_xticks(np.arange(len(all_info['ens_size'])))
        ax.set_xticklabels(all_info['ens_size'])
        # ax.set_title( 'Num Obs: %s\nObs Err %s' % (all_info['n_obs'][iobs], all_info['obs_err_var'][ierr]) )
    else:
        ax.set_xticks(np.arange(len(all_info['ens_size'])))
        ax.set_xticklabels(["" for i in np.arange(len(all_info['ens_size']))])

# ------- End of loop over filters


# Add super title
fig.subplots_adjust(left=0.15,bottom=0.2,top=0.9, right=0.99, wspace=0.15, hspace=0.4)
plt.suptitle( 'PESE impacts on RMSE' )

# Add colorbars
cbar_ax_good = fig.add_axes([0.175,0.08,0.35,0.025])
cbar_good = plt.colorbar(sc_good, cax=cbar_ax_good, orientation='horizontal')
cbar_ax_good.set_xlabel('$\%$ improvement in RMSE')

cbar_ax_bad = fig.add_axes([0.6,0.08,0.35,0.025])
cbar_bad = plt.colorbar(sc_bad, cax=cbar_ax_bad, orientation='horizontal')
cbar_ax_bad.set_xlabel('$\%$ degradation in RMSE')


# Save figure
plt.savefig('figures/impact_of_pese_on_rmse.svg' )


