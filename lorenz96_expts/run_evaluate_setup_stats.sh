#!/bin/bash

# =============================================================================
# SCRIPT TO COMPUTE THE PERFORMACE OF EVERY TRIAL INDICATDE BY CONFIG_FILE
# ------------------------------------------------------------------------------


# Load configuration files
CONFIG_FILE=$1
. $CONFIG_FILE



# Translate list of ROI values from space-separated to comma-separated
csv_roi=""
for roi in $loc_halfradius_list; do
    csv_roi=$roi,$csv_roi
done



# Iterate over all cycling intervals
for cyc_interval_seconds in $cyc_interval_seconds_list; do 

    # Iterate over ensemble sizes
    for ens_size in $ens_size_list; do

        # Iterate over obs error variances
        for obs_err_var in $obs_err_var_list; do

            # Iterate over filter types
            for filter_type in $filter_type_list; do

                # Iterate over observation counts
                for n_obs in $nObs_list; do

                    # Iterate over PESE usage
                    for use_pese in $use_pese_list; do

                        results_string=""

                        # Iterate over trials
                        for itrial in `seq -f %03g 1 $nTrials`; do

                            # Generate path of trial directory
                            trial_dir=$dir_expts/$obs_name"/filter"$filter_type"/ObsErrVar"$obs_err_var"/CycInt"$cyc_interval_seconds"/EnsSize"$ens_size"/ObsNum"$n_obs"/PESE"$use_pese/Trial$itrial
                            
                            # Enter trial directory
                            cd $trial_dir

                            results_string=`python $dir_bash/evaluate_setup_stats.py $csv_roi`" "$results_string


                        done # --- End of loop over trials

                        # Generate label & output results
                        printf "%-20s %-10s %-15s %-20s %-20s %-10s\n" "filter_type:$filter_type" "n_obs:$n_obs" "ens_size:$ens_size" "obs_err_var:$obs_err_var" "cyc_int:$cyc_interval_seconds" "PESE:$use_pese"
                        echo $results_string

                    done # -- End of loop over PESE settings
                done # ---- End of loop over number of obs
            done # ---- End of loop over all filters
        done # ---- End of loop over all obs error variances
    done # ---- End of loop over all ensemble sizes
done # ---- End of loop over all cycling intervals