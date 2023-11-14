#!/bin/bash

# ======================================================
# RUN EXPERIMENTS WITH VARIOUS SETTINGS
# ======================================================

# Configuration file to work with
CONFIG_FILE=NAME_OF_CONFIG_FILE

# Generate job name
jobname=`echo $CONFIG_FILE | cut -c 29-`

. $CONFIG_FILE

# Number of active processes
nprocs=0



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

                    # Iterate over PESE settings
                    for use_pese in $use_pese_list; do

                        # Iterate over all trials
                        for trial in `seq 1 $nTrials`; do

                            # Run experiments with tuning
                            ./run_expt_with_roi_tuning.sh $CONFIG_FILE $filter_type $n_obs $ens_size \
                                $obs_err_var $cyc_interval_seconds $use_pese $trial                  \
                                >> $dir_bash"/proc_logs/log."$jobname"_proc_"`printf "%03d" $(($nprocs+1))` &
                            
                            # Update number of active processes
                            nprocs=$(($nprocs+1))

                            # Output info for user
                            msgstring=`date`"  --  Proc "`printf %03d $nprocs`": processing filter_type:$filter_type"
                            msgstring=$msgstring" n_obs:$n_obs ens_size:$ens_size obs_err_var:$obs_err_var"
                            msgstring=$msgstring" cyc_int:$cyc_interval_seconds PESE:$use_pese trial_no:$trial"
                            echo $msgstring

                            # Wait if there are no idle processes
                            if [[ $nprocs -ge $tot_nprocs ]]; then
                                echo `date`"  --  Waiting for all processes to clear "
                                wait
                                echo `date`"  --  All processes cleared!"
                                echo ""
                                nprocs=0
                            fi

                        done # ---- End of loop over trials
                    done # ---- End of loop over PESE flags
                done # ---- End of loop over number of obs
            done # ---- End of loop over all filters
        done # ---- End of loop over all obs error variances
    done # ---- End of loop over all ensemble sizes
done # ---- End of loop over all cycling intervals



echo `date`"  --  Waiting for all processes to clear "
wait
echo `date`"  --  All processes cleared!"
echo ""