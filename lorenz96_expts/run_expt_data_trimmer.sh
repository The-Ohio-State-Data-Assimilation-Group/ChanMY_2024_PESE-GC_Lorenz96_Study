#!/bin/bash

# ===========================================================
# SCRIPT TO REMOVE UNNEEDED FILES AND COMPRESS LOG FILES
# ===========================================================



CONFIG_FILE=$1

. $CONFIG_FILE


# Useful stuff for parallelization
nprocs_total=4

nprocs_active=0




# Function to trim down on log files
function trim_files {

    
    # # Adding dart_log files to log_files
    # zip -ur log_files.zip dart_log.nml dart_log.out  

    # rm dart_log.nml dart_log.out

    # Compressing log files
    zip log_files.zip log.* dart_log*

    # Save first 500 lines and last 500 lines of log.filter 
    head -n 500 log.filter > short.filter.log
    tail -n 500 log.filter >> short.filter.log

    # Remove unneeded log files
    rm log.* dart_log*

    # # Remove unneeded nc files
    # if [[ -e filter_input.nc ]]; then
    #     rm filter_input.nc filter_output.nc
    # fi

}
export -f trim_files





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

                            # Determine directory trial name
                            trial_padded=`printf %03d $trial`
                            expt_name=$obs_name"/filter"$filter_type"/ObsErrVar"$obs_err_var"/CycInt"$cyc_interval_seconds"/EnsSize"$ens_size"/ObsNum"$n_obs"/PESE"$use_pese"/Trial"$trial_padded

                            date
                            echo Trimming $expt_name space usage
                            echo ""
                            expt_name=$dir_expts/$expt_name

                            # Removing all forecast.nc files from tuning runs
                            for sub_expt in `ls --color=none -d $expt_name/halfROI*`; do
                                ncdump -h $sub_expt/forecast.nc >& $sub_expt/log.ncdump_forecast
                                rm $sub_expt/forecast.nc  &
                            done

                            
                            # Compressing all log files
                            for sub_expt in `ls --color=none -d $expt_name/halfROI*` $expt_name/long; do
                                
                                cd $sub_expt

                                if [[ -e log.filter ]]; then

                                    # Run file trimming in background
                                    trim_files &
                                    
                                    # Increment number of active processes
                                    nprocs_active=$(($nprocs_active+1))
                                    if [[ $nprocs_active -ge $nprocs_total ]]; then
                                        wait
                                        nprocs_active=0
                                    fi

                                fi

                            done # --- End of loop over sub directories


                        done # ---- End of loop over trials
                    done # ---- End of loop over PESE flags
                done # ---- End of loop over number of obs
            done # ---- End of loop over all filters
        done # ---- End of loop over all obs error variances
    done # ---- End of loop over all ensemble sizes
done # ---- End of loop over all cycling intervals


wait
