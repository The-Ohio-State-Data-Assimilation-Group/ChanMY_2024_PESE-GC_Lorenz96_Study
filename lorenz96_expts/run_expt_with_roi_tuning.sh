#!/bin/bash

# ================================================
# SCRIPT TO SET UP AND RUN ROI TUNING EXPERIMENTS
# ================================================


# READ IN SETTINGS FROM COMMAND LINE
CONFIG_FILE=$1
filter_type=$2
n_obs=$3
ens_size=$4
obs_err_var=$5
cyc_interval_seconds=$6
use_pese=$7
trial=$8


# LOAD CONFIG FILE
. $CONFIG_FILE


# Number-padded trial id
trial_padded=`printf %03d $trial`


# =====================================
# RUN TESTS WITH DIFFERENT ROI VALUES
# -------------------------------------

# Access directory of the trial
trial_dir=$dir_expts/$obs_name"/filter"$filter_type"/ObsErrVar"$obs_err_var"/CycInt"$cyc_interval_seconds"/EnsSize"$ens_size"/ObsNum"$n_obs"/PESE"$use_pese"/Trial"$trial_padded
mkdir -p $trial_dir
cd $trial_dir




# Loop over all ROIs
for half_roi in $loc_halfradius_list; do 
    
    echo ""
    echo "RUNNING TEST FOR HALF-ROI: "$half_roi

    # Directory name to hold current tuning experiment
    expt_name="halfROI"$half_roi
    dir_expt=$expt_name   #$dir_expts/$expt_name

    # Make and setup directory for experiment
    if [[ -e $dir_expt ]]; then
        rm -r $dir_expt
    fi
    mkdir -p $dir_expt
    cd $dir_expt

    # Link in useful stuff from DART compilation
    ln -s $dir_dart/create_obs_sequence .
    ln -s $dir_dart/create_fixed_network_seq .
    ln -s $dir_dart/perfect_model_obs .
    ln -s $dir_dart/filter .
    cp $dir_dart/filter_input_list.txt .
    cp $dir_dart/filter_output_list.txt .

    # copy in useful bash scripts
    cp $dir_bash/run_create_obs_sequence.sh .
    cp $dir_bash/run_create_fixed_network_seq.sh .
    cp $dir_bash/setup_namelist.sh .

    # Construct namelist
    pese_size=$(( $ens_size*$pese_ens_factor ))
    ./setup_namelist.sh $ens_size $use_pese 1 $pese_size $half_roi $filter_type $model_marginal_dist $obs_marginal_dist

    # Symlink in the netcdf files
    ln -s $dir_init_conditions/id$trial_padded"_filter_input.nc" filter_input.nc
    ln -s $dir_init_conditions/id$trial_padded"_perfect_input.nc" perfect_input.nc

    # Create observation network
    bash ./run_create_obs_sequence.sh $n_obs $obs_type $obs_power $obs_err_var >& log.create_obs_sequence
    bash ./run_create_fixed_network_seq.sh $n_cycle_full $cyc_interval_seconds >& log.create_fixed_network_seq

    echo FINISHED SETTING UP DIRECTORY




    echo CREATE OBSERVATION SEQUENCE
    # Generate observation sequence file
    ./perfect_model_obs >& log.perfect_model_obs  
    echo FINISHED CREATING OBSERVATION SEQUENCE



    echo RUN DA EXPERIMENT
    ./filter >& log.filter
    
    # Remove obs_seq files to save space
    rm obs_seq.*

    # Compress log files to save space
    zip log_files.zip log.* dart_log*

    # Save first 500 lines and last 500 lines of log.filter 
    head -n 500 log.filter > short.filter.log
    tail -n 500 log.filter >> short.filter.log

    # Remove unneeded log files and netcdf files
    rm log.* dart_log*
    
    echo FINISHED RUNNING DA EXPERIMENT

    
    echo COMPUTE USEFUL STATS AND REMOVE NCFILES
    python -u $dir_bash/determine_useful_stats.py
    rm *nc
    

    cd ..


done
# ------ End of loop over half-ROIs




# echo ""
# echo ""
# echo "DETERMINING BEST ROI IN THE SET OF ROIs TESTED"



# # ==============================
# # DETERMINE BEST PERFORMING ROI
# # ------------------------------
# all_roi_list=""
# for roi in $loc_halfradius_list; do
#     all_roi_list=$roi,$all_roi_list
# done

# best_dir=`python -u $dir_bash/determine_optimal_roi.py $all_roi_list`
# ln -s $best_dir long

# echo ""