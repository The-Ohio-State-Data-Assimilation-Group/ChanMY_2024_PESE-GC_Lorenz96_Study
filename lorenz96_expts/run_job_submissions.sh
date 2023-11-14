#!/bin/bash


# Script to chain-submit a bunch of jobs with diff expt configs








# =======================================================================
# SECTION 1: FUNCTION DEFINITION
# =======================================================================




# FUNCTION TO RUN SPECIFIED EXPT SETTINGS ON SLURM
# LEVEL 1 FUNCTION (MOST RUDIMENTARY LEVEL)
# ------------------------------------------------
function submit_jobs () {

    # User inputs
    expt_set_name=$1
    obs_power=$2
    obs_err_var=$3
    filt_type=$4
    halfroi=$5
    pese_factor=$6
    marginal_distribution=$7

    # Hard-coded number of observations
    num_obs=40
    padded_num_obs=`printf "%03d" "$num_obs"`

    # Iterate over different settings
    for ens_size in 10 20 40 80; do
        for cyc_int in 3600 7200 10800 21600 43200; do

            # Padding everything that needs padding
            padded_filt_type=`printf "%02d" "$filt_type"`
            padded_ens_size=`printf "%03d" "$ens_size"`
            
            # Generate suffix for the config file
            fname_suffix=$expt_set_name
            fname_suffix=$fname_suffix"_obsPow"$obs_power
            fname_suffix=$fname_suffix"_obsVar"$obs_err_var
            fname_suffix=$fname_suffix"_obsNum"$padded_num_obs
            fname_suffix=$fname_suffix"_filt"$padded_filt_type
            fname_suffix=$fname_suffix"_ens"$padded_ens_size
            fname_suffix=$fname_suffix"_cycInt"$cyc_int
            fname_suffix=$fname_suffix"_peseFac"$pese_factor
            fname_suffix=$fname_suffix"_margDist"$marginal_distribution
            fname_suffix=$fname_suffix"_halfROI"$halfroi


            # Copy template config file
            config_fname="runtime_config_files/config."$fname_suffix
            cp TEMPLATE_config $config_fname

            # Write settings into config file
            sed -i "s|SED_HALF_ROI|$halfroi|g" $config_fname
            sed -i "s|SED_CYCLING_INTERVAL|$cyc_int|g" $config_fname
            sed -i "s|SED_ENS_SIZE|$ens_size|g" $config_fname
            sed -i "s|SED_OBS_NUM|$num_obs|g" $config_fname
            sed -i "s|SED_OBS_NAME|$expt_set_name|g" $config_fname
            sed -i "s|SED_OBS_POWER|$obs_power|g" $config_fname
            sed -i "s|SED_OBS_ERR_VAR|$obs_err_var|g" $config_fname
            sed -i "s|SED_FILTER_TYPE|$filt_type|g" $config_fname
            sed -i "s|SED_PESE_FACTOR|$pese_factor|g" $config_fname
            sed -i "s|SED_MARGINAL_DISTRIBUTION|$marginal_distribution|g" $config_fname

            # Copy template experiment script
            cp TEMPLATE_run_loop_over_expts.sh tmp.sh

            # Overwrite script with config file name
            sed -i "s|NAME_OF_CONFIG_FILE|$config_fname|g" tmp.sh

            # Submit experiment 
            job_name=$fname_suffix
            qsub -l select=1:ncpus=36:mpiprocs=36 -l walltime=12:00:00 -q economy -A $DAV_PROJECT -j oe -N $job_name tmp.sh                   
                
        done
    done

    #rm tmp.sh
}
# ----------------------------------------------- End of function to submit jobs











# FUNCTION TO LOOP OVER JOB SUBMISSIONS WITH DIFFERENT ROI VALUES 
# LEVEL 2 FUNCTION 
# ---------------------------------------------------------------
function run_roi_trials () {

    # User inputs
    expt_set_name=$1
    obs_power=$2
    obs_err_var=$3
    filt_type=$4
    pese_factor=$5
    marginal_distribution=$6


    # Localization half-radii to test
    halfroi_list="0.075"
    halfradius=0.075
    for ii in `seq 1 15`; do
        halfradius=`perl -e "printf  $halfradius * 1.3"`
        halfroi_list=$halfroi_list" "$halfradius
    done
    halfroi_list=$halfroi_list" "999.9

    
    # Submit jobs to test different ROIs
    for halfroi in $halfroi_list; do

        # Actual jobs submissions
        submit_jobs $expt_set_name $obs_power $obs_err_var $filt_type $halfroi $pese_factor $marginal_distribution

        # Jobs checker (need to prevent the number of jobs from exceeding 32)
        check_active_jobs 
    
    done 

}
# ----------------------------------------------- End of function to try ROIs









# FUNCTION TO RUN EXPERIMENTS WITH DIFFERENT PESE FACTORS
# LEVEL 3 FUNCTION (ENCAPSULATES ROI-TRIALS AND SUBMISSION)
function run_pese_factor_trials () {

    # User inputs
    expt_set_name=$1
    obs_power=$2
    obs_err_var=$3
    filt_type=$4
    marginal_distribution=$5

    # Anchoring path 
    anchor_path=`pwd`

    # Run with different PESE factors
    for pese_fac in 100; do #5 10 20; do

        # Run expts with specified pese factors
        run_roi_trials $expt_set_name $obs_power $obs_err_var $filt_type $pese_fac $marginal_distribution
        wait

        # Rename directory
        old_path=$anchor_path/all_expts/$expt_set_name/filter"$filt_type"
        new_path=$anchor_path/all_expts/$expt_set_name/filter"$filt_type"_pese"$pese_fac"X
        mv $old_path $new_path

    done


}
# ----------------------------------------------- End of functions to try PESE factors












# FUNCTION TO CHECK IF ALL JOBS ON SLURM ARE DONE
# UTILITY FUNCTION
# -----------------------------------------------
function check_active_jobs () {

    # Waiting for all existing jobs to clear
    flag_pause=true
    while $flag_pause; do 
        sleep 30
        njobs=`qstat -u manyau`
        njobs=`echo "$njobs" | wc -l`
        if [[ $njobs -gt 13 ]]; then
            flag_pause=true
            echo `date` ---- No. of active jobs: $njobs
        else
            flag_pause=false
            echo `date` ---- All jobs cleared. Proceeding with next round.
        fi
    done # ----- End of job-waiting loop
}
# ----------------------------------------------- End of function to check jobs

















































# =======================================================================
# SECTION 2: ACTUALLY RUN EXPERIMENTS
# =======================================================================


# SQRT obs with Gaussian PESE & linear regression
expt_set_name=SQRT_xGauss_yGauss
obs_power=0.5
obs_err=0.25
marginal_distribution=1

for filt_type in 1 2 8 ; do
    run_pese_factor_trials $expt_set_name $obs_power $obs_err $filt_type $marginal_distribution
done









# IDEN obs with Gaussian PESE & linear regression
expt_set_name=IDEN_xGauss_yGauss
obs_power=1.0
obs_err=1.0
marginal_distribution=1

for filt_type in 1 2 8 ; do
    run_pese_factor_trials $expt_set_name $obs_power $obs_err $filt_type $marginal_distribution
done



# IDEN obs with Gaussian PESE & linear regression
expt_set_name=SQUARE_xGauss_yGauss
obs_power=2.0
obs_err=16.0
marginal_distribution=1

for filt_type in 1 2 8 ; do
    run_pese_factor_trials $expt_set_name $obs_power $obs_err $filt_type $marginal_distribution
done



wait


