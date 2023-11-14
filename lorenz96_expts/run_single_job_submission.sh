#!/bin/bash

# Script to submit a single job to redo experiments

expt_set_name=$1
filt_type=$2
obs_err_var=$3
ens_size=$4
cyc_int=$5
num_obs=40


# Section to submit job
# ----------------------
# Generate fomatted suffix
padded_filt_type=`printf "%02d" "$filt_type"`
padded_ens_size=`printf "%03d" "$ens_size"`
padded_num_obs=`printf "%03d" "$num_obs"`
fname_suffix="filt"$padded_filt_type"_ens"$padded_ens_size"_obs"$padded_num_obs"_cycInt"$cyc_int

# Inject obs error info
fname_suffix=$fname_suffix"_obsVar"$obs_err_var


# Generate config script
config_fname="config."$expt_set_name"."$fname_suffix
cp TEMPLATE_config."$expt_set_name" $config_fname
sed -i "s|SED_FILTER_TYPE|$filt_type|g" $config_fname
sed -i "s|SED_ENS_SIZE|$ens_size|g" $config_fname
sed -i "s|SED_OBS_NUM|$num_obs|g" $config_fname
sed -i "s|SED_CYCLING_INTERVAL|$cyc_int|g" $config_fname

# TEMPORARY
sed -i "s|SED_OBS_ERR_VAR|$obs_err_var|g" $config_fname

# Copy template experiment script
cp TEMPLATE_run_loop_over_expts.sh tmp.sh

# Overwrite script with config file name
sed -i "s|NAME_OF_CONFIG_FILE|$config_fname|g" tmp.sh

# Submit experiment 
job_name=`echo $config_fname | cut -c 8-`
qsub -l select=1:ncpus=36:mpiprocs=36 -l walltime=12:00:00 -q economy -A $DAV_PROJECT -j oe -N $job_name tmp.sh


rm tmp.sh