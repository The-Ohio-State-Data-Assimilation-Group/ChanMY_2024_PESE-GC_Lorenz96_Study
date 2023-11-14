#!/bin/bash

# Simple script to feed info into create_obs_sequence


# USER SETTINGS
# --------------

# Number of observations
nObs=$1

# Observation type
obs_type=$2             # 1 for IDEN, 5 for power obs

# Observation power
obs_pow=$3

# Observation error variance
obs_err_variance=$4







# ACTUAL CODE
# ------------

# Generate command prompt feeder
inputs="$nObs \n 0 \n 0 \n"        # ---- Number of obs, number of data copies, number of QC values
for iob in `seq 1 $nObs`; do

    #obs_pos=`perl -e "printf $((iob-1))*1./$nObs"`

    if [[ $obs_type == 5 ]]; then
        inputs="$inputs 0 \n 5 \n $obs_pow \n -1 \n 0 0 \n $obs_err_variance \n"
    fi

    if [[ $obs_type == 1 ]]; then
        inputs="$inputs 0 \n 1 \n -1 \n 0 0 \n $obs_err_variance \n"
    fi

done
inputs="$inputs \n "

# Feed inputs into code
printf "$inputs" | ./create_obs_sequence 
