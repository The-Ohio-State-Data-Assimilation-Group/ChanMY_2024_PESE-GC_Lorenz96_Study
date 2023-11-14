#!/bin/bash

# Simple script to feed info into create_obs_sequence


# USER SETTINGS
# --------------

# Number of cycles
nCycle=$1

# Interval between cycles (in seconds)
cyc_interval_seconds=$2


# ACTUAL CODE
# ------------

# Generate command prompt feeder
inputs="\n 1 \n $nCycle \n 50 0 \n 0 $cyc_interval_seconds \n \n"

printf "$inputs" | ./create_fixed_network_seq