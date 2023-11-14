#!/bin/bash

# Script to make namelist for DART
# ------------------------------------

path_to_template_namelist=$0/TEMPLATE_input.nml

# User controls
# --------------
ens_size=$1
use_pese=$2
pese_marginal=$3
ens_size_expanded=$4
localization_halfradius=$5
filter_type=$6
model_marginal_dist=$7
obs_marginal_dist=$8




# Write input.nml
# ----------------

# Copy over template namelist
cp $path_to_template_namelist  input.nml

# Put user inputs into namelist
sed -i "s|ENS_SIZE_DYNAMIC|$ens_size|g" input.nml
sed -i "s|USE_PESE|$use_pese|g" input.nml
sed -i "s|PESE_MARGINAL|$pese_marginal|g" input.nml
sed -i "s|ENS_SIZE_EXPANDED|$ens_size_expanded|g" input.nml
sed -i "s|LOCALIZATION_HALFRADIUS|$localization_halfradius|g" input.nml
sed -i "s|FILTER_TYPE|$filter_type|g" input.nml
sed -i "s|MODEL_DISTRIBUTION_TYPE|$model_marginal_dist|g" input.nml
sed -i "s|OBSERVATION_DISTRIBUTION_TYPE|$obs_marginal_dist|g" input.nml