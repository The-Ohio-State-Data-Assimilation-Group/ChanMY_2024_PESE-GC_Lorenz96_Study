#!/bin/bash

# Script to check if PESE is working correctly

# Run default stuff
rm input.nml
ln -sf input_default.nml input.nml
./filter >& log.filter_default

# Run PESE stuff
rm input.nml
ln -sf input_pese.nml input.nml
./filter >& log.filter_pese