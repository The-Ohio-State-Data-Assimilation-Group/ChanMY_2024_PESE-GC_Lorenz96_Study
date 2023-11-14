#!/bin/bash

# Run default settings (no PESE)
rm input.nml
ln -sf input_default.nml input.nml
./filter >& log.filter_default


# Run filter with PESE
rm input.nml
ln -sf input_pese.nml input.nml
./filter >& log.filter_pese