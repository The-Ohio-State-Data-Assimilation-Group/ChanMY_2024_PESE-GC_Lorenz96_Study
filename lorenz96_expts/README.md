# Directory of scripts regulating Lorenz 96 experiments

Written by: Man-Yau (Joseph) Chan



## Instructions to replicate experiments in Chan's manuscript ("Improving EnsDA with PESE-GC")


### Important notes
This collection of experiments is computationally intensive to replicate! You will need to run this collection on a computing cluster.

If you have difficulty replicating these experiments, please contact Joseph Chan at `chan.1063@osu.edu`.

&nbsp; &nbsp;


### Step 1: Compile DART
Compile DART with PESE-GC (`../DART_w_PESE-GC`) for the Lorenz 1996 model. Compilation instructions can be found on the [DART documentation](https://docs.dart.ucar.edu/en/latest/).

Note that you will need to compile DART with either Intel's Math Kernel Library (MKL) or LAPACK. PESE-GC relies on matrix operations found in MKL/LAPACK.

&nbsp; &nbsp;



### Step 2: Adjust `run_job_submissions.sh` to work with your cluster

Currently, `run_job_submissions.sh` is designed to work with NCAR's Cheyenne supercomputer.

You will need to replace all the `qsub` and `qstat` commands with relevant commands for your cluster.

Note that each configuration (see manuscript for definition) is trialed 36 times! So, you will need to use 36 cores for each configuration.

&nbsp; &nbsp;




### Step 3: Compute RMSEs for the various experiments

For experiments with IDEN observations, run
```
bash run_evaluate_setup_stats.sh config.IDEN_consolidated >> performance_logs/IDEN_performance.txt
```

&nbsp;

For experiments with SQRT observations, run
```
bash run_evaluate_setup_stats.sh config.SQRT_consolidated >> performance_logs/SQRT_performance.txt
```

&nbsp;

For experiments with SQUARE observations, run
```
bash run_evaluate_setup_stats.sh config.SQUARE_consolidated >> performance_logs/SQUARE_performance.txt
```

&nbsp; &nbsp;



### Step 4: Plot impacts of PESE-GC on RMSEs

For experiments with IDEN observations, run
```
python plot_impacts_of_pese_on_rmse.py performance_logs/IDEN_performance.txt
```

&nbsp;

For experiments with SQRT observations, run
```
python plot_impacts_of_pese_on_rmse.py performance_logs/SQRT_performance.txt
```

&nbsp;

For experiments with SQUARE observations, run
```
python plot_impacts_of_pese_on_rmse.py performance_logs/SQUARE_performance.txt
```


