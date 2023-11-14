'''
    SCRIPT TO PLOT VARIATIONS IN PROBIT-SPACE VARIANCE W.R.T. RANK HISTOGRAM PDF ENS SIZE
'''

import numpy as np
from matplotlib import use
use('agg')
import matplotlib.pyplot as plt
from scipy.stats import norm




# User controls
ens_size_list = np.arange(5,1000+1)


# place to store probit space variances
probit_variance_list = ens_size_list * 0.0


# Compute probit-space variances for every ens size
for isize in range(len(ens_size_list)):

    ens_size = ens_size_list[isize]

    # Compute RH distribution quantiles
    quantiles = (np.arange(ens_size)+1.)/(ens_size+1.)

    # Map quantiles to probit space
    probits = norm.ppf( quantiles )

    # Compute variance and store
    probit_variance_list[isize] = np.var( probits, ddof=1 )

# --- End of probit computaiton



# Plot out probit-space variances
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4,3))
ax.plot(ens_size_list,probit_variance_list,'-r')
ax.axhline(1., color='r',linestyle='--')
ax.set_ylabel('Un-adjusted\nProbit Space Variance')
ax.set_xlabel('Forecast Ensemble Size')
ax.set_title('Gaussian-tailed Rank Histogram\nUn-adjusted Probit Space Variances')
fig.subplots_adjust(left=0.2, bottom=0.2, top=0.8)
plt.savefig('EXPLORE_rank_histogram_probit_variance/probit_space_variances.svg')
