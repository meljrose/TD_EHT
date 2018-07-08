EHT_timedomain_variability
==============================
This uses python3.6 and jupyter notebooks. 

Package requirements:

jupyter notebook (in conda, or can pip install)

eht-imaging (https://github.com/achael/eht-imaging)

h5py (can pip install)

emcee (can pip install)

ffmpeg (optional) (https://ffmpeg.org/)

==============================

Time domain variability analysis using EHT simulated observations

This code takes a hotspot simulation, a 'quiescent' flux simulation, and a 'noisy' simulation and combines them with variable ratios. It also adjusts the hotspot orbital period. It is meant to generate toy data for developing statistical tests that can retrieve the hotspot orbital period parameter and/or the 'quiescent' flux variabilty for different levels of noise. It uses the EHT-imaging package (https://github.com/achael/eht-imaging, specifically commit f8cf49e) which has updates that are not backward compatible/are buggy (just a warning).  

Installing eht-imaging: 

Once you've cloned eht-imaging, you can use "git checkout f8cf49e397c99e04de971b29e69270fee76dacdf . " to roll back to this commit. 
This version of eht-imaging requires "emcee" which you can get via pip. You'll get a warning that you haven't installed "NFFT" but that isn't required for these notebooks. 

EHT-imaging uses the ffmpeg backend (https://ffmpeg.org/) for their visuals, so I do too. You don't need this if you don't need exported movies but you won't see any of the animations in the notebook without it. 

Within EHT-imaging I replaced SITES.txt (eht-imaging/arrays/SITES.txt) with my own version that has differing SEFD values.

0_data_prep.ipynb needs to be run before data generation. The last few cells are examples of how to interact with the data object and are not required to run.

1_walkthrough.ipynb is not needed to generate data, but it has notes for the method.
