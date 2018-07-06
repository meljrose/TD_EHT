EHT_timedomain_variability
==============================
This uses python3.6!


Time domain variability analysis using EHT simulated observations

This code takes a hotspot simulation, a 'quiescent' flux simulation, and a 'noisy' simulation and combines them with variable ratios. It also adjusts the hotspot orbital period. It is meant to generate toy data for developing statistical tests that can retrieve the hotspot orbital period parameter and/or the 'quiescent' flux variabilty for different levels of noise. It uses the EHT-imaging package (https://github.com/achael/eht-imaging, commit f8cf49e) which has updates that are not backward compatible/are buggy (just a warning).  

EHT-imaging uses the ffmpeg backend (https://ffmpeg.org/) for their visuals, so I do too. You don't need this if you don't need exported movies.

Within EHT-imaging I replaced SITES.txt (eht-imaging/arrays/SITES.txt) with my own version that has differing SEFD values.

0_data_prep.ipynb needs to be run before data generation. The last few cells are examples of how to interact with the data object and are not required to run.

1_walkthrough.ipynb is not needed to generate data, but it has notes for the method.
