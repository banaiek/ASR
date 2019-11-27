# ASR
Adaptive Spike Removal from LFP

Spikes introduce artificial phase lag to the LFP when they are recorded from the same channel. The scripts implement the adaptive spike removal method (ASR) that succeeds removing spike current leakage to the LFP across a broad range of frequencies. This study is published in J. Neuroscience Methods https://doi.org/10.1016/j.jneumeth.2019.108485

%% The data set is named Dataset_example.mat should be load directly into the Matlab path. % There are 3 functions and one main script in the folder. Adaptive Spike Removal (ASR) method %% ASR method is supposed to remove spike-related transient from Wideband data

The dataset contains parameters listed below:

%in.Spike_trial is the trial number that spikes were recorded %in.Spike_cell is the corresponded to the trough of the spike (indices of %spikes which their time is in another parameter named in.Spike_time) %in.data_wb is wide band extracellularly recorded signal from one channel of recording %in.trial_time: shows the time series for time of each trial

%This main script finds frequencies to decompose Wideband data and finds zero crossing. %Using function BP_filter, ASR_Decomp_Func, & ASR_Reconst_Func to respectively filter WB data for decomposition, decomposing spikes and find average of their first derivative, & removing the average we found in ASR_Decomp_Func from each single spikes in ASR_Reconst_Func.

Trough alignment is done for each spike and in each decomposed frequency band zero crossing finds the length of removal for each frequency band average of first derivative of spikes.

To download a Dataset_Example you can visit: http://accl.psy.vanderbilt.edu/resources/adaptive-spike-removal-code/

Added Nov. 2019 #Primary spike trough alignment is added in case that “in.Spike_cell” (which is the indices of spikes aligned into their trough) contains spikes miss-aligned to their trough.

Copyright (C) 2019, Kianoush Banaie Boroujeni.
