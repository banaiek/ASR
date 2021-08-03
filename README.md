# ASR
Adaptive Spike Removal from LFP

Spikes introduce artificial phase lag to the LFP when they are recorded from the same channel. The scripts implement the adaptive spike removal method (ASR) that succeeds removing spike current leakage to the LFP across a broad range of frequencies. This study is published in J. Neuroscience Methods https://doi.org/10.1016/j.jneumeth.2019.108485

%% The data set is named Dataset_example.mat should be load directly into the Matlab path. % There are 4 functions and one main script in the folder. Adaptive Spike Removal (ASR) method %% ASR method is supposed to remove spike-related transient from Wideband data

The dataset contains an structure of cells and doubles listed below:

%in.data_wb (cell): is wide band exteracellularly recorded signal from one channel of recording, each cell corresponds to 1 trial of wideband signal.

%in.trial_time (cell): shows the time series for time of each trial (corresponds to the timeseries of wideband data)

%in.Spike_cell (cell): is a cell structure where each cell shows spikes in one trial which has a column of values corresponding to the location of the spike troughs on the wideband data timeseries in that trial (it is location indice not time which time itself is in.Spike_time(unused here))).

%in.Spike_trial (double): is the trial number corresponding to the cells in.Spike_cell. (for example, if the first trial that we have spike values is trial 10 on the in.data_wb then the in.Spike_trial(1)=10.


%The main script finds frequencies to decompose Wideband data and finds zero crossing. %Using function BP_filter, ASR_Decomp_Func, & ASR_Reconst_Func to respectively filter WB data for decomposition, decomposing spikes and find average of their first derivative, & removing the average we found in ASR_Decomp_Func from each single spikes in ASR_Reconst_Func.

Edit: Nov. 2019 #Primary spike trough alignment is added in case that “in.Spike_cell” (which is the indices of spikes aligned into their trough) contains spikes miss-aligned to their trough.

Edit: Aug. 2021: Removing the constraint on the margin needed for removal. Still a sanity check for not having a spike indice on the first or the last few (e.g. 5) data points (which in practice and in a fine dataset should not) of a wideband data timeseries prevents any potentional problem.
This Edit was very extensive and resulted in improvment and ease of using the method



Trough alignment is done for each spike and in each decomposed frequency band zero crossing finds the length of removal for each frequency band average of first derivative of spikes.

To download a Dataset_Example you can visit: http://accl.psy.vanderbilt.edu/resources/adaptive-spike-removal-code/



Copyright (C) 2019, Kianoush Banaie Boroujeni.
