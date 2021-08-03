%% Adaptive Spike Removal method
%% this method is supposed to remove spike-related transient from Wideband data
%The example data set should be load directly which is a struct named "in".


%in.data_wb: is wide band exteracellularly recorded signal from one channel
%of recording, each cell corresponds to 1 trial of wideband signal.

%in.trial_time: shows the time series for time of each trial (corresponds
%to the timeseries of wideband data)

%in.Spike_cell: is a cell structure where each cell shows spikes in one trial which
%has a column of values corresponding to the location of the spike troughs on
%the wideband data timeseries in that trial (it is location indice not time which time
%itself is in.Spike_time(unused here))).

%in.Spike_trial: is the trial number corresponding to the cells in.Spike_cell.
%(for example, if the first trial that we have spike values is trial 10 on
%the in.data_wb then the in.Spike_trial(1)=10.


%The length of the three varibles above is the number of spikes to be processed




%This script finds the frequencies to decompose Wideband data
%using function BP_filter, ASR_Decomp_Func, & ASR_Reconst_Func to respectively filter
%WB data for decomposition, decomposing spikes and find their
%average of first derivative, & removing the average we found in ASR_Decomp_Func
%from each single spikes in ASR_Reconst_Func.
%trough alignment  is done for each spike and in each decomposed frequency
%band
% zero crossing finds the length of removal for each frequency band average
% of first derivative of spikes
% The code needs the MATLAB signal processing toolbox for function (pspectrum).

%Inputs:
% hp contains parameters to be specified by user
% in contains the entire data struct
%Outputs
% out.stLFP_ASR : Spike-triggered average after ASR
% out.stLFP_raw : Spike-triggered average on raw wideband signal
% out.data_wb   : widband signal with spikes removed
% out.Orig_data : original data-structure

%%----------
%written by Kia. Banaie Boroujeni 2018
%adjusted on Aug. 2021

clc
close all
clear all

load('Dataset_example.mat') %loding the dataset example

hp.Fs=32000; %sampling rate Hz
hp.a1ms=400; %time (msec.) of extraction after the spike trough
hp.b1ms=400; %time (msec.) of extraction before the spike trough
hp.StDur=500; %time (msec.) of spike triggered average
hp.Fpl=2;%lower limit for pspectrum to find the minimum frequency of decomposition
hp.Fpu=200;%upper limit for pspectrum to find the minimum frequency of decomposition
hp.MAIa= 5; %miss-alignment time interval (msec.) after the current spike indices
hp.MAIb= 5; %miss-alignment time interval (msec.) before the current spike indices
hp.l_cutoff=2; %setting a cut off frequency for removing low freq. decomposition depending on how noisy the data is
%% Spike trough alignment Function
[hp,in]=Spike_align(hp,in);
%% Spike removal Function
tic;  out=ASR_Func(hp,in); toc;


%% Plotting out.stLFP for raw and ASR
figure,
t=[(1:length(out.stLFP_raw))/(hp.Fs/1000)-find(out.stLFP_raw==min(out.stLFP_raw))/(hp.Fs/1000)];

plot(t,out.stLFP_ASR,'r','LineWidth',1)
hold on;
plot(t,out.stLFP_raw,'b','LineWidth',1)

xlabel({'Time (mS)'});
ylabel({'Amplitude'});
