%% Adaptive Spike Removal method
%% this method is supposed to remove spike-related transient from Wideband data
%The example data set should be load directly which is a struct named in
%in.trial_time: shows the time series for time of each trial 
%in.Spike_trial is the trial number that spikes were recorded
%in.Spike_cell is the corresponded to the trough of the spike (indicies of
%spikes which their time is in another parameter named in.Spike_time)
%The length of the three varibles above is the number of spikes to be processed

%in.data_wb is wide band exteracellularly recorded signal from one channel
%of recording
%in.trial_time: shows the time series for time of each trial
%The length of the two varibles above is the number of trials that were recorded

%This script finds the frequencies to decompose Wideband data
%find zero crossing
%using function BP_filter, ASR_Decomp_Func, & ASR_Reconst_Func to respectively filter
%WB data for decomposition, decomposing spikes and find their
%average of first derivative, & removing the average we found in ASR_Decomp_Func
%from each single spikes in ASR_Reconst_Func.
%trough alignment  is done for each spike and in each decomposed frequency
%band
% zero crossing finds the length of removal for each frequency band average
% of first derivative of spikes
% The code needs the MATLAB signal processing toolbox for function (pspectrum).

%%----------
%written by Kia. Banaie Boroujeni 2018
clc
close all
clear all

load('Dataset_example.mat') %loding the dataset example

 hp.Fs=32000; %sampling rate Hz
 hp.a1ms=400; %time (msec.) of extraction after the spike trough
 hp.b1ms=400; %time (msec.) of extraction before the spike trough
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
          plot([(1:length(out.stLFP_raw))/(hp.Fs/1000)-find(out.stLFP_raw==min(out.stLFP_raw))/(hp.Fs/1000)],out.stLFP_ASR,'r','LineWidth',1)
          hold on;
          plot([(1:length(out.stLFP_raw))/(hp.Fs/1000)-find(out.stLFP_raw==min(out.stLFP_raw))/(hp.Fs/1000)],out.stLFP_raw,'b','LineWidth',1)

          xlabel({'Time (mS)'});
          ylabel({'Amplitude'});