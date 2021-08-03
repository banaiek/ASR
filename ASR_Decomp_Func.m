%% This function gets each decomposed trial and find the average of 1st drivative decomposed and trough aligned spikes
%outputs: 
%mean_De_f :  average of first derivative of
%decomposed spikes
%of spikes in the trial
%n is number of decomposed spikes in the frequency band (it migth vary
%frequency to frequency depend on whether there is NaNs.
%Inputs:
%De_f: decomposed signal
%st: spike trough
%a1,b1: initial length
%n:counter of number of spikes from previous trials
%mean_De :average of first derivative of
%decomposed spikes in previous trials
%fbw : frequency band 
%%----------
%written by Kia. Banaie Boroujeni 2018
%adjusted on Aug. 2021

function [mean_De_f]=ASR_Decomp_Func(De_f,St,a1,b1,fbw,Fs)
%% trough alignment
 mean_De_f = nan(1,a1+b1+1);
 fbw2=round(fbw/2)+1;
 L=length(De_f);
 l1=min(fbw2,St-1);
 l2=min(fbw2,L-St);
 [~,l_m]=min(De_f(St-l1:St+l2));
 
 St = St-l1+l_m-1;
 
 fl1 = min(St-2,b1);
 fl2 = min(a1,L-St);
 
 ml1 = b1+1-fl1;
 ml2 = b1+fl2+1;
    

%differentiation 
 d_De_f = (diff(De_f(St-fl1-1:St+fl2))/(1/Fs));


 mean_De_f(1,ml1:ml2) = d_De_f; %sum of  decomposed spikes


end
