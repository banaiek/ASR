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

function [mean_De_f,n]=ASR_Decomp_Func(De_f,St,a1,b1,n,mean_De,fbw,Fs)
%% trough alignment
 fbw=round(fbw/2)+1;
 [~,l_m]=min(De_f(St-fbw:St+fbw));
 l_m=(fbw+1)-l_m;

%differentiation 
 d_De_f=(diff(De_f(St-a1-l_m:St+b1-l_m+1))/(1/Fs));

if ~isnan(d_De_f)
 mean_De_f=mean_De+d_De_f; %sum of  decomposed spikes
n=n+1; %counter
else    
 mean_De_f=mean_De;   
end

end
