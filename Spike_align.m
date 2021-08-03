%%This function is added to align spikes into their trough before feeding the ASR algorithm.
%%hp.MAIa is the interval for looking toward a local minima after the
%%current spike time onset
%%hp.MAIb is the interval for looking toward a local minima before the
%%current spike time onset
% the miss-aligned spike onsets will be saved in a cell structure,
% named"in.Spike_cell_missaligned", with the same data format as
% in.Spike_cell.

%written by Kia. Banaie Boroujeni (Nov. 2019)
%adjusted on Aug. 2021

function [hp,in]=Spike_align(hp,in)

Fs=hp.Fs;
mai_a=round(hp.MAIa*(Fs/1000)); %after spike indicies 
mai_b=round(hp.MAIb*(Fs/1000)); %before spike indicies 
Spike_trial=in.Spike_trial; %trial number we have spike
Spike_cell=in.Spike_cell; %spike index (index refers to the spike time or spike trough)

%% Trough Alignment
     
    for st=1:length(Spike_trial)
      spike_wb=Spike_cell{st}; %spike indices of each trial
      y=in.data_wb{Spike_trial(st)};
      L=length(y);
      for vr=1:length(Spike_cell{st})
          sp1 = min(spike_wb(vr)-1,mai_b); % borders around spikes
          sp2 = min(mai_a,L-spike_wb(vr)); 
          [~,loc]=min(in.data_wb{Spike_trial(st)}(spike_wb(vr)-sp1:spike_wb(vr)+sp2));
         in.Spike_cell{st}(vr,1)=spike_wb(vr)+loc-sp1-1; 
      end
      
 
    end
    in.Spike_cell_missaligned=Spike_cell;
end
