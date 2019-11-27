function out=ASR_Func(hp,in)
% hp contains some parameters
% in contains the entire data struct

Fs=hp.Fs;
b1=hp.b1ms*(Fs/1000); %length of extraction before the spike trough
a1=hp.a1ms*(Fs/1000); %length of extraction after the spike trough
Fpl=hp.Fpl;%lower limit for pspectrum to find the minimum frequency of decomposition
Fpu=hp.Fpu;%upper limit for pspectrum to find the minimum frequency of decomposition
L_Cutoff=hp.l_cutoff;
Hfr_dur=round(((1/Fpu)*(Fs))/2);
%% this method is supposed to remove spike-related transient from wideband data
%The example data set should be load directly which is a struct named in
%in.Spike_trial is the trial number that spikes were recorded
%in.Spike_cell is the corresponded to the trough of the spike (indecies of
%spikes which their time is in anoter parameter named in.Spike_time)
%in.data_wb is wide band exteracellularly recorded signal from one channel
%of recording
%in.trial_time: shows the time series for time of each trial
%This script finds the frequencies to decompose wideband data
%find zero crossing
%using function BP_filter, ASR_Decomp_Func, & ASR_Reconst_Func to respectively filter
%WB data for decomposition, decomposing spikes and find their normalized
%average of first derivative, & removing the average we found in ASR_Decomp_Func
%from each single spikes in ASR_Reconst_Func.
%trough alignment is done for each spike and in each decomposed frequency
%band
% zero crossing finds the length of removal for each frequency band average
% of first derivative of spikes


%%----------
%written by Kia. Banaie Boroujeni 2018


Spike_trial=in.Spike_trial; %trial number we have spike
Spike_cell=in.Spike_cell; %spike index (index refers to the spike time or spike trough)

%% computing spike triggered LFP average
STA_LFP=zeros(1,a1+b1+1); %spike triggered LFP average initialization.
vs=0; %spike counter
% these intervals can be changed, depending on the minimum frequency
% for artefact removal.
for st=1:length(Spike_trial)
    spike_wb=Spike_cell{st}; %spike indices of each trial
    for vr=1:length(Spike_cell{st})
        STA_LFP=STA_LFP+in.data_wb{Spike_trial(st)}(spike_wb(vr)-b1:spike_wb(vr)+a1);
        vs=vs+1;
    end
    
end
STA_LFP=STA_LFP/(vs); %Spike triggered LFP average

[Orig_pwr,Orig_fr]=pspectrum((STA_LFP),Fs,'FrequencyLimits',[Fpl Fpu]);
[~,Freq_l]=findpeaks((Orig_pwr.*Orig_fr),'SortStr','descend');   %Power/frequency peaks in the STA_LFP

Freq_l=Orig_fr(Freq_l)';
if isempty(Freq_l)
    temp=Fpl*sqrt(2);
else
    temp=Freq_l(1); %smallest starting frequency
end
frp=temp;
% construction of frequency array for filtering

while 1
    if sqrt(2)*temp>Fs/2
        break
    end
    temp=sqrt(2)*temp;
    frp=[frp,temp];
end

Freq_peaks=round(frp)'; %Decomposing Frequencies
display('Decomposing Frequencies are:');
if ~isempty(L_Cutoff)
    display('Decomposing Frequencies are:');
    Freq_peaks(Freq_peaks<=L_Cutoff)=[] % depending on how noisy is the...
    %data this limitation migth be set
else
    display('Decomposing Frequencies are:');
    
    Freq_peaks
end

fbw=round(Fs./Freq_peaks); %cycle length (for the center of each frequency band)


Steps=length(Freq_peaks)-1; % # of filters
Freq_seq=Freq_peaks(2:end); % Frequencies for filtering without the first frequency
n(1:Steps)=0;

for st=1:length(Spike_trial)
    [Lb,Hb]=BP_filter((Freq_peaks(1)),in.data_wb{Spike_trial(st)},Fs);%first low pass filter
    rYF{1,st}=Lb; %saving the first low passed signal for the final summation
    for i=1:Steps
        % finding the last step of filtering
        if i==Steps
            Lb=Hb; %the last high passed filtered signal
        else
            [Lb,Hb]=BP_filter(Freq_seq(i),Hb,Fs);
        end
        LbY{i,st}=Lb; %all band-passed signals & the last high-passed
        
        %% finding averageof 1st derivative for each decomposed signal
        
        if st==1
            
            mean_De=zeros(1,a1+b1+1);
            
        else
            mean_De= mean_De_f{1,i};
            
        end
        
        spike_wb=Spike_cell{st};
        for vr=1:length(Spike_cell{st}) %finding trough of each spike to match with LFP
            St=spike_wb(vr);
            [mean_De,n(i)]=ASR_Decomp_Func(Lb,St,b1,a1,n(i),mean_De,fbw(i),Fs);
        end
        mean_De_f{1,i}=mean_De;
        
    end
end

% Initilization of points to search for zero-crossing of each d_stLFP
for i=1:Steps
    [pval,ploc]=findpeaks((cumsum(mean_De_f{1, i}/n(i))));
    LocP=ploc(pval<mean(pval)+std(pval));
    loc_b=LocP((LocP-a1-1)<0);
    loc_a=LocP((LocP-a1-1)>0);
    cycle_Freq_b(i)=a1+1-loc_b(end);%limits from the left side
    cycle_Freq_a(i)=loc_a(1)-a1+1;%limits from the rigth side
end
%set the maximum limits not to pass our boundaries
fbw(fbw>100*(Fs/1000))=100*(Fs/1000);
cycle_Freq_b(cycle_Freq_b>100*(Fs/1000))=100*(Fs/1000);
cycle_Freq_a(cycle_Freq_a>100*(Fs/1000))=100*(Fs/1000);
display('Left init. points are:');
cycle_Freq_a
display('Rigth init. points are:');
cycle_Freq_b

dt=in;
%% finding zero-crossing for each band average 1st derivative
for i=1:Steps
    mean_De_f{1,i}=mean_De_f{1,i}/n(i);
    mean_De=(mean_De_f{1,i});
    mll1=find(fliplr((mean_De(1:b1-cycle_Freq_b(i)-1))).*(fliplr((mean_De(2:b1-cycle_Freq_b(i)))))<0);
    mlr1=find((mean_De(b1+2+cycle_Freq_a(i):end)).*(mean_De(b1+1+cycle_Freq_a(i):end-1))<0);
    mll(i)=mll1(1);
    mlr(i)=mlr1(1);
end

%% trough alignment & normalization & removing & integration
for st=1:length(Spike_trial)
    spike_wb=Spike_cell{st};
    yrF=0;
    
    for i=1:Steps
        LPF_Rmv=ASR_Reconst_Func(mean_De_f{1,i},b1,a1,cycle_Freq_a(i),cycle_Freq_b(i),mll(i),mlr(i),spike_wb,LbY{i,st},fbw(i),Fs);
        yrF=yrF+LPF_Rmv;
    end
    in.data_wb{Spike_trial(st)}=(yrF)+rYF{1,st}; % adding spike removed signal and the first low passed filter
    YRF{1,st}=yrF;
end


%% computing stLFP for raw and ASR
stLFP_ASR=0;
stLFP_raw=0;
ASR_rem=0;
%high_freq remenant removal
for st=1:length(Spike_trial)
    [~,c_yrem]=BP_filter(hp.Fpu,in.data_wb{Spike_trial(st)},Fs);
    spike_wb=Spike_cell{st};
    
    for vr=1:length(Spike_cell{st})
        ASR_rem=ASR_rem+c_yrem(spike_wb(vr)-Hfr_dur:spike_wb(vr)+Hfr_dur);
    end
end
ASR_rem=ASR_rem/vs;

for st=1:length(Spike_trial)
    c_yrF=in.data_wb{Spike_trial(st)};
    spike_wb=Spike_cell{st};
    
    for vr=1:length(Spike_cell{st})
        c_yrF(spike_wb(vr)-Hfr_dur:spike_wb(vr)+Hfr_dur)=c_yrF(spike_wb(vr)-Hfr_dur:spike_wb(vr)+Hfr_dur)-ASR_rem;
        stLFP_raw=stLFP_raw+dt.data_wb{Spike_trial(st)}(spike_wb(vr)-(225*(Fs/1000)):spike_wb(vr)+(225*(Fs/1000)));
        stLFP_ASR=stLFP_ASR+c_yrF(spike_wb(vr)-(225*(Fs/1000)):spike_wb(vr)+(225*(Fs/1000)));
    end
    in.data_wb{Spike_trial(st)}=c_yrF;
end

stLFP_ASR=stLFP_ASR/vs; %stLFP average after applying adaptive removal
stLFP_raw=stLFP_raw/vs;

out.stLFP_ASR=stLFP_ASR;
out.stLFP_raw=stLFP_raw;
out.data_wb=in.data_wb;
out.dt=dt;
return
