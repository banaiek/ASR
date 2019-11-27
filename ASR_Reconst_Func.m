%% SPike removal and integration in each frequency band
%output:
%LPF is the decomposed signal after removal
%Inputs:
%D_Dec is the dstLFP of each frequency band
%a1 &b1 the initial length
%a&b minimum limits
%mll & mlr are zero crossing pints in laft and rigth sides of trough
%SPT is spike in a particular trial
%LPF decomposed signal before removal
%fbW center of frequency band

%%----------
%written by Kia. Banaie Boroujeni 2018

function [LPF]=ASR_Reconst_Func(D_Dec,a1,b1,a,b,mll,mlr,SPT,LPF,fbw,Fs)  
    St=length(SPT); %Number of spikes
    fbw=round(fbw/2)+1; %half of the cycle
    tD_Dec=(D_Dec(a1-a-mll+1:a1+b+mlr+1)); % length adjusted average of decomposed 1st d-stLFP 

    for vSt=1:St
          %%trough alignment
          Dec=(LPF(SPT(vSt)-fbw:SPT(vSt)+fbw));
          [~,l_m]=min(Dec);
          l_m=(fbw+1)-l_m;
          td_Dec=(diff(LPF(SPT(vSt)-a1-l_m:SPT(vSt)+b1-l_m+1))/(1./Fs));%first derivative of each single spike trough aligned 
          d_Dec=td_Dec(a1-a-mll+1:a1+b+mlr+1); %length adjustment
          
          %integral of each individual and the average of 1st d-stLFP
          integ_M=[cumsum([LPF(SPT(vSt)-a-l_m-mll),tD_Dec*(1./Fs)])];
          integ_s=[cumsum([LPF(SPT(vSt)-a-l_m-mll),d_Dec*(1./Fs)])];
          cm2=abs(min((integ_s))); %trough value
          cm1=abs(min(integ_M)); %trough value of the average of 1st d-stLFP
          
          %Root mean square computation
          %integral of difference after normalization 
          integ_d=[cumsum([LPF(SPT(vSt)-a-l_m-mll),(1./Fs)*((d_Dec)/(cm2)-(tD_Dec)/(cm1))])];
          % RMS value of each individual and the average of 1st d-stLFP
          c1=sqrt(mean(integ_M.^2));
          c2=sqrt(mean(integ_s.^2));
          
          %RMS value of difference after normalization          
          c3=sqrt(mean(integ_d.^2));
          
          %RMS matched integral of the difference
          d_wf=cumsum([LPF(SPT(vSt)-a-l_m-mll),((d_Dec)/(cm2)-(tD_Dec)/(cm1))*(1./Fs)])*abs((c2-c1)/c3);
          d_wf(isnan(d_wf))= 0;

          LPF(SPT(vSt)-a-l_m-mll:SPT(vSt)+b-l_m+mlr+1)=d_wf; %Replacement in the decomposed portion of signal
      end
    end
    
