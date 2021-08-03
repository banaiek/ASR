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
%adjusted on Aug. 2021

function [LPF]=ASR_Reconst_Func(D_Dec,a1,b1,a,b,mll,mlr,SPT,LPF,fbw,Fs)  
    St=length(SPT); %Number of spikes
    fbw=round(fbw/2)+1; %half of the cycle
    L=length(LPF);
          
          
    for vSt=1:St
          %%trough alignment

           l1=min(fbw,SPT(vSt)-1);
            l2=min(fbw,L-SPT(vSt));
            [~,l_m]=min(LPF(SPT(vSt)-l1:SPT(vSt)+l2));
 
           SPT(vSt) = SPT(vSt)-l1+l_m-1;
 
%          l1=max(SPT(vSt)-fbw,1);
%          l2=min(SPT(vSt)+fbw,L);
%          
%           Dec=(LPF(l1:l2));
%           [~,l_m]=min(Dec);
%           l_m = SPT(vSt)-l1-l_m+1;
%           
%          SPT(vSt)=SPT(vSt)-l_m;

         ls1=min([b+mll,b1,SPT(vSt)-2]);
         ls2=min([a+mlr,a1,L-SPT(vSt)-1]);

         fl1=min([b1,SPT(vSt)-2]);
         fl2=min([a1,L-SPT(vSt)-1]);
         
          tD_Dec=(D_Dec(b1-ls1:b1+ls2)); % length adjusted average of decomposed 1st d-stLFP 
          td_Dec=(diff(LPF(SPT(vSt)-fl1-1:SPT(vSt)+fl2))/(1./Fs));%first derivative of each single spike trough aligned 
          

          
          d_Dec=(diff(LPF(SPT(vSt)-ls1-1:SPT(vSt)+ls2))/(1./Fs));
          %integral of each individual and the average of 1st d-stLFP
          integ_M=[cumsum([LPF(SPT(vSt)-ls1-1),tD_Dec*(1./Fs)])];
          integ_s=[cumsum([LPF(SPT(vSt)-ls1-1),d_Dec*(1./Fs)])];
          cm2=abs(min((integ_s))); %trough value
          cm1=abs(min(integ_M)); %trough value of the average of 1st d-stLFP
          
          %Root mean square computation
          %integral of difference after normalization 
          integ_d=[cumsum([LPF(SPT(vSt)-ls1-1),(1./Fs)*((d_Dec)/(cm2)-(tD_Dec)/(cm1))])];
          % RMS value of each individual and the average of 1st d-stLFP
          c1=sqrt(mean(integ_M.^2));
          c2=sqrt(mean(integ_s.^2));
          
          %RMS value of difference after normalization          
          c3=sqrt(mean(integ_d.^2));
          
          %RMS matched integral of the difference
          d_wf=cumsum([LPF(SPT(vSt)-ls1-1),((d_Dec)/(cm2)-(tD_Dec)/(cm1))*(1./Fs)])*abs((c2-c1)/c3);
          d_wf(isnan(d_wf))= 0;

          LPF(SPT(vSt)-ls1:SPT(vSt)+ls2+1)=d_wf; %Replacement in the decomposed portion of signal
      end
    end
    
