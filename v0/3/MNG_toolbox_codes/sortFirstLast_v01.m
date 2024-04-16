function [FISRT LAST] = sortFirstLast(t_ECG,t_SPIKES,ind_SPIKES)
% t_ECG=t_events,
% t_SPIKES,
% ind_SPIKES=1:length(t_SPIKES);



t_CORTOpeak=t_SPIKES(ind_SPIKES)';%ECG_RESULTS.t_events;  %BELT_Rpeak;

t_LUNGOpeak=t_ECG;%t_maxima(1:end);
N_ipol=50; %interpolation

%RR interval 
rrint=[];
tint=[];
for i=1:length(t_CORTOpeak)-1 
    rrint(i)=t_CORTOpeak(i+1)-t_CORTOpeak(i); 
    tint(i)=t_CORTOpeak(i+1);
end

 
Tcnt=1;
Tmat=[]; rrmat=[];
Tmat(:,:)=0; rrmat(:,:)=0; tmct=1; 

for cct1=1:length(t_LUNGOpeak)-1%:2:2:2:2
    for cct2=Tcnt:length(tint)
        if tint(cct2)<t_LUNGOpeak(cct1) 
         continue;
        elseif tint(cct2)>=t_LUNGOpeak(cct1)... 
         & tint(cct2)<t_LUNGOpeak(cct1+1)
            Tmat(tmct,cct1)=tint(cct2); 
            rrmat(tmct,cct1)=rrint(cct2); 
            tmct=tmct+1;
        else
        Tcnt=cct2; 
        tmct=1; 
        break;

        end
    end
end
 


for i=1:length(t_LUNGOpeak)-1
    try
        DTresp=(t_LUNGOpeak(i+1)-t_LUNGOpeak(i));
    for j=1:length(Tmat(:,i))
        radX(j,i)=(Tmat(j,i)-t_LUNGOpeak(i))/DTresp*2*pi;
        if  radX(j,i)<=0
            radX(j,i)=0;
        end
    end
    catch
        radX(j,i)=0;
    end
end

for i=1:length(t_LUNGOpeak)-1
    try
    
        DTresp=(t_LUNGOpeak(i+1)-t_LUNGOpeak(i));
    for j=1:length(Tmat(:,i))
        SpikeTimeX(j,i)=(Tmat(j,i)-t_LUNGOpeak(i));
        if  SpikeTimeX(j,i)<=0
            SpikeTimeX(j,i)=0;
        end
    end
    catch
        radX(j,i)=0;
    end
end

FISRT=Tmat(1,:)';
LAST=Tmat(end,:)'';

end


