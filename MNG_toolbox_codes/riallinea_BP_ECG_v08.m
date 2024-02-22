function [AAAt_ECGpeak222 AAAt_BPfoot222 AAAt_footIndex222]=riallinea_BP_ECG_v04(AAAt_ECGpeak,AAAt_BPfoot,AAAt_footIndex) 

% figure

AAAt_ECGpeak1=unique(AAAt_ECGpeak);
AAAt_BPfoot1=unique(AAAt_BPfoot);
AAAt_footIndex1=unique(AAAt_footIndex);

if AAAt_BPfoot1(1)>AAAt_ECGpeak1(2)
    AAAt_ECGpeak1(1)=[];
end

imin3=diff(AAAt_BPfoot1)<0.5;
AAAt_BPfoot1(imin3)=[]
AAAt_footIndex1(imin3)=[];

imin4=diff(AAAt_ECGpeak1)<0.5;
AAAt_ECGpeak1(imin4)=[]

for i=1:5
    
    if AAAt_BPfoot1(1)>AAAt_ECGpeak1(2)
        AAAt_ECGpeak1(1)=[];
    end

    if AAAt_BPfoot1(end)+0.3<AAAt_ECGpeak1(end-1)
        AAAt_ECGpeak1(end)=[];
    end

end
    
    
% subplot(311)
% for i=1:length(AAAt_ECGpeak1)
% SP3=AAAt_ECGpeak1(i); %your point goes here 
% line([SP3 SP3],[-50 50],'LineStyle','-.','LineWidth',0.7,'Color','r')%colorPlot{w},'LineStyle','--')
% hold on
% end
% 
% for i=1:length(AAAt_BPfoot1)
% SP3=AAAt_BPfoot1(i); %your point goes here 
% line([SP3 SP3],[-50 50],'LineStyle','-.','LineWidth',0.7,'Color','b')%colorPlot{w},'LineStyle','--')
% hold on
% end
% 
% xlim([0 AAAt_ECGpeak1(end)])

AAA_artifact=[];

length(AAAt_BPfoot1)

for i=1:min(length(AAAt_BPfoot1),length(AAAt_ECGpeak1))-1
    if abs(AAAt_BPfoot1(i)-AAAt_ECGpeak1(i))>= AAAt_ECGpeak1(i+1)-AAAt_ECGpeak1(i)
        AAA_artifact(i)=1;
    end
    
end


artifact=min(find(AAA_artifact==1))-1;

if artifact>0

    if length(artifact)>0 & artifact>=round(length(AAAt_ECGpeak1)/2)-4
    
        AAAt_ECGpeak222=AAAt_ECGpeak1(1:artifact);
        AAAt_BPfoot222=AAAt_BPfoot1(1:artifact);
        AAAt_footIndex222=AAAt_footIndex1(1:artifact);
        
        AAAt_ECGpeak222=[AAAt_ECGpeak222 AAAt_ECGpeak1(min(find(AAAt_ECGpeak1>AAAt_BPfoot222(end))))];
    
    end

    if length(artifact)>0 & artifact<round(length(AAAt_ECGpeak1)/2)-4
    
        AAAt_ECGpeak222=[AAAt_ECGpeak1(max(find(AAAt_ECGpeak1<AAAt_BPfoot1(artifact+4))):end)];
        AAAt_BPfoot222=AAAt_BPfoot1(artifact+4:end);
        AAAt_footIndex222=AAAt_footIndex1(artifact+4:end);
        
        if length(AAAt_BPfoot222)>=length(AAAt_ECGpeak222)
            AAAt_BPfoot222(end)=[];
            AAAt_footIndex222(end)=[];
        end
    
    end



% subplot(312)
% 
% for i=1:length(AAAt_ECGpeak222)
% SP3=AAAt_ECGpeak222(i); %your point goes here 
% line([SP3 SP3],[-50 50],'LineStyle','-.','LineWidth',0.7,'Color','r')%colorPlot{w},'LineStyle','--')
% hold on
% end
% 
% for i=1:length(AAAt_BPfoot222)
% SP3=AAAt_BPfoot222(i); %your point goes here 
% line([SP3 SP3],[-50 50],'LineStyle','-.','LineWidth',0.7,'Color','b')%colorPlot{w},'LineStyle','--')
% hold on
% end



% for i=1:length(AAAt_BPfoot222)%floor(length(t_CORTOpeak)/2
%     if AAAt_BPfoot222(1)<AAAt_ECGpeak222(1)
%         AAAt_BPfoot222(1)=[];
%        AAAt_footIndex222(1)=[];
%     end
%     
%     if AAAt_ECGpeak222(end)>AAAt_BPfoot222(end)
%         AAAt_ECGpeak222(end)=[];
%     end    
% end


% 
% xlim([0 AAAt_ECGpeak1(end)])
% 
% 
% subplot(313)
% 
% for i=1:length(AAAt_ECGpeak222)
% SP3=AAAt_ECGpeak222(i); %your point goes here 
% line([SP3 SP3],[-50 50],'LineStyle','-.','LineWidth',0.7,'Color','r')%colorPlot{w},'LineStyle','--')
% hold on
% end
% 
% for i=1:length(AAAt_BPfoot222)
% SP3=AAAt_BPfoot222(i); %your point goes here 
% line([SP3 SP3],[-50 50],'LineStyle','-.','LineWidth',0.7,'Color','b')%colorPlot{w},'LineStyle','--')
% hold on
% end
% 
% xlim([0 AAAt_ECGpeak1(end)])

    else
        AAAt_ECGpeak222=AAAt_ECGpeak1;
        AAAt_BPfoot222=AAAt_BPfoot1(1:end-1);
        AAAt_footIndex222=AAAt_footIndex1(1:end-1);
    end
    for i=1:5
        
        if AAAt_BPfoot222(1)>AAAt_ECGpeak222(2)
            AAAt_ECGpeak222(1)=[];
        end
        
        if AAAt_BPfoot222(end)<AAAt_ECGpeak222(end-1)
            AAAt_ECGpeak222(end)=[];
        end
    
    end

end