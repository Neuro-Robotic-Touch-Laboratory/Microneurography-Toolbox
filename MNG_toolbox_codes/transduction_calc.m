function [RRR7, RRR8 ] = transduction_calc(data, t_events_ecg, t_bpFoot, bpValues, t_spikes, footIndex)

% GDAplot_DynamicIndexes_v26(namefileXXX(1:end-9), t_spikes,                  num2str(i),      'ecg',                t_bpFoot*1000,                        colorVec(i,:), thre);
% plottingPhasesGDA(         nome_file,            t_spikes,                  ref,             rif,                  footIndex,                            colorBG,       thre) 
%                            for loading files     spiketimes of the cluster  clusternumber?   label for plotting ?  foot indx? from BP(check if fs =1khz)                3.5 
% 
% load([nome_file '_msna']) %raw msna only data (1 x N) 
% load([nome_file '_ecg'])% rpeaks timestamps (t_events_ecg)
% 
% load([nome_file '_bpSyst']) %% timestamps  systolic peaks (t_bpSyst)
% load([nome_file '_bpFoot']) %% timestamps of foot idx (t_bpFoot)
% load([nome_file '_bpValues']) %% bloodpressure data (t_bpFoot)

AAAt_BPfoot=t_bpFoot;
AAAt_footIndex=footIndex;
AAAt_ECGpeak=t_events_ecg';

for i = 1 : length(AAAt_BPfoot)
    if AAAt_BPfoot(1) < AAAt_ECGpeak(1)
        AAAt_BPfoot(1) = [];
        AAAt_footIndex(1) = [];
    end
   
end

[t_ecgCUT, t_bpCUT, footCUT] = riallinea_BP_ECG_v08(AAAt_ECGpeak, AAAt_BPfoot, AAAt_footIndex);

AAAt_ECGpeak = [];
AAAt_ECGpeak = t_ecgCUT;

AAAt_BPfoot = [];
AAAt_BPfoot = t_bpCUT;

AAAt_footIndex = [];
AAAt_footIndex = footCUT;

% if rif == 'ecg'
    t_events = AAAt_ECGpeak;
% end
% 
% if rif == 'bpF'
%     t_events = AAAt_BPfoot;
% end

Ts = 1e-4;
t = (0:length(data)-1)*Ts;

t_SPIKES = t_spikes/1000;

% id_spikes = floor(t_spikes*10);
% spikeTrigger = false(1, length(data));
% spikeTrigger(id_spikes) = true;
%  
% h = 1e-4; 
% time_window = 500E-3;
% time_slide = 10E-4;

% slide_steps = floor(time_slide/h);
% length_steps = floor(length(spikeTrigger)/slide_steps);
% windows_steps = floor(floor(time_window/h)/slide_steps);
% spike_hist_tmp = spikeTrigger(1 : (length_steps*slide_steps));
% spike_hist_tmp = reshape(spike_hist_tmp, slide_steps, length_steps);
% spike_hist = sum(spike_hist_tmp);
% kernel = zeros(1, length_steps);
% kernel(1:windows_steps) = 1;
% rate_tmp = conv(spike_hist, kernel);
% rate = rate_tmp(windows_steps : length_steps)/time_window;
% y = 0 : time_slide : length(data);
% y = y(1 : length(rate));


[radX, SpikeTimeX] = computeRadFiring_v02(t_events, t_SPIKES, 1 : length(t_SPIKES));

[FIRSTtime, LASTtime] = sortFirstLast_v01(t_events, t_SPIKES, 1 : length(t_SPIKES));

ifirst111 = [];
ilast111 = [];
for i = 1 : length(FIRSTtime)
    if FIRSTtime(i) > 0
        ifirst111(i) = find(t_SPIKES == FIRSTtime(i));
        if LASTtime(i) > 0
            ilast111(i) = find(t_SPIKES == LASTtime(i));
        else
            ilast111(i) = ifirst111(i);    
        end
    end
end
ifirst111 = nonzeros(ifirst111);
ilast111 = nonzeros(ilast111);

puntiX=10;

for i = 1 : length(ifirst111)
    daXF = t_SPIKES(ifirst111(i))/Ts-puntiX;
    aXF = t_SPIKES(ifirst111(i))/Ts+puntiX;
end

for i = 1 : length(t_SPIKES)
    daXF = t_SPIKES(i)/Ts-puntiX;
    aXF = t_SPIKES(i)/Ts+puntiX;
end

for i = 1 : length(ifirst111)
    daXL = t_SPIKES(ilast111(i))/Ts-puntiX;
    aXL = t_SPIKES(ilast111(i))/Ts+puntiX;
end

FR_tot = length(t_SPIKES)/t(end);

for i = 1 : length(radX(1,:))
    DTcycle(i) = (t_events(i+1)-t_events(i));
    FR_cycle(i) = length(nonzeros(radX(:,i)))/DTcycle(i);
end

Xrad_mean_plus = [];
Xrad_min_plus = [];
Xrad_max_plus = [];
DTcycle_plus = [];
FR_cycle_plus = [];

for w = 1 : length(t_events)
    for i = 1
        try
            t_events_plus = [];
            t_events_plus = [t_events(w) t_events(w+i)];
            
            [radX_plus, SpikeTimeX_plus] = computeRadFiring_v02(t_events_plus, t_SPIKES, 1 : length(t_SPIKES));
            
            [Xrad_mean_plus_i, Xrad_min_plus_i, Xrad_max_plus_i] = ComputeLatency(radX_plus);
            [Xtime_mean_plus_i, Xtime_min_plus_i, Xtime_max_plus_i] = ComputeLatency(SpikeTimeX_plus);
            
            Xrad_mean_plus(w, i) = mean(Xrad_mean_plus_i);
            Xtime_mean_plus(w, i) = mean(Xtime_mean_plus_i);
            
            Xrad_min_plus(w, i) = mean(Xrad_min_plus_i);
            Xtime_min_plus(w, i) = mean(Xtime_min_plus_i);
            
            Xrad_max_plus(w, i) = mean(Xrad_max_plus_i);
            Xtime_max_plus(w, i) = mean(Xtime_max_plus_i);
            
            if i == 1
                for k = 1 : 10
                    try
                        XtimeFirst_min_plus(w, k) = Xtime_min_plus(w, 1)+(t_events(w)-t_events(w-k+1));
                        XtimeAVG_mean_plus(w, k) = Xtime_mean_plus(w, 1)+(t_events(w)-t_events(w-k+1));
                        XtimeLAST_max_plus(w, k) = Xtime_max_plus(w, 1)+(t_events(w)-t_events(w-k+1));
                    end
                end
            end
            
            Xrad_max_plus(w, i) = Xrad_max_plus_i(1);
            Xtime_max_plus(w, i) = mean(Xtime_max_plus_i);
            
            for j = 1 : length(radX_plus(1, :))
                DTcycle_plus(i, j) = (t_events_plus(j+1)-t_events_plus(j));
                FR_cycle_plus(i, j) = length(nonzeros(radX_plus(:,j)))/DTcycle_plus(i,j);
            end
    
        end
    end
end


for w=1:length(t_events)%-1
    for i=1%:10
        try
            t_events_minus = [];
            t_events_minus = unique(-[t_events(w), t_events(w+i)]);
                        
            [radX_minus, SpikeTimeX_minus] = computeRadFiring_v02(t_events_minus, unique(-t_SPIKES), 1 : length(t_SPIKES));
            
            [Xrad_mean_minus_i, Xrad_min_minus_i, Xrad_max_minus_i] = ComputeLatency(radX_minus);
            [Xtime_mean_minus_i, Xtime_min_minus_i, Xtime_max_minus_i] = ComputeLatency(SpikeTimeX_minus);
            
            Xrad_mean_minus(w, i) = mean(Xrad_mean_minus_i);
            Xtime_mean_minus(w, i) = mean(Xtime_mean_minus_i);
            
            Xrad_min_minus(w, i) = mean(Xrad_min_minus_i);
            Xtime_min_minus(w, i) = mean(Xtime_min_minus_i);
            
            Xrad_max_minus(w, i) = Xrad_max_minus_i(1);
            Xtime_max_minus(w, i) = mean(Xtime_max_minus_i);
            
            if i == 1
                for k = 1 : 10
                    try
                        XtimeFirst_min_minus(w, k) = Xtime_min_minus(w,1)+(-t_events(w)+t_events(w+k));
                        XtimeAVG_mean_minus(w, k) = Xtime_mean_minus(w,1)+(-t_events(w)+t_events(w+k));
                        XtimeLAST_max_minus(w, k) = Xtime_max_minus(w,1)+(-t_events(w)+t_events(w+k));
                    end
                end
            end
            for j = 1 : length(radX_minus(1,:))
                DTcycle_minus(i, j) = (t_events_minus(j+1)-t_events_minus(j));
                FR_cycle_minus(i, j) = length(nonzeros(radX_minus(:, j)))/DTcycle_minus(i, j);
            end
        end
    end
end

for i = 1 : 10
    XtimeFirst_min_plus_sorted(:, i) = XtimeFirst_min_plus(:, 11-i); 
    XtimeAVG_mean_plus_sorted(:, i) = XtimeAVG_mean_plus(:, 11-i);
    XtimeLAST_max_plus_sorted(:, i) = XtimeLAST_max_plus(:,  11-i); 
end

FIRSTx10 = [XtimeFirst_min_plus_sorted, XtimeLAST_max_minus];
AVGx10 = [XtimeAVG_mean_plus_sorted, XtimeAVG_mean_minus];
LASTx10 = [XtimeLAST_max_plus_sorted, XtimeFirst_min_minus];

try
    BP_FR = [bpValues(int32(AAAt_footIndex))', FR_cycle', DTcycle'];
catch
    try
        aqws = bpValues(int32(AAAt_footIndex))';  
        swed = FR_cycle(1 : length(aqws))';
        derf = DTcycle(1 : length(aqws))';
        BP_FR = [bpValues(int32(AAAt_footIndex))', swed, derf];
    catch
        try
            aqws1 = bpValues(int32(AAAt_footIndex(1 : end-1)))';  
            swed1 = FR_cycle(1 : length(aqws1))';
            derf1 = DTcycle(1 : length(aqws1))';
            BP_FR = [aqws1, swed1, derf1];   
        catch
            try
                aqws2 = bpValues(int32(AAAt_footIndex(1 : end-2)))';  
                swed2 = FR_cycle(1 : length(aqws2))';
                derf2 = DTcycle(1 : length(aqws2))';
                BP_FR = [aqws2, swed2, derf2];       
            catch
                aqws3 = bpValues(int32(AAAt_footIndex(1 : end-3)))';  
                swed3 = FR_cycle(1 : length(aqws3))';
                derf3 = DTcycle(1 : length(aqws3))';
                BP_FR = [aqws3, swed3, derf3];  
            end
        end
    end
end


FIRSTx10 = [FIRSTx10; zeros(3, length(FIRSTx10(1, :)))];
AVGx10 = [AVGx10; zeros(3, length(AVGx10(1, :)))];
LASTx10 = [LASTx10; zeros(3, length(LASTx10(1, :)))];


RRsenzaSpike = [];

for r = 1 : length(FIRSTx10(:, 1))
    if sum(isnan(FIRSTx10(r, :))) > 0
        FIRSTx10(r, :) = NaN; 
        RRsenzaSpike(r) = 1;
    end
end


for r = 1 : length(AVGx10(:, 1))
    if sum(isnan(AVGx10(r, :))) > 0
       AVGx10(r, :) = NaN; 
    end
end

for r = 1 : length(LASTx10(:, 1))
    if sum(isnan(LASTx10(r, :))) > 0
       LASTx10(r, :) = NaN; 
    end
end

for r = 1 : length(FIRSTx10(:, 1))
    for c = 1 : length(FIRSTx10(1, :))
        if isnan(FIRSTx10(r, c)) > 0
             if isnan(FIRSTx10(r+1, c))==0 & FIRSTx10(r+1, c)>0
           FIRSTx10(r, c) = DTcycle(r)+FIRSTx10(r+1, c); 
             else
                 if isnan(FIRSTx10(r+2, c))==0 & FIRSTx10(r+2, c)>0
                FIRSTx10(r, c) = DTcycle(r)+ DTcycle(r+1)+FIRSTx10(r+2, c);
                 else
                     try
                     if isnan(FIRSTx10(r+3, c))==0 & FIRSTx10(r+3, c)>0
                        FIRSTx10(r, c) =DTcycle(r)+ DTcycle(r+1)+ DTcycle(r+2)+FIRSTx10(r+3, c);
                     else
                          FIRSTx10(r, c) = 0;
                     end
                     catch
                         FIRSTx10(r, c) = 0;
                     end
                 end
             end
        end
    end
end

for r = 1 : length(AVGx10(:, 1))
    for c = 1 : length(AVGx10(1, :))
        if isnan(AVGx10(r, c)) > 0
             if isnan(AVGx10(r+1, c))==0 & AVGx10(r+1, c)>0
           AVGx10(r, c) = DTcycle(r)+AVGx10(r+1,c); 
             else
                 if isnan(AVGx10(r+2, c))==0 & AVGx10(r+2, c)>0
                AVGx10(r, c)=DTcycle(r)+ DTcycle(r+1)+AVGx10(r+2, c);
                 else
                     try
                     if isnan(AVGx10(r+3, c))==0 & AVGx10(r+3, c)>0
                        AVGx10(r, c)=DTcycle(r)+ DTcycle(r+1)+ DTcycle(r+2)+AVGx10(r+3, c);
                     else
                          AVGx10(r, c)=0;
                     end
                     catch
                          AVGx10(r, c)=0;
                     end
                 end
             end
        end
    end
end
        
for r = 1 : length(LASTx10(:, 1))
    for c = 1 : length(LASTx10(1, :))
        if isnan(LASTx10(r, c)) > 0
             if isnan(LASTx10(r+1, c))==0 & LASTx10(r+1, c)>0
           LASTx10(r, c) = DTcycle(r)+LASTx10(r+1, c); 
             else
                 if isnan(LASTx10(r+2, c))==0 & LASTx10(r+2, c)>0
                LASTx10(r, c) = DTcycle(r)+ DTcycle(r+1)+LASTx10(r+2, c);
                 else
                     try
                     if isnan(LASTx10(r+3,c))==0 & LASTx10(r+3,c)>0
                        LASTx10(r, c) = DTcycle(r)+ DTcycle(r+1)+ DTcycle(r+2)+LASTx10(r+3, c);
                     else
                          LASTx10(r, c) = 0;
                     end
                     catch
                        LASTx10(r, c) = 0; 
                     end
                 end
             end
        end
    end
end

for i = 1 : 9
    FIRSTx10(1 : 10-i, i) = 0;
    AVGx10(1 : 10-i, i) = 0;
    LASTx10(1 : 10-i, i) = 0;  
end

i_RRsenzaSpike = find(RRsenzaSpike == 1)';

FIRSTx10_daSaR = FIRSTx10; 
% AVGx10_daSaR = AVGx10; 
% LASTx10_daSaR = LASTx10; 

FIRSTx10_daSaR(i_RRsenzaSpike, :) = []; 
% AVGx10_daSaR(i_RRsenzaSpike, :) = []; 
% LASTx10_daSaR(i_RRsenzaSpike, :) = []; 

for plotDaSaR=1
    for i = 1 : 10
        BP_FR_daSaR = [];   
        BP_FR_daSaR = BP_FR;

        if length(i_RRsenzaSpike) >= 1
            if i+i_RRsenzaSpike(end)-1 <= length(BP_FR_daSaR)
                BP_FR_daSaR(i+i_RRsenzaSpike-1, :) = [];
            end
        end
    
        fineBP = min(length(BP_FR_daSaR), i+length(nonzeros(FIRSTx10_daSaR(:, 10+i)))-1);
        
        A = BP_FR_daSaR(i : fineBP, 1);
        B = nonzeros(FIRSTx10_daSaR(:, 10+i));
%         C=nonzeros(AVGx10_daSaR(:,10+i));
%         D=nonzeros(LASTx10_daSaR(:,10+i));
        
        [~, iirmOA] = rmoutliers(A);
        [~, iirmOB] = rmoutliers(B);
        
        iirmO = zeros(length(A), 1);
        iirmO(iirmOA) = 1;
        iirmO(iirmOB) = 1;
    
        if max(iirmO) > 0
            A1 = A(not(iirmO));
            B1 = B(not(iirmO));
%             C1 = C(not(iirmO));
%             D1 = D(not(iirmO));
        else
            A1 = A;
            B1 = B;
%             C1 = C;
%             D1 = D;
        end
        
        [R, ~] = corrcoef(A1, B1);
        RRR7(10+i, 1)=R(1, 2);

        try            
            A = BP_FR_daSaR(1: length(nonzeros(FIRSTx10_daSaR(:, 11-i))), 1);
            B = -nonzeros(FIRSTx10_daSaR(:, 11-i));
%             C = -nonzeros(AVGx10_daSaR(:, 11-i));
%             D = -nonzeros(LASTx10_daSaR(:, 11-i));
            
            [~, iirmOA] = rmoutliers(A);
            [~, iirmOB] = rmoutliers(B);
            iirmO = zeros(length(A), 1);
            iirmO(iirmOA) = 1;
            iirmO(iirmOB) = 1;
            
            if max(iirmO) > 0
                A1 = A(not(iirmO));
                B1 = B(not(iirmO));
%                 C1 = C(not(iirmO));
%                 D1 = D(not(iirmO));                
            else
                A1 = A;
                B1 = B;
%                 C1 = C;
%                 D1 = D;
            end
            
            [R, ~] = corrcoef(A1, B1);
            RRR7(11-i, 1) = R(1, 2);
        end
    end
    
    RRR8=[];

   
    for i=1:10
        try
%             [R_fr1,PValue,H] = corrplot([BP_FR_daSaR(i+1:end,1),BP_FR_daSaR(1:end-i,2)],'varNames', VariableNames,'type','Pearson','testR','on','alpha',0.05)
            [R_fr1,~] = corr([BP_FR_daSaR(i+1:end,1),BP_FR_daSaR(1:end-i,2)],'type','Pearson','rows','pairwise','tail','both');           
            RRR8(10+i,1)=R_fr1(1,2);
        end
    
        try
%             [R_fr2,PValue,H] = corrplot([BP_FR_daSaR(1:end-i,1),BP_FR_daSaR(i+1:end,2)],'varNames', VariableNames,'type','Pearson','testR','on','alpha',0.05)
            [R_fr2,~] = corr([BP_FR_daSaR(1:end-i,1),BP_FR_daSaR(i+1:end,2)],'type','Pearson','rows','pairwise','tail','both');
            RRR8(11-i,1)=R_fr2(1,2);
    
        end
    end

    try
        for i = 1 : 10 
            A = BP_FR_daSaR(i : i+length(nonzeros(FIRSTx10_daSaR(:, 10+i)))-1, 3);
            B = nonzeros(FIRSTx10_daSaR(:, 10+i));
%             C = nonzeros(AVGx10_daSaR(:, 10+i));
%             D = nonzeros(LASTx10_daSaR(:, 10+i));
            
            [~, iirmOA] = rmoutliers(A);
            [~, iirmOB] = rmoutliers(B);
            iirmO = zeros(length(A), 1);
            iirmO(iirmOA) = 1;
            iirmO(iirmOB) = 1;
            
            if max(iirmO) > 0
                A1 = A(not(iirmO));
                B1 = B(not(iirmO));
%                 C1 = C(not(iirmO));
%                 D1 = D(not(iirmO));
            else
                A1 = A;
                B1 = B;
%                 C1 = C;
%                 D1 = D;
            end
            
            [R, ~] = corrcoef(A1, B1);
            RRR7(10+i, 2) = R(1, 2);

            try
                A = BP_FR_daSaR(1 : length(nonzeros(FIRSTx10_daSaR(:, 11-i))), 3);
                B = -nonzeros(FIRSTx10_daSaR(:, 11-i));
%                 C = -nonzeros(AVGx10_daSaR(:, 11-i));
%                 D = -nonzeros(LASTx10_daSaR(:, 11-i));
                
                [~, iirmOA] = rmoutliers(A);
                [~, iirmOB] = rmoutliers(B);
                iirmO = zeros(length(A), 1);
                iirmO(iirmOA) = 1;
                iirmO(iirmOB) = 1;
                
                if max(iirmO) > 0
                    A1 = A(not(iirmO));
                    B1 = B(not(iirmO));
%                     C1 = C(not(iirmO));
%                     D1 = D(not(iirmO));
                else
                    A1 = A;
                    B1 = B;
%                     C1 = C;
%                     D1 = D;
                end
                
                [R,~] = corrcoef(A1, B1);
                
                RRR7(11-i, 2) = R(1, 2);                  
            end
        end            
    end

    for i = 1 : 10
        try
%             VariableNames ={'DT','FRate'};
%             [R_fr3,PValue,H] = corrplot([BP_FR_daSaR(i+1:end,3),BP_FR_daSaR(1:end-i,2)],'varNames', VariableNames,'type','Pearson','testR','on','alpha',0.05)
            [R_fr3,PValue] = corr([BP_FR_daSaR(i+1:end,3),BP_FR_daSaR(1:end-i,2)],'type','Pearson','rows','pairwise','tail','both');
            RRR8(10+i, 2) = R_fr3(1, 2);
        end
        
        try
%             [R_fr4,PValue,H] = corrplot([BP_FR_daSaR(1:end-i,3),BP_FR_daSaR(i+1:end,2)],'varNames', VariableNames,'type','Pearson','testR','on','alpha',0.05)
            [R_fr4,PValue] = corr([BP_FR_daSaR(1:end-i,3),BP_FR_daSaR(i+1:end,2)],'type','Pearson','rows','pairwise','tail','both');
            RRR8(11-i, 2) = R_fr4(1, 2);
        end
    end        
end

[Xtime_mean, Xtime_min, Xtime_max] = ComputeLatency(SpikeTimeX);

for i = 2 : length(Xtime_min)
    try
        TTT_000(i) = t_events(i)-t_events(i-1);
        FIFIFI_min(i) = Xtime_min(i)/TTT_000(i);
        FIFIFI_mean(i) = Xtime_mean(i)/TTT_000(i);
        FIFIFI_max(i) = Xtime_max(i)/TTT_000(i);
        for j = 1 : 10
            TTT_post(j) = t_events(i+j)-t_events(i);
            DELTA_FI(j, i) = (j*TTT_000(i)-TTT_post(j))/TTT_000(i);
            GIGIGI_FI_min(j, i) = FIFIFI_min(i) + (j*TTT_000(i)-TTT_post(j))/TTT_000(i);
            GIGIGI_FI_mean(j, i) = FIFIFI_mean(i) + (j*TTT_000(i)-TTT_post(j))/TTT_000(i);
            GIGIGI_FI_max(j, i) = FIFIFI_max(i) + (j*TTT_000(i)-TTT_post(j))/TTT_000(i);
        end
    end
end

for i = 1 : 10
    try
        A = BP_FR_daSaR(i : i+length(nonzeros(GIGIGI_FI_min(i, :)))-1, 3);
    catch
        A = BP_FR_daSaR(i : i+length(nonzeros(GIGIGI_FI_min(i, :)))-10, 3);
    end
    
    B = nonzeros(GIGIGI_FI_min(i, :));
%     C = nonzeros(GIGIGI_FI_mean(i, :));
%     D = nonzeros(GIGIGI_FI_max(i, :));
    
    [~, iirmOA] = rmoutliers(A);
    [~, iirmOB] = rmoutliers(B);
    iirmO = zeros(length(A), 1);
    iirmO(iirmOA) = 1;
    iirmO(iirmOB) = 1;
    
    if max(iirmO) > 0
        try
            A1 = A(not(iirmO));
            B1 = B(not(iirmO));
%             C1 = C(not(iirmO));
%             D1 = D(not(iirmO));
        catch
            A1 = A;
            B1 = B;
%             C1 = C;
%             D1 = D;   
        end
    else
        A1 = A;
        B1 = B;
%         C1 = C;
%         D1 = D;
    end
    
    lllab = min(length(A1), length(B1));
    A2 = A1(1 : lllab);
    B2 = B1(1 : lllab);
    A1 = [];
    B1 = [];
    A1 = A2;
    B1 = B2;

    [R, ~] = corrcoef(A1, B1);
    
    RRR7(10+i, 4) = R(1, 2);  
end

for i = 1 : 10
    try
        A = BP_FR_daSaR(i : i+length(nonzeros(GIGIGI_FI_min(i, :)))-1, 1);
    catch
        A = BP_FR_daSaR(i : i+length(nonzeros(GIGIGI_FI_min(i, :)))-10, 1);
    end
    
    B = nonzeros(GIGIGI_FI_min(i, :));
%     C = nonzeros(GIGIGI_FI_mean(i, :));
%     D = nonzeros(GIGIGI_FI_max(i, :));
    
    [~, iirmOA] = rmoutliers(A);
    [~, iirmOB] = rmoutliers(B);
    iirmO = zeros(length(A), 1);
    iirmO(iirmOA) = 1;
    iirmO(iirmOB) = 1;
    
    if max(iirmO) > 0
        try
            A1 = A(not(iirmO));
            B1 = B(not(iirmO));
%             C1 = C(not(iirmO));
%             D1 = D(not(iirmO));
        catch
            A1 = A;
            B1 = B;
%             C1 = C;
%             D1 = D;   
        end
    else
        A1 = A;
        B1 = B;
%         C1 = C;
%         D1 = D;
    end        
    
    lllab = min(length(A1), length(B1));
    A2 = A1(1 : lllab);
    B2 = B1(1 : lllab);
    A1 = [];
    B1 = [];
    A1 = A2;
    B1 = B2;
    [R, ~] = corrcoef(A1, B1); 
    RRR7(10+i, 3) = R(1, 2);
end 
% RRR7

% RRR8

end

