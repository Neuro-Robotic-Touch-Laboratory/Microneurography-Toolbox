function burst_res = burst_analysis(app)
%BURST_ANALYSIS Summary of this function goes here
%   Detailed explanation goes here
cla(app.ax_msna_raw), cla(app.ax_burst_interval_his), cla(app.ax_burst_duration_his),...
cla(app.ax_burst_amplitude_his),cla(app.ax_burst_integral_his), cla(app.ax_burst_res),...
cla(app.ax_msna_int)
Soglia = app.edt_threshold.Value;
Finestra = app.edt_window.Value;
t_i = app.settings.interval(1,1);
t_f = app.settings.interval(1,2);
t_RMVi = [];
t_RMVf = [];

MSNA = app.data(app.settings.channel_idx.msna).data;
sel = app.data(app.settings.channel_idx.msna).data(app.settings.interval(1,1)/app.data(app.settings.channel_idx.msna).ts(1):app.settings.interval(1,2)/app.data(app.settings.channel_idx.msna).ts(1)); 
sel(1) = []; 
units_msna = app.data(app.settings.channel_idx.msna).unit;
t_msna = (1:length(sel))*app.data(app.settings.channel_idx.msna).ts(1);
Ts_msna = app.data(app.settings.channel_idx.msna).ts(1);
CalcolaBurst=1;
app.lbl_working.Text = 'Working 10%'; drawnow
if CalcolaBurst==1
    sel=bandpass(sel,[300 5000],1e4);


    tmp = sqrt(movmean(sel.^2,2001));
    rmsa1 = tmp';
    

    rmsa2=rmsa1;

    if app.chkbx_detrent_rms.Value
        x=detrend(rmsa2,2);
    else
        x=rmsa2;
    end
    t_msna = (1:length(x))*Ts_msna;
    
    WindowSize = 10.0; % [s]
    adaptiveWindowSize = round(WindowSize/Ts_msna);
    sampleLength = length(x);
    Ksoglia= Soglia;
    PeakWidth = Finestra; % [seconds]
    PeakWidth = round(PeakWidth/Ts_msna);
    app.lbl_working.Text = 'Working 20%';drawnow

    basmean = zeros(1, sampleLength);
    basStd = zeros(1, sampleLength);


    HThreshold=[];
    tbursts=[];

    
    debug=0;


app.lbl_working.Text = 'Working 30%'; drawnow
%     if debug==0
%         x_windowValues = x(1:end);
%         basmed = median(x_windowValues);
%         basStd=std(x_windowValues);
%         HThreshold=  basmed+ Ksoglia*basStd;  
%         HThresholdOutlier=  basmed+ app.edt_outlier_threshold.Value*Ksoglia*basStd;    %% outlier threshold     
%     end

    if debug==0
        x_windowValues = x(1:end);
        basmed = median(x_windowValues);
        basStd=std(x_windowValues);
    
        if app.chkbx_manual.Value %gggggggggg
            HThreshold=app.edt_threshold_abs.Value;
        else
            HThreshold=  basmed+ Ksoglia*basStd; 
        end
        %HThreshold=  basmed+ Ksoglia*basStd;  
        
        if app.chkbx_th_on_silent.Value
            tic%%
            x2 = x;
            bidx = find (x2 >= HThreshold);
            strt = [bidx(1),1];
            for i=2:length(bidx)
                if bidx(i)~=bidx(i-1)+1
                    if (bidx(i-1)-strt(1)+1) < PeakWidth
                        bidx(strt(2):i-1) = nan; 
                    end
                    strt= [bidx(i),i];
                end
            end
            
            bidx(isnan(bidx)) = [];
            x2(bidx) = [];
            sil_med = median(x2);
            sil_Std = std(x2);
            HThresholdOutlier = sil_med +app.edt_outlier_threshold.Value *Ksoglia *sil_Std;
            clear x2
            toc%%
        else
            HThresholdOutlier = basmed +app.edt_outlier_threshold.Value *Ksoglia *basStd;    %% outlier threshold     
        end
    end
   
 %% fix outlier rejection 
    all_areas = ((x - HThreshold) > 0 & (x - HThresholdOutlier) < 0);     
%     tic %%
%     peak_areas2 = imopen( all_areas > 0, strel('line', PeakWidth, 0)); % what the hell does it do ?
%     toc%%
%%  
%     tic
    peak_areas = all_areas;

    idx = find (peak_areas);
    strt = [idx(1),1];
    for i=2:length(idx)
        if idx(i)~=idx(i-1)+1
            if (idx(i-1)-strt(1)) < PeakWidth
                peak_areas(idx(strt(2)):idx(i-1)) = false;
                %idx(strt(2):i-1) = false; 
            end
            strt= [idx(i),i];
        end
    end
    
%     isequal (peak_areas2,peak_areas)
%     toc
    %%
    tmp_idx = find(peak_areas == 0,1);          % remove bursts that are not enirely in the intervall
    peak_areas(1:tmp_idx) = 0;                  %   
    tmp_idx = find(peak_areas == 0,1,'last');   %
    peak_areas(tmp_idx:end) = 0;                %
    areaStart_idx = find((peak_areas(1:end-1)~=1 & peak_areas(2:end)==1)) + 1;
    areaEnd_idx = find((peak_areas(1:end-1)==1 & peak_areas(2:end)~=1));

    if length(areaEnd_idx)>length(areaStart_idx)
        areaEnd_idx(1)=[];
    end
    app.lbl_working.Text = 'Working 40%'; drawnow
    if length(areaStart_idx)>length(areaEnd_idx)
        areaStart_idx(end)=[];
    end
    
    bursttMask = diff(peak_areas(areaStart_idx)) == 0;
    burst_idx = round((areaStart_idx(bursttMask) + areaEnd_idx(bursttMask))/2);

    burstTrigger = false(1, sampleLength);
    burstTrigger(burst_idx) = true;
    tbursts(1:length(burst_idx)) = t_msna(burst_idx)';
    n_Peaks=sum(bursttMask);

    try
        t_Peak_mean= mean(areaEnd_idx-areaStart_idx)*Ts_msna;
        t_Peaks=sum(areaEnd_idx-areaStart_idx)*Ts_msna;
        NOT_t_Peaks=t_msna(end)-sum(areaEnd_idx-areaStart_idx)*Ts_msna;
    end



    IpeakAreas=find(peak_areas==1);
    HThreshold=mean(HThreshold);
    x_shift=x-HThreshold;
    app.lbl_working.Text = 'Working 50%'; drawnow
    IpeakAreas2=[areaStart_idx';areaEnd_idx'];
    
    x_shift_corretto=x_shift;
    x_shift_corretto(1)=0;
    x_shift_corretto(end)=0;
    x_shift_corretto(unique(IpeakAreas2))=0;
    

 
    IpeakAreasLoc=[areaStart_idx',areaEnd_idx'];
    
    tbursts=[tbursts mean(IpeakAreasLoc(end,:)*Ts_msna)];


    
    AA_t_I=[tbursts' IpeakAreasLoc];
%% remove intervall ?
    index_remove=[];
    try
        for j=1:length(IpeakAreasLoc)
            if (IpeakAreasLoc(j,1)*Ts_msna >  t_RMVi-t_i & IpeakAreasLoc(j,1)*Ts_msna < t_RMVf-t_i)
    
                index_remove=[index_remove j];
            end
    
            for i=1:length(cursor_remove_burst)
                if (cursor_remove_burst(i).Position(1)/Ts_msna >  IpeakAreasLoc(j,1) ...
                        & cursor_remove_burst(i).Position(1)/Ts_msna < IpeakAreasLoc(j,2)) ...
    
                    index_remove=[index_remove j];
    
                end
    
            end
    
        end
    end
app.lbl_working.Text = 'Working 60%'; drawnow
    index_remove=unique(index_remove);
    AA_t_I(index_remove,:)=[];
    DD_t_I= sortrows(AA_t_I,1);
    I_Z_DD=find(DD_t_I(:,2)<=0);
    DD_t_I(I_Z_DD,:)=[];
    %%

    tbursts=DD_t_I(:,1)';
    IpeakAreasLoc=DD_t_I(:,2:3);
    assignin('base','DD_t_I',DD_t_I) 


    %ECG_RESULTS.t_events(find(ECG_RESULTS.t_events> t_RMVi-t_i & ECG_RESULTS.t_events< t_RMVf-t_i))=[]; %thas for hb removal ???
%% plotting top graph (msna)



    burst_res = struct('ts',[],'x',[], 'burst_loc',[],'burst_int',[],...
                       'burst_amp',[],'burst_dur',[],'htresh',[],'xshift',[],...
                       't_burst',[],'dt_burst',[]);

    for ttt=1
%%  plotting  2nd graph (msna rms , rpeaks...) 
        
        burst_res.ts = [t_msna(1), mean(diff(t_msna)), t_msna(end)];
        burst_res.x = x;

        app.lbl_working.Text = 'Working 70%'; drawnow
        for j=1:length(IpeakAreasLoc) 
       
            AUC_IpeakAreasLoc(j)=trapz(t_msna(IpeakAreasLoc(j,1):IpeakAreasLoc(j,2)),x_shift_corretto(IpeakAreasLoc(j,1):IpeakAreasLoc(j,2)));

            AMPL_IpeakAreasLoc(j)=max(x(IpeakAreasLoc(j,1):IpeakAreasLoc(j,2)-50)); 

        end
   
        [ffffrrr iirmOA]=rmoutliers(AUC_IpeakAreasLoc,'percentiles',[0 100]);
        iirmO=zeros(length(AUC_IpeakAreasLoc),1);
        iirmO(iirmOA)=1;
        
        
        tAUC1=[];
        AUC1=[];
    
        if max(iirmO)>0
            tAUC1=tbursts(not(iirmO));
            AUC1=1000*AUC_IpeakAreasLoc(not(iirmO));
        else
            tAUC1=tbursts;
            AUC1=1000*AUC_IpeakAreasLoc;
        end
        app.lbl_working.Text = 'Working 80%'; drawnow
        tbursts=tAUC1;
        AUC_IpeakAreasLoc = AUC1;

        burst_res.t_burst = tbursts;
        burst_res.dt_burst = diff(tbursts);
        burst_res.burst_int = AUC_IpeakAreasLoc;
        burst_res.burst_amp = AMPL_IpeakAreasLoc;
        burst_res.htresh = HThreshold;
        burst_res.xshift = x_shift_corretto;
        burst_res.use_burst = ones(size(tbursts,2),2);


%% save msna results

%         debugUUU=1;
%         if debugUUU==1
%             for i=3;%[2 5]
%                 MSNA_RESULTS=[];
%                 window4mean=i;
%                 
%                 MSNA_RESULTS=getINDEXES(t_msna, tbursts, window4mean); %% might be improved TOOO SLOW !!
%                 
%             end
%         end
%         app.lbl_working.Text = 'Working 90%'; drawnow
%         try
%         MSNA_RESULTS.nPeaks= n_Peaks;
%         MSNA_RESULTS.tPeaks= t_Peaks;
%         MSNA_RESULTS.Not_tPeaks= NOT_t_Peaks;
%         MSNA_RESULTS.t_Peak_mean= t_Peak_mean;
%         MSNA_RESULTS.AVGintegral=mean(100*AUC_IpeakAreasLoc);
%         end
    
        %imin3=diff(tbursts)<0.3;
        %tbursts_cl=tbursts;

    %%lower left graph (burst integral)

        %CC.burst_integral=AUC_IpeakAreasLoc;

    %% lower middle left graph (burst amplitude)

        %FF.burst_amplitude=AMPL_IpeakAreasLoc;


%%lower middle right graph (burst duration)
    
        DD.burst_duration=(IpeakAreasLoc(:,2)-IpeakAreasLoc(:,1))*Ts_msna;
        burst_res.burst_dur = DD.burst_duration;
        burst_res.burst_loc = IpeakAreasLoc;


        %EE.inter_burst_interval=MSNA_RESULTS.dt_instantaneous;

    end

end

