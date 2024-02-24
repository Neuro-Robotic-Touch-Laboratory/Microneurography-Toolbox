function [INDEXES] = getINDEXES(t, tevents, window4mean)
    INDEXES = struct();
    INDEXES.t_events = tevents;
    INDEXES.dt_instantaneous = diff(tevents);
    INDEXES.rate_instantaneous = 60./diff(tevents);
    
    if length(tevents) > 1
        
      %INTERVALLO
        INDEXES.t = t(t > tevents(2)); 
        windowIntervalMean_fun = @(t) mean(diff(tevents(t-window4mean<tevents & tevents<t)));
        INDEXES.dt_mean = arrayfun(windowIntervalMean_fun, INDEXES.t);
        
      %RATE
        INDEXES.rate_mean = 60./INDEXES.dt_mean;
        
        rate_mean_corrected=INDEXES.rate_mean;
        rate_mean_corrected(isnan(INDEXES.rate_mean))=[];
        INDEXES.rate_meanTot = mean(rate_mean_corrected);
    
      %VARIABILITà
        INDEXES.dtDifferences = diff(INDEXES.dt_instantaneous);
        INDEXES.rateVariability = sqrt(mean(INDEXES.dtDifferences.^2));  %RMSSD = “root mean square of successive differences” An Overview of Heart Rate Variability Metrics and Norms Fred Shaffer and J. P. Ginsberg
        INDEXES.rateVariability2 = mean(abs(INDEXES.dtDifferences.^2) > 0.05); %pNN50
    
    else
        INDEXES.t = []; 
        INDEXES.dt_mean = [];
        INDEXES.rate_mean = [];
        INDEXES.rate_meanTot =[];
        INDEXES.dtDifferences = [];
        INDEXES.rateVariability = []; 
        INDEXES.rateVariability2 = [];
    
    end
    
end