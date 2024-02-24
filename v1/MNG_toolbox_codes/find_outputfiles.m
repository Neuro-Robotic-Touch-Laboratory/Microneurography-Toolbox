function [figs, path] = find_outputfiles(app)

path = uigetdir;

listing = dir(path);

f_names = {[]};
f_idx = 1;
for i = 1: length(listing)
    if ~listing(i).isdir 
        if ~contains(listing(i).name,'xls' )
            tmp_idx = [];
            tmp_idx = strfind(listing(i).name,'_INT');
            if ~isempty(tmp_idx)
                f_names{f_idx,1} = listing(i).name(1:tmp_idx(1)-1);
                listing(i).name(1:tmp_idx+4)= [];
                tmp_idx = find(listing(i).name== '_',1);
                f_names{f_idx,2} = listing(i).name(1:tmp_idx-1);
                f_names{f_idx,3} = listing(i).name(tmp_idx+1:end);
                f_idx = f_idx+1;
            end
        end
    end
end

% figs = struct('name', [], 'tmp', {[]},'ints', struct('name', [], 'tmp', [],...
%                              'anas',struct('burst', struct('done', false, 'amp', struct('format', [false,false,false]),...
%                                                                           'dur', struct('format', [false,false,false]),...
%                                                                           'int', struct('format', [false,false,false]),...
%                                                                           'rr',  struct('format', [false,false,false])),...
%                                            'spike', struct('done', false, 'subint', struct('name',[], 'unit',struct('name',[],'format', [false,false,false]))),... 
%                                            'spectral', struct('done', false,'subint', struct('name', [], 'signal', struct('name', [], 'format', [false,false,false]))),...
%                                            'wavelet', struct('done', false,'subint', struct('name', [], 'signal', struct('name', [], 'format', [false,false,false]))),...
%                                            'entropy', struct(),...
%                                            'correlation',struct('done', false, 'subint', struct('name',[], 'signal', struct('name', [], 'format', [false,false,false]))),...
%                                            'stepwise', struct('done', false, 'subint', struct('name', [],'format', [false,false,false])),...
%                                            'annotate', struct('done', false, 'subint', struct('name', [], 'format',[false,false,false])),...
%                                            'transduction', struct('done', false, 'subint', struct('name', [], 'format',[false,false,false])))));

figs = struct('name', [], 'tmp', {[]},'ints', struct('name', [], 'tmp', [],...
                             'anas',[]));
if ~isempty(f_names{1})
    tmp_names = unique(f_names(:,1));
    for i =1:length(tmp_names)
        figs(i).name = tmp_names{i};
        tmp_idx = strcmp(f_names(:,1),figs(i).name);
        figs(i).tmp = f_names(tmp_idx,2:3);
        tmp_ints = unique(figs(i).tmp(:,1));
        for j = 1:length(tmp_ints)
            figs(i).ints(j).name = tmp_ints{j};
            figs(i).ints(j).anas = struct('burst', struct('done', false, 'amp', struct('format', [false,false,false]),...
                                                                              'dur', struct('format', [false,false,false]),...
                                                                              'int', struct('format', [false,false,false]),...
                                                                              'rr',  struct('format', [false,false,false])),...
                                               'spike', struct('done', false, 'subint', struct('name',[],'format', [false,false,false], 'unit',struct('name',[],'format', [false,false,false]))),... 
                                               'spectral', struct('done', false,'subint', struct('name', [], 'signal', struct('name', [], 'format', [false,false,false]))),...
                                               'wavelet', struct('done', false,'subint', struct('name', [], 'signal', struct('name', [], 'format', [false,false,false]))),...
                                               'entropy', struct('done',false, 'subint',struct('name',[],'signal',struct('name',[],'format',[false,false,false]))),...
                                               'correlation',struct('done', false, 'subint', struct('name',[], 'signal', struct('name', [], 'format', [false,false,false]))),...
                                               'stepwise', struct('done', false, 'subint', struct('name', [],'format', [false,false,false])),...
                                               'annotate', struct('done', false, 'subint', struct('name', [], 'format',[false,false,false])),...
                                               'transduction', struct('done', false, 'subint', struct('name', [], 'format',[false,false,false])));
            tmp_idx = strcmp(figs(i).tmp,figs(i).ints(j).name);
            figs(i).ints(j).tmp = figs(i).tmp(tmp_idx,2);
            for k = 1 : length(figs(i).ints(j).tmp)
                if strfind(figs(i).ints(j).tmp{k}, 'burst_analysis_')
                    figs(i).ints(j).anas.burst.done = true; 
                    figs(i).ints(j).tmp{k}(1:15) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    
                    switch tmp_sint
                        case 'burst_amplitude'
                            figs(i).ints(j).anas.burst.amp.format(m) = true;
                        case 'burst_duration'
                            figs(i).ints(j).anas.burst.dur.format(m) = true;
                        case 'burst_integral'
                            figs(i).ints(j).anas.burst.int.format(m) = true;
                        case 'rr-interval'
                            figs(i).ints(j).anas.burst.rr.format(m) = true;
                    end
                end
                
                if strfind(figs(i).ints(j).tmp{k}, 'annotate_')
                    figs(i).ints(j).anas.annotate.done = true;
                    figs(i).ints(j).tmp{k}(1:9) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
                    if isempty(figs(i).ints(j).anas.annotate.subint(1).name)
                        figs(i).ints(j).anas.annotate.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.annotate.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.annotate.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.annotate.subint);
                            figs(i).ints(j).anas.annotate.subint(l).format = [false,false,false];
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.annotate.subint(l).format(m) = true;
                end
    
                if strfind(figs(i).ints(j).tmp{k}, 'correlation_')
                    figs(i).ints(j).anas.correlation.done = true;
                    figs(i).ints(j).tmp{k}(1:12) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    sig_idx = strfind(figs(i).ints(j).tmp{k}, '_SIG_');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:sig_idx-1);
                    tmp_sigs = figs(i).ints(j).tmp{k}(sig_idx+5:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
                    if isempty(figs(i).ints(j).anas.correlation.subint(1).name)
                        figs(i).ints(j).anas.correlation.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.correlation.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.correlation.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.correlation.subint);
                            figs(i).ints(j).anas.correlation.subint(l).signal = struct('name', [], 'format', [false,false,false]);
                        end
                    end
                    if isempty(figs(i).ints(j).anas.correlation.subint(l).signal(1).name)
                        figs(i).ints(j).anas.correlation.subint(l).signal(1).name = string(tmp_sigs);
                        n = 1;
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.correlation.subint(l).signal.name);
                        n =  find(strcmp(tmp_nms, tmp_sigs));
                        if isempty(n)
                            figs(i).ints(j).anas.correlation.subint(l).signal(end+1).name = string(tmp_sigs);
                            n = length(figs(i).ints(j).anas.correlation.subint(l).signal);
                            figs(i).ints(j).anas.correlation.subint(l).signal(n).format = [false,false,false];
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.correlation.subint(l).signal(n).format(m) = true;
                end
    
                if strfind(figs(i).ints(j).tmp{k}, 'spectral_analysis_')
                    figs(i).ints(j).anas.spectral.done = true;
                    figs(i).ints(j).tmp{k}(1:18) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    sig_idx = strfind(figs(i).ints(j).tmp{k}, '_SIG_');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:sig_idx-1);
                    tmp_sigs = figs(i).ints(j).tmp{k}(sig_idx+5:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
    
                    if isempty(figs(i).ints(j).anas.spectral.subint(1).name)
                        figs(i).ints(j).anas.spectral.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.spectral.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.spectral.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.spectral.subint);
                            figs(i).ints(j).anas.spectral.subint(l).signal = struct('name', [], 'format', [false,false,false]);
                        end
                    end
                    if isempty(figs(i).ints(j).anas.spectral.subint(l).signal(1).name)
                        figs(i).ints(j).anas.spectral.subint(l).signal(1).name = string(tmp_sigs);
                        n = 1;
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.spectral.subint(l).signal(1).name);
                        n =  find(strcmp(tmp_nms, tmp_sigs));
                        if isempty(n)
                            figs(i).ints(j).anas.spectral.subint(l).signal(end+1).name = string(tmp_sint);
                            n = length(figs(i).ints(j).anas.spectral.subint(l).signal);
                            figs(i).ints(j).anas.spectral.subint(l).signal(n) = struct('name', [], 'format', [false,false,false]);
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.spectral.subint(l).signal(n).format(m) = true;
                end
    
                if strfind(figs(i).ints(j).tmp{k}, 'spike_analysis_')
                    figs(i).ints(j).anas.spike.done = true;
                    figs(i).ints(j).tmp{k}(1:15) = [];
                    if strfind(figs(i).ints(j).tmp{k}, 'overview_')
                        figs(i).ints(j).tmp{k}(1:9) = [];
                        tmp_unit = nan;
                    else
                        figs(i).ints(j).tmp{k}(1:5) = [];
                        us_idx = find( figs(i).ints(j).tmp{k} == '_',1);
                        tmp_unit = figs(i).ints(j).tmp{k}(1:us_idx-1); %int16(str2num(figs(i).ints(j).tmp{k}(1:us_idx-1)));
                        figs(i).ints(j).tmp{k}(1:us_idx) = [];
                    end
    
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:pnt_idx -1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx +1 : end);
    
                    if isempty(figs(i).ints(j).anas.spike.subint(1).name)
                        figs(i).ints(j).anas.spike.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.spike.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.spike.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.spike.subint);
                            figs(i).ints(j).anas.spike.subint(l).format = [false,false,false];
                            figs(i).ints(j).anas.spike.subint(l).unit = struct('name', [], 'format', [false,false,false]);
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    if ~isnan(tmp_unit)
    
                        if isempty(figs(i).ints(j).anas.spike.subint(l).unit(1).name)
                            figs(i).ints(j).anas.spike.subint(l).unit(1).name = tmp_unit;
                            n = 1;
                        else
%                             tmp_nms = vertcat(figs(i).ints(j).anas.spike.subint(l).unit.name);
%                             n =  find(strcmp(tmp_nms, tmp_unit));
                            n = [];
                            for o = 1 :length(figs(i).ints(j).anas.spike.subint(l).unit)
                                if strcmp(figs(i).ints(j).anas.spike.subint(l).unit(o).name,tmp_unit)
                                    n = o;
                                end
                            end
                            if isempty(n)
                                figs(i).ints(j).anas.spike.subint(l).unit(end+1).name = tmp_unit;
                                n = length(figs(i).ints(j).anas.spike.subint(l).unit);
                                figs(i).ints(j).anas.spike.subint(l).unit(n).format = [false,false,false];
                            end
                        end
                        
                        figs(i).ints(j).anas.spike.subint(l).unit(n).format(m) = true;
                    else
                        figs(i).ints(j).anas.spike.subint(l).format(m) = true;
                    end
                end
                
    
                if strfind(figs(i).ints(j).tmp{k}, 'stepwise_regression_')
                    figs(i).ints(j).anas.stepwise.done = true;
                    figs(i).ints(j).tmp{k}(1:20) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
                    if isempty(figs(i).ints(j).anas.stepwise.subint(1).name)
                        figs(i).ints(j).anas.stepwise.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.stepwise.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.stepwise.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.stepwise.subint);
                            figs(i).ints(j).anas.stepwise.subint(l).format = [false,false,false];
                        end
                    end
                    
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.stepwise.subint(l).format(m) = true;
                end
    
                if strfind(figs(i).ints(j).tmp{k}, 'transduction_')
                    figs(i).ints(j).anas.transduction.done = true;
                    figs(i).ints(j).tmp{k}(1:13) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
                    if isempty(figs(i).ints(j).anas.transduction.subint(1).name)
                        figs(i).ints(j).anas.transduction.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.transduction.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.transduction.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.transduction.subint);
                            figs(i).ints(j).anas.transduction.subint(l).format = [false,false,false];
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.transduction.subint(l).format(m) = true;
                end
    
                if strfind(figs(i).ints(j).tmp{k}, 'wavelet_')
                    figs(i).ints(j).anas.wavelet.done = true;
                    figs(i).ints(j).tmp{k}(1:17) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    sig_idx = strfind(figs(i).ints(j).tmp{k}, '_SIG_');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:sig_idx-1);
                    tmp_sigs = figs(i).ints(j).tmp{k}(sig_idx+5:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
    
                    if isempty(figs(i).ints(j).anas.wavelet.subint(1).name)
                        figs(i).ints(j).anas.wavelet.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.wavelet.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.wavelet.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.wavelet.subint);
                            figs(i).ints(j).anas.wavelet.subint(l).signal = struct('name', [], 'format', [false,false,false]);
                        end
                    end
                    if isempty(figs(i).ints(j).anas.wavelet.subint(l).signal(1).name)
                        figs(i).ints(j).anas.wavelet.subint(l).signal(1).name = string(tmp_sigs);
                        n = 1;
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.wavelet.subint(l).signal.name);
                        n =  find(strcmp(tmp_nms, tmp_sigs));
                        if isempty(n)
                            figs(i).ints(j).anas.wavelet.subint(l).signal(end+1).name = string(tmp_sigs);
                            n = length(figs(i).ints(j).anas.wavelet.subint(l).signal);
                            figs(i).ints(j).anas.wavelet.subint(l).signal(n).format = [false,false,false];
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.wavelet.subint(l).signal(n).format(m) = true;
                end

                if strfind(figs(i).ints(j).tmp{k}, 'entropy_analysis_')
                    figs(i).ints(j).anas.entropy.done = true;
                    figs(i).ints(j).tmp{k}(1:17) = [];
                    pnt_idx = find(figs(i).ints(j).tmp{k} == '.',1,'last');
                    sig_idx = strfind(figs(i).ints(j).tmp{k}, '_SIG_');
                    tmp_sint = figs(i).ints(j).tmp{k}(1:sig_idx-1);
                    tmp_sigs = figs(i).ints(j).tmp{k}(sig_idx+5:pnt_idx-1);
                    tmp_form = figs(i).ints(j).tmp{k}(pnt_idx+1:end);
    
                    if isempty(figs(i).ints(j).anas.entropy.subint(1).name)
                        figs(i).ints(j).anas.entropy.subint(1).name = string(tmp_sint);
                        l = 1;                    
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.entropy.subint.name);
                        l =  find(strcmp(tmp_nms, tmp_sint));
                        if isempty(l)
                            figs(i).ints(j).anas.entropy.subint(end+1).name = string(tmp_sint);
                            l = length(figs(i).ints(j).anas.entropy.subint);
                            figs(i).ints(j).anas.entropy.subint(l).signal = struct('name', [], 'format', [false,false,false]);
                        end
                    end
                    if isempty(figs(i).ints(j).anas.entropy.subint(l).signal(1).name)
                        figs(i).ints(j).anas.entropy.subint(l).signal(1).name = string(tmp_sigs);
                        n = 1;
                    else
                        tmp_nms = vertcat(figs(i).ints(j).anas.entropy.subint(l).signal.name);
                        n =  find(strcmp(tmp_nms, tmp_sigs));
                        if isempty(n)
                            figs(i).ints(j).anas.entropy.subint(l).signal(end+1).name = string(tmp_sigs);
                            n = length(figs(i).ints(j).anas.entropy.subint(l).signal);
                            figs(i).ints(j).anas.entropy.subint(l).signal(n).format = [false,false,false];
                        end
                    end
                    switch tmp_form
                        case 'fig'
                            m = 1;
                        case 'jpeg'
                            m = 2;
                        case 'epsc'
                            m = 3;
                    end
                    figs(i).ints(j).anas.entropy.subint(l).signal(n).format(m) = true;
                end

            end
        end
    end
end
if ~isempty(figs(1).name)
    lstbxstr = {[]};
    for i = 1:length(figs)
        for j = 1: length(figs(i).ints)
            if isempty(lstbxstr{1})
                lstbxstr{1} = ['[' num2str(i) '-' num2str(j) '] ' figs(i).name ' ' figs(i).ints(j).name ];
            else
                lstbxstr{end+1} = ['[' num2str(i) '-' num2str(j) '] ' figs(i).name ' ' figs(i).ints(j).name ];
            end
        end
    end
else
    lstbxstr = {'no result files found in selected folder'};
end
app.lstbx_report.Items = lstbxstr;


end