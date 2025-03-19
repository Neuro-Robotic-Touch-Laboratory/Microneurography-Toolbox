function [brst_files, int_templ, path] = find_xls(app)

path = uigetdir(pwd, 'please select master folder that contain all *.xls files to combine');
files = dir(fullfile(path,'**','*.xls'));
tmp = string([]); 
for i = 1: length(files)
    tmp(i)= string(files(i).name);
end
[~, tmp_idx] = sort(tmp);
files(tmp_idx) = files;
[file, path] = uigetfile('templ_name.mat','Please select the templ_name.mat file');
tmp = load([path  file]);
int_templ = tmp.int_templates;
brst_files = struct('recording',[],'interval',[],'filename',[], 'subintervals', {[]},'signals',{[]});
strc_idx = 1;
lb_str = {[]};

sint_str = {[]};
sigs_str = {[]};

for i = 1:length(files)
    if ~isempty(strfind(files(i).name, '_burst_results.xls'))
        %skip = false;%%
        tmp_sints = {[]};
        brst_files(strc_idx).filename = [files(i).folder '\' files(i).name];
        tmp_idx = strfind(files(i).name, '_burst_results.xls');
        tmp = files(i).name(1:tmp_idx-1);
        us_idx = strfind(tmp ,'_INT_');
        brst_files(strc_idx).recording = tmp(1:us_idx-1);
        brst_files(strc_idx).interval = tmp(us_idx+5:end);
        lb_str{strc_idx} = tmp;
        
        tmp_cell = readcell(brst_files(strc_idx).filename,'Sheet', 'general');
        tmp_sint_cell = tmp_cell(1,2:end);
        tmp_sig_cell = tmp_cell(4:end,1);
        is_tmpl = zeros(size(tmp_sint_cell));
        for j = 1 : length(tmp_sint_cell)
            if find(strcmp(tmp_sint_cell{j},int_templ))
                if isempty(tmp_sints{1,1})
                    tmp_sints{1,1} = tmp_sint_cell{j};
                    tmp_sints{1,2} = j+1;
                else
                    tmp_sints{end+1,1} = tmp_sint_cell{j};
                    tmp_sints{end,2} = j+1;
                end
            else
                if ~isempty(strfind(tmp_sint_cell{j},'§$1>')) || ~isempty(strfind(tmp_sint_cell{j},'§$2>'))
                    if isempty(tmp_sints{1,1})
                        tmp_sints{1,1} = tmp_sint_cell{j};
                        tmp_sints{1,2} = j+1;
                    else
                        tmp_sints{end+1,1} = tmp_sint_cell{j};
                        tmp_sints{end,2} = j+1;
                    end
                end
            end
        end
        tmp_sig = {[]};
        for j = 1:length(tmp_sig_cell)
            if isempty(tmp_sig{1,1})
                tmp_sig{1,1} = tmp_sig_cell{j};
                tmp_sig{1,2} = j+3;
            else
                tmp_sig{end+1,1} = tmp_sig_cell{j};
                tmp_sig{end,2} = j+3;
            end
        end
    
        brst_files(strc_idx).subintervals = tmp_sints;
        brst_files(strc_idx).signals = tmp_sig;
        if isempty(sint_str{1})
            sint_str = tmp_sints(:,1);
        else
            %disp(num2str(i))
            if ~isempty(tmp_sints{1,1})
                sint_str = unique(vertcat(sint_str, tmp_sints(:,1)),'stable');
                %skip = true;%%
            end
        end
        %if ~skip
            if isempty(sigs_str{1})
                sigs_str = tmp_sig(:,1);
            else
                sigs_str = unique(vertcat(sigs_str, tmp_sig(:,1)),'stable');
            end
            strc_idx = strc_idx+1;
        %end
    end
end

app.lstbx_rec_combine.Items = lb_str;
app.lstbx_intervals.Items = sint_str; 
app.lstbx_signals.Items = sigs_str; 
end