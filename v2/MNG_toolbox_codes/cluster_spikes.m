
function [cluster_res] = cluster_spikes(index, spikes, par,threshold)

par_inp = par;     %%

min_spikes4SPC = 16; % if are less that this number of spikes, clustering won't be made.

%default config
par_input = struct;
parallel = false;
make_times = true;
make_plots = true;%make_plots = true
resolution = '-r150';
save_spikes = true;
%search for optional inputs
par_input = par;
run_par_for = parallel; 
current = pwd;
par_input.folder = GetWritableFolder;
cd(par_input.folder)
input = 'temp_spikes.mat';%input = [ par_input.folder '\temp.mat'];

save(input, 'par','spikes','index','threshold')% save(input, 'par','spikes','index','psegment','sr_psegment','threshold')

filenames = {input}; 
tic
par_file = set_parameters(par.sr);


initial_date = now;
Nfiles = length(filenames);


filename = filenames{1};

par = struct;
par = update_parameters(par,par_file,'clus');
par = update_parameters(par,par_file,'batch_plot');
par = update_parameters(par,par_input,'clus');
par = update_parameters(par,par_input,'batch_plot');
par.filename = filename;
par.reset_results = true;
par.sr = par_input.sr;
par.segments_length = par_input.segments_length;
par.tmax = par_input.tmax;
par.tmin = par_input.tmin;
par.w_pre = par_input.w_pre;
par.w_post = par_input.w_post;
par.alignment_window = par_input.alignment_window;
par.stdmin = par_input.stdmin;
par.stdmax = par_input.stdmax;
par.detect_fmin = par_input.detect_fmin;
par.detect_fmax = par_input.detect_fmax;
par.sort_fmin = par_input.sort_fmin;
par.sort_fmax = par_input.sort_fmax;
par.ref_ms = par_input.ref_ms;
par.detection =  par_input.detection;
par.int_factor = par_input.int_factor;
par.interpolation = par_input.interpolation;
par.sort_order = par_input.sort_order;
par.detect_order = par_input.detect_order;
par.detection_date = par_input.detection_date;
par.ref = floor(par.ref_ms*par.sr/1000);

check_WC_params(par)

par.fname_in = ['tmp_data_wc' num2str(1)];%par.fname_in = [par.folder '\tmp_data_wc' num2str(fnum)];                       % temporary filename used as input for SPC
[~, par.nick_name, ~] = fileparts(par.filename);
par.nick_name = par.nick_name(1:end-7);
par.fname = ['data_' par.nick_name]; %par.fname = [par.folder '\data_' data_handler.nick_name];

par.fnamespc = ['data_wc' num2str(1)];%par.fnamespc = [par.folder '\data_wc' num2str(fnum)];


% LOAD SPIKES
nspk = size(spikes,1);
naux = min(par.max_spk,size(spikes,1));
disp('after loading')                                                                           %%
disp([num2str(nspk) ' spikes ' num2str(min_spikes4SPC) ' needed'])                              %%
disp(['size spikes: ' num2str(size(spikes)) ' size index: ' num2str(size(index)) ' needed'])    %%
if nspk < min_spikes4SPC
    warning('MyComponent:noValidInput', 'Not enough spikes in the file');
    return
end

% CALCULATES INPUTS TO THE CLUSTERING ALGORITHM.
inspk = wave_features(spikes,par);     %takes wavelet coefficients.
par.inputs = size(inspk,2);                       % number of inputs to the clustering

if par.permut == 'n'
    % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
    if size(spikes,1)> par.max_spk
        % take first 'par.max_spk' spikes as an input for SPC
        inspk_aux = inspk(1:naux,:);
    else
        inspk_aux = inspk;
    end
else
    % GOES FOR TEMPLATE MATCHING IF TOO MANY SPIKES.
    if size(spikes,1)> par.max_spk;
        % random selection of spikes for SPC
        ipermut = randperm(length(inspk));
        ipermut(naux+1:end) = [];
        inspk_aux = inspk(ipermut,:);
    else
        ipermut = randperm(size(inspk,1));
        inspk_aux = inspk(ipermut,:);
    end
end
%INTERACTION WITH SPC
save(par.fname_in,'inspk_aux','-ascii');
pwd
bla = dir;
for i = 3:length(bla)
    disp([bla(i).name ' ' num2str(bla(i).bytes/1000) 'kB'])
end
try
    [clu, tree] = run_cluster_loc(par,true);
    if exist([par.fnamespc '.dg_01.lab'],'file')
        movefile([par.fnamespc '.dg_01.lab'], [par.fname '.dg_01.lab'], 'f');
        movefile([par.fnamespc '.dg_01'], [par.fname '.dg_01'], 'f');
    end
catch ME
    warning('MyComponent:ERROR_SPC', 'Error in SPC');
    return
end
cluster_res.clu = clu;
cluster_res.tree = tree;
[clust_num temp auto_sort] = find_temp(tree,clu,par);

if par.permut == 'y'
    clu_aux = zeros(size(clu,1),2 + size(spikes,1)) -1;  %when update classes from clu, not selected elements go to cluster 0
    clu_aux(:,ipermut+2) = clu(:,(1:length(ipermut))+2);
    clu_aux(:,1:2) = clu(:,1:2);
    clu = clu_aux;
    clear clu_aux
end

classes = zeros(1,size(clu,2)-2);
for c =1: length(clust_num)
    aux = clu(temp(c),3:end) +1 == clust_num(c);
    classes(aux) = c;
end

if par.permut == 'n'
    classes = [classes zeros(1,max(size(spikes,1)-par.max_spk,0))];
end

Temp = [];
% Classes should be consecutive numbers
classes_names = nonzeros(sort(unique(classes)));
for i= 1:length(classes_names)
   c = classes_names(i);
   if c~= i
       classes(classes == c) = i;
   end
   Temp(i) = temp(i);
end

% IF TEMPLATE MATCHING WAS DONE, THEN FORCE
if (size(spikes,1)> par.max_spk || ...
        (par.force_auto))
    f_in  = spikes(classes~=0,:);
    f_out = spikes(classes==0,:);
    class_in = classes(classes~=0);
    class_out = force_membership_wc(f_in, class_in, f_out, par);
    forced = classes==0;
    classes(classes==0) = class_out;
    forced(classes==0) = 0;
else
    forced = zeros(1, size(spikes,1));
end

gui_status = struct();
gui_status.current_temp =  max(temp);
gui_status.auto_sort_info = auto_sort;
gui_status.original_classes = zeros(size(classes));

for i=1:max(classes)
    gui_status.original_classes(classes==i) = clust_num(i);
end

current_par = par;
par = struct;
par = update_parameters(par, current_par, 'relevant');
par = update_parameters(par,current_par,'batch_plot');

par.sorting_date = datestr(now);
cluster_class = zeros(nspk,2);
cluster_class(:,2)= index';
cluster_class(:,1)= classes';

vars = {'cluster_class', 'par','inspk','forced','Temp','gui_status'};
cluster_res.cluster_class = cluster_class;
cluster_res.par = par;
cluster_res.inspk = inspk;
cluster_res.forced = forced;
cluster_res.Temp = Temp;
cluster_res.gui_status = gui_status;

if exist('ipermut','var')
    vars{end+1} = 'ipermut';
    cluster_res.ipermut = ipermut;
end
if save_spikes
    vars{end+1} = 'spikes';
else
    spikes_file = filename;
    vars{end+1} = 'spikes_file';
end

tocaux = toc;
disp(['Computations Done (' num2str(tocaux,'%2.2f') 's).'])


tmp = dir;
for i = 3:length(tmp)
    delete (tmp(i).name)
end
cd(current)
toc

end



