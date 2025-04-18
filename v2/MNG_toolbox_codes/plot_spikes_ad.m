function plot_spikes_ad(app)

drawnow;

data = app.user_data;


par = data{1};
spikes = data{2};
classes = data{6};
classes = classes(:)';
inspk = data{7};
temp = data{8};
ls = size(spikes,2);
minclus = app.settings.minclus;
clustering_results = data{10};

% Extract spike features if needed

if strcmp(app.dd_shape.Value, app.dd_shape.Items{1})
    if isempty(inspk) || (length(inspk)~=size(spikes,1))
        [inspk] = wave_features(spikes,app.settings.par);
        data{7} = inspk;
    end
end

% Classes should be consecutive numbers
classes_names = sort(unique(classes));
classes_names = classes_names(classes_names>0);

% updates 'clustering_results_bk'
data{11} = clustering_results; 

% Forcing
if app.settings.force==1
    for i = classes_names
        ind = find(clustering_results(:,2)==i); % get index of GUI class
        oclass = clustering_results(ind(1),4); % get original class
        otemp = clustering_results(ind(1),3); % get original temperature
        ind2 = find(classes==i); % get index of forced class
        clustering_results(ind2,3) = otemp; % update original temperatures with forced class
        clustering_results(ind2,4) = oclass; % update original class with forced class
    end
end

for i = 1:length(classes_names)
   c = classes_names(i);
   if c~= i
       classes(classes == c) = i;
   end
end

% Defines nclusters
cluster_sizes = zeros(1,33);
ifixflag = zeros(1,33);
for i=1:size(cluster_sizes,2)
    cluster_sizes(i) = nnz(classes==i);
end

if app.settings.setclus == 0 && app.settings.undo==0 && app.settings.merge==0 && app.settings.force==0  
    sizemin_clus = minclus;
elseif app.settings.setclus == 1
	sizemin_clus = 1;
else
    sizemin_clus = minclus;
end

clusn = find(cluster_sizes >= sizemin_clus);
nclusters = length(clusn);

% Get fixed clusters
fix_class2 = [];
nfix_class = [];

if ~isfield(app.settings,'new_manual')
	for i = 1 : 33
        if app.settings.cluster(i).fix && ~isempty(data{19+i})
            nclusters = nclusters +1;
            fix_class = data{19+i}';
            classes(classes==nclusters)=0;
            classes(fix_class)=nclusters;
            ifixflag(nclusters) = 1;
			fix_class2 = [fix_class2 fix_class];
			nfix_class = [nfix_class i];
			clusn = [clusn nclusters];
        end
    end
end

% Merge operations
if app.settings.merge == 1 && ~isempty(nfix_class)
    imerge = 1;% index for the original temperature that will represent all the fixed classes
    bigger_fix = 0;
    for i = 1:length(nfix_class)
        aux = nnz(clustering_results(:,2) == nfix_class(i));
        if aux > bigger_fix
            imerge = i;
            bigger_fix = aux;
        end
    end
    imerge = find(clustering_results(:,2) ==nfix_class(imerge),1);
    mtemp = clustering_results(imerge,3); % temperature that represents all the fixed classes
    classes(fix_class2) = nfix_class(1); % labels all the fixed classes with the new number
end

if ~isfield(app.settings,'unforce')
    if app.settings.force==0  &&  app.settings.setclus==0
        forced = data{13};
        data{14} = forced;
        new_forced = false(size(forced));
        new_forced(fix_class2) = forced(fix_class2);
        clear forced
        data{13} = new_forced;
    elseif app.settings.force==0
        data{14} = data{13};
    end
end
% Defines classes
non_clustered = ones(1,size(spikes,1));
nclusters = 0;
for i = clusn
    class_temp = find(classes == i);
    if ((ifixflag(i)==1) && (~isempty(class_temp)))
        ifixflagc = 1;
    else
        ifixflagc = 0;
    end
    if ((length(class_temp) >= sizemin_clus) || (ifixflagc == 1))
        nclusters = nclusters+1;
        eval(['class' num2str(nclusters) '= class_temp;'])
        non_clustered(class_temp) = 0;
    end
end
rejected = data{15};
class0 = find(non_clustered & ~rejected);
clear non_clustered rejected

% Redefines classes
classes = zeros(size(spikes,1),1);
for i = 0:nclusters
    if ~ (isempty(class0) && i==0)
        eval(['classes(class' num2str(i) ') = ' num2str(i) ';']);
    end
end

% Saves new classes
data{6} = classes;


% new temperature when merge
if app.settings.merge == 1 && ~isempty(nfix_class)
    clustering_results(fix_class2,3) = mtemp;
    clustering_results(fix_class2,4) = clustering_results(imerge,4);
end 
clustering_results(:,1) = temp; % GUI temperature
clustering_results(:,5) = minclus; % GUI minimum cluster
 
% Saves new classes and keep fixed classes in 'clustering_results'. 
% Keep the original temperature and cluster number in the fixed spikes.
% The temperature of the non-fixed spikes will be 
% the GUI temperature (temp) and cluster number will be 
% the GUI cluster number (classes)
if (~isempty(fix_class2)) && app.settings.merge==0 && app.settings.undo==0 && app.settings.force==0
    % selects the index of the non-fixed spikes
    % since those are the ones which are going to be updated
    ind_non_fix = 1:length(classes); 
    ind_non_fix(fix_class2) = []; 
    if isfield(app.settings,'new_spc_classes')
            clustering_results(ind_non_fix,4) = app.settings.new_spc_classes(ind_non_fix);
    end
    if app.settings.setclus == 0
        clustering_results(ind_non_fix,3) = temp; % temperature of the non-fixed spikes in the original temperature column
    end
end

% update new classes
clustering_results(:,2) = classes;

% If there are no fix and rejected clusters and undo operations, 
% original classes are the same as current classes
% or 0 if they are rejected or manual selected
if isempty(fix_class2) && app.settings.undo==0 && app.settings.merge==0 && app.settings.force==0
    if isfield(app.settings,'new_manual')
		clustering_results(:,4) = clustering_results(:,4); % same as before
		clustering_results(:,3) = clustering_results(:,3);
		clustering_results(app.settings.new_manual|classes==0,4) = 0;
		clustering_results(app.settings.new_manual,3) = temp;
	elseif app.settings.setclus == 0
        if isfield(app.settings,'new_spc_classes')
            clustering_results(:,4) = app.settings.new_spc_classes;
        else
            clustering_results(:,4) = clustering_results(:,2); % clusters
        end
		clustering_results(:,3) = temp; % temperatures
    end
end

clear classes
% Updates clustering_results in USER_DATA
data{10} = clustering_results; 

for i=20:52
    data{i} = [];
end

for i=1:33
    app.settings.cluster(i).fix = false;
    app.settings.cluster(i).isclu = false;
end

%% plot cluster spikes

app.user_data = data;
cla(app.ax_clu0_isi)
cla(app.ax_clu0,'reset')
cla(app.ax_clusters)
cla(app.ax_clus_a)
cla(app.ax_clus_b)
cla(app.ax_clus_c)
cla(app.ax_clus_d)
cla(app.ax_clus_a_isi)
cla(app.ax_clus_b_isi)
cla(app.ax_clus_c_isi)
cla(app.ax_clus_d_isi)
hold(app.ax_clusters, 'on')
ylimit = [nan,nan];
colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
    [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
    [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
maxc = size(colors,1);

forced = data{13};

for i = 0 : nclusters
    tmpy  = []; %as a flag to don't make the same vector twice
    if ~ (isempty(class0) && i==0)
        class_i = eval(['class' num2str(i)]);
        sup_spikes = length(class_i);
        max_spikes = min(sup_spikes, par.max_spikes_plot);
        permut = randperm(sup_spikes);
        permut = permut(1:max_spikes);
        xlim(app.ax_clusters,'manual')
        if find(contains(app.dd_shape.Items,app.dd_shape.Value)) == 1 && find(contains(app.dd_plot.Items,app.dd_plot.Value)) == 1 && ~strcmp(par.all_classes_ax, 'mean') %supposed for batchplotting?
            tmpy=spikes(class_i(permut),:);                                                                             %
            tmpy=spikes(class_i(permut),:);                                                                             % need
            tmpn=size(tmpy,1);                                                                                          % to
            tmpx=repmat([1:ls NaN]',1,tmpn);                                                                            % be
            tmpx=reshape(tmpx,numel(tmpx),1);                                                                           % checked
            tmpy=[tmpy'; repmat(NaN,1,tmpn)];                                                                           %
            tmpy=reshape(tmpy,numel(tmpy),1);                                                                           %
            line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',app.ax_clusters,'Visible','off');      %
			xlim(app.ax_clusters, [1 ls])                                                                           %
        elseif find(contains(app.dd_shape.Items,app.dd_shape.Value)) == 1
            av = mean(spikes(class_i,:),1);
            plot(app.ax_clusters,1:ls,av,'Color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2);
            xlim(app.ax_clusters,[1 ls])
        else
            plot(app.ax_clusters,inspk(class_i,1),inspk(class_i,2),'.','Color',colors(mod(i-1,maxc)+1,:)*(i~=0),'markersize',.5);
            axis(app.ax_clusters,'auto');
        end
        
        av = mean(spikes(class_i,:),1);
        avup = av + par.to_plot_std * std(spikes(class_i,:),0,1);
        avdw = av - par.to_plot_std * std(spikes(class_i,:),0,1);
        ylimit(1) = min(ylimit(1),min(avdw));
        ylimit(2) = max(ylimit(2),max(avup));
        if find(contains(app.dd_plot.Items,app.dd_plot.Value)) == 1
            % optimizing for speed:
            if isempty(tmpy)
                tmpy=spikes(class_i(permut),:);
                tmpn=size(tmpy,1);
                tmpx=repmat([1:ls NaN]',1,tmpn);
                tmpx=reshape(tmpx,numel(tmpx),1);
                tmpy=[tmpy'; repmat(NaN,1,tmpn)];
                tmpy=reshape(tmpy,numel(tmpy),1);
            end
            if i == 0
                line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',app.ax_clu0);
                line(1:ls,av,'color','c','linewidth',2,'Parent',app.ax_clu0)
                line(1:ls,avup,'color','c','linewidth',0.5,'Parent',app.ax_clu0)
                line(1:ls,avdw,'color','c','linewidth',0.5,'Parent',app.ax_clu0)
            else
                app.settings.cluster(i).av = av;
                app.settings.cluster(i).avup = avup;
                app.settings.cluster(i).avdw = avdw;
                app.settings.cluster(i).spikes = [tmpx,tmpy];
            end
        else
            if i == 0
                plot(app.ax_clu0,1:ls,av,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2)
                plot(app.ax_clu0,1:ls,avup,1:ls,avdw,'color',[.65 .65 .65],'linewidth',.5)
            else
                app.settings.cluster(i).av = av;
                app.settings.cluster(i).avup = avup;
                app.settings.cluster(i).avdw = avdw;
                app.settings.cluster(i).spikes = [];
            end
        end
        eval(['aux=num2str(length(class' num2str(i) '));']);
        if i == 0
            app.lbl_clu0.Text =['Cluster ' num2str(i) ':  # ' aux];
            xlim(app.ax_clu0, [1 ls])
        else
            app.settings.cluster(i).ttl = ['Cluster ' num2str(i) ':  # ' aux ' (' num2str(nnz(clustering_results(:,2)==i & ~forced(:))) ')'];
%             app.settings.cluster(i).ylimit = [ylimit; ylim(clus_ax)]; %manually link y-axes
        end
    end
end
app.settings.cluster_0.yl = [ylimit(1)-(diff(ylimit)*0.2), ylimit(2)+(diff(ylimit)*0.2)];
ylim(app.ax_clu0, app.settings.cluster_0.yl)
app.settings.cluster_0.n_clus = nclusters;
%% plot isi

spk_times = data{3};
classes = data{6};
for i = 0:nclusters
    if classes == 0
        rejected = data{15};
        times = diff(spk_times(classes(:)==i & ~rejected(:)));
        clear rejected 
    else
        times = diff(spk_times(classes==i));
        if i ~=0
            app.settings.cluster(i).isclu = true;
        end
    end
    % Calculates # ISIs < 3ms  
    multi_isi = nnz(times < 3); 
    % Builds and plots the histogram
    if i == 0
        bin_step = app.settings.cluster_0.bin_step;
        nbins = app.settings.cluster_0.nbins;
    else
        bin_step = app.settings.cluster(i).bin_step;
        nbins = app.settings.cluster(i).nbins;
    end

    [N,X]=hist(times,0:bin_step:nbins);

    ttl = [num2str(multi_isi) ' in < 3ms'];
    if i == 0
        xlim (app.ax_clu0_isi,'manual')
        bar(app.ax_clu0_isi,X(1:end-1),N(1:end-1))
        xlim(app.ax_clu0_isi,[0 nbins]);
        app.lbl_clus0_isi.Text = ttl;
    else
        app.settings.cluster(i).N = N;
        app.settings.cluster(i).X = X;
        app.settings.cluster(i).ttl_isi = ttl;
    end
end

%% plot temp 

if ~isempty(data{5})
    tree = data{5};
    app.settings.par.min.clus = clustering_results(1,5);
    nclasses = max(clustering_results(:,2));
    class_plot = [];
    for i=1:nclasses
        ind = find(clustering_results(:,2)==i);
        classgui_plot(i) = clustering_results(ind(1),2);
        class_plot(i) = clustering_results(ind(1),4);
        if class_plot(i) == 0 %null original cluster
	        class_plot(i) =1; %plot like they were from cluster 1
        end
        temp_plot(i) = clustering_results(ind(1),3); 
    end
    
    num_temp = floor((app.settings.par.maxtemp -app.settings.par.mintemp)/app.settings.par.tempstep);     % total number of temperatures 
    
    tree(num_temp+1,2) = app.settings.par.mintemp+(num_temp)*app.settings.par.tempstep; %added for handle selection of max temp
    
    temperature = tree(clustering_results(1,1)+1,2);
    
    colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
        [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
        [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
    maxc = size(colors,1);
    auto_sort_info = app.autosort_nfo;
    % draw temperature diagram and mark clusters 
    cla(app.ax_temp);
    
    switch app.settings.par.temp_plot
        case 'lin'
            % draw diagram
            hold(app.ax_temp, 'on');
            if ~isempty(auto_sort_info)
                [xp, yp] = find(auto_sort_info.peaks);
                ptemps = app.settings.par.mintemp+(xp)*app.settings.par.tempstep;
                psize = tree(sub2ind(size(tree), xp,yp+4));
                plot(app.ax_temp,ptemps,psize,'xk','MarkerSize',7,'LineWidth',0.9);
                area(app.ax_temp,app.settings.par.mintemp+app.settings.par.tempstep.*[auto_sort_info.elbow,num_temp],max(ylim(app.ax_temp)).*[1 1],'LineStyle','none','FaceColor',[0.9 0.9 0.9]);
            end
            plot(app.ax_temp, [app.settings.par.mintemp, app.settings.par.maxtemp-app.settings.par.tempstep],[app.settings.par.min.clus2, app.settings.par.min.clus2],'k:',...
                app.settings.par.mintemp+(1:num_temp)*app.settings.par.tempstep, ...
                tree(1:num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
            % mark clusters
            for i=1:min(size(tree,2)-4,length(class_plot))
                tree_clus = tree(temp_plot(i),4+class_plot(i));
                tree_temp = tree(temp_plot(i)+1,2);
                plot(app.ax_temp, tree_temp,tree_clus,'.','color',colors(mod(classgui_plot(i)-1,maxc)+1,:),'MarkerSize',20);
            end
            set(get(app.ax_temp,'ylabel'),'vertical','Baseline');%set(get(gca,'ylabel'),'vertical','Baseline');
        case 'log'
            % draw diagram
            set(app.ax_temp,'yscale','log');
            hold(app.ax_temp, 'on');
            if ~isempty(auto_sort_info)
                [xp, yp] = find(auto_sort_info.peaks);
                ptemps = app.settings.par.mintemp+(xp)*app.settings.par.tempstep;
                psize = tree(sub2ind(size(tree), xp,yp+4));
                semilogy(app.ax_temp,ptemps,psize,'xk','MarkerSize',7,'LineWidth',0.9);
                area(app.ax_temp,app.settings.par.mintemp+app.settings.par.tempstep.*[auto_sort_info.elbow,num_temp],max(ylim(app.ax_temp)).*[1 1],'LineStyle','none','FaceColor',[0.9 0.9 0.9],'basevalue',1);
            end
            semilogy(app.ax_temp, [app.settings.par.mintemp app.settings.par.maxtemp-app.settings.par.tempstep], ...
                [app.settings.par.min.clus app.settings.par.min.clus],'k:',...
                app.settings.par.mintemp+(1:num_temp)*app.settings.par.tempstep, ...
                tree(1:num_temp,5:size(tree,2)),[temperature temperature],[1 tree(1,5)],'k:')
            % mark clusters
            for i=1:length(class_plot)
                if class_plot(i)+4>size(tree,2)
                    continue
                end
                tree_clus = tree(temp_plot(i),4+class_plot(i));
                tree_temp = tree(temp_plot(i)+1,2);
                semilogy(app.ax_temp, tree_temp,tree_clus,'.','color',colors(mod(classgui_plot(i)-1,maxc)+1,:),'MarkerSize',20);
            end
            set(get(app.ax_temp,'ylabel'),'vertical','Baseline');% set(get(handles.temperature_plot,'ylabel'),'vertical','Baseline');
    end
    
    % xlim(handles.temperature_plot, [0 handles.par.maxtemp])
    xlabel(app.ax_temp, 'Temperature','FontSize',8); 
    ylabel(app.ax_temp, 'Clusters size','FontSize',8);
    
    set(allchild(app.ax_temp),'Visible','on')
end 
    
%Resize axis
if ~isempty(ylimit)
    
  
    ylim(app.ax_clu0,app.settings.cluster_0.yl);
    if find(contains(app.dd_shape.Items,app.dd_shape.Value)) ==1
        ylim(app.ax_clusters,app.settings.cluster_0.yl);
    end
%     linkaxes(ax_v,'xy'); %drawnow inside
%     ylim(ax_v(1),[ymin ymax]);
else
    drawnow
end

update_cluster_plots(app)



end









        
        
%         %% apply all changes to cluster axis,btns,etc
% %update cluster label
% update_cluster_plots(app) 
% ylimit = [];
% colors = [[0.0 0.0 1.0];[1.0 0.0 0.0];[0.0 0.5 0.0];[0.620690 0.0 0.0];[0.413793 0.0 0.758621];[0.965517 0.517241 0.034483];
%     [0.448276 0.379310 0.241379];[1.0 0.103448 0.724138];[0.545 0.545 0.545];[0.586207 0.827586 0.310345];
%     [0.965517 0.620690 0.862069];[0.620690 0.758621 1.]]; 
% maxc = size(colors,1);
% 
% 
% if i>0 
%                 ylim(clus_ax,'auto');
%                 ylimit = [ylimit; ylim(clus_ax)];
%                 title(clus_ax,['Cluster ' num2str(i) ':  # ' aux ' (' num2str(nnz(clustering_results(:,2)==i & ~forced(:))) ')'],'Fontweight','bold');
%             else            
% %%
% 
% for i = 0:nclusters
%     tmpy = []; %as a flag to don't make the same vector twice
%     if ~ (isempty(class0) && i==0)
%         %PLOTS SPIKES OR PROJECTIONS                    %% 
%         class_i = eval(['class' num2str(i)]);
%         sup_spikes = length(class_i);
%         max_spikes = min(sup_spikes, par.max_spikes_plot);
%         permut = randperm(sup_spikes);
%         permut = permut(1:max_spikes);
%         xlim(handles.projections,'manual');
%         if get(handles.spike_shapes_button,'value') ==1 && (get(handles.plot_all_button,'value') ==1) && ~strcmp(par.all_classes_ax,'mean')
%             % optimizing for speed:
%             tmpy=spikes(class_i(permut),:);
%             tmpn=size(tmpy,1);
%             tmpx=repmat([1:ls NaN]',1,tmpn);
%             tmpx=reshape(tmpx,numel(tmpx),1);
%             tmpy=[tmpy'; repmat(NaN,1,tmpn)];
%             tmpy=reshape(tmpy,numel(tmpy),1);
%             line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',handles.projections,'Visible','off');
% 			xlim(handles.projections, [1 ls])
%         elseif get(handles.spike_shapes_button,'value') ==1
%             av   = mean(spikes(class_i,:));
%             plot(handles.projections,1:ls,av,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2); %% Plot cluster shape
%             xlim(handles.projections,[1 ls])
%         else
%             plot(handles.projections,inspk(class_i,1),inspk(class_i,2),'.','Color',colors(mod(i-1,maxc)+1,:)*(i~=0),'markersize',.5);
%             axis(handles.projections,'auto');
%         end
%                                                                                     
%                                     %% plot cluster
%         if i < 4
%             clus_ax = eval(['handles.spikes' num2str(i)]); 
%             xlim(clus_ax,'manual');
%             xlim(clus_ax,[1 ls]);
%             hold(clus_ax,'on')
%             
%             av   = mean(spikes(class_i,:));
%             avup = av + par.to_plot_std * std(spikes(class_i,:));
%             avdw = av - par.to_plot_std * std(spikes(class_i,:));
%                       
%             if get(handles.plot_all_button,'value') ==1
%                 % optimizing for speed:
%                 if isempty(tmpy)
%                     tmpy=spikes(class_i(permut),:);
%                     tmpn=size(tmpy,1);
%                     tmpx=repmat([1:ls NaN]',1,tmpn);
%                     tmpx=reshape(tmpx,numel(tmpx),1);
%                     tmpy=[tmpy'; repmat(NaN,1,tmpn)];
%                     tmpy=reshape(tmpy,numel(tmpy),1);
%                 end
%                 line(tmpx,tmpy,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'Parent',clus_ax); %% plot spikes
% 				
% 				if i==0
%                     line(1:ls,av,'color','c','linewidth',2,'Parent',clus_ax)
%                     line(1:ls,avup,'color','c','linewidth',0.5,'Parent',clus_ax)
%                     line(1:ls,avdw,'color','c','linewidth',0.5,'Parent',clus_ax)
%                 else
%                     line(1:ls,av,'color','k','linewidth',2,'Parent',clus_ax) % plot mean spike
%                     line(1:ls,avdw,'color',[.4 .4 .4],'linewidth',0.5,'Parent',clus_ax) % plot std+
%                     line(1:ls,avup,'color',[.4 .4 .4],'linewidth',0.5,'Parent',clus_ax) % plot std-
%                 end
%             else
%                 plot(clus_ax,1:ls,av,'color',colors(mod(i-1,maxc)+1,:)*(i~=0),'linewidth',2)
%                 plot(clus_ax,1:ls,avup,1:ls,avdw,'color',[.65 .65 .65],'linewidth',.5)
%             end
%             
%             eval(['aux=num2str(length(class' num2str(i) '));']);
%             if i>0 
%                 ylim(clus_ax,'auto');
%                 ylimit = [ylimit; ylim(clus_ax)];
%                 title(clus_ax,['Cluster ' num2str(i) ':  # ' aux ' (' num2str(nnz(clustering_results(:,2)==i & ~forced(:))) ')'],'Fontweight','bold');
%             else            
%                 title(clus_ax,['Cluster ' num2str(i) ':  # ' aux],'Fontweight','bold');
%                 xlim(clus_ax, [1 ls])
%             end
%             
%         else
%             par.axes_nr = i+1;
%             par.ylimit = ylimit;
%             eval(['par.class_to_plot = class' num2str(i) ';']);
%             par.plot_all_button = get(handles.plot_all_button,'value');
%             USER_DATA{1} = par;
%             set(handles.wave_clus_figure,'userdata',USER_DATA)
% 
%             if i < 9
%                 opened_figs{1} = wave_clus_aux('Visible', 'off');
%             elseif i < 14
%                 opened_figs{2} = wave_clus_aux1('Visible', 'off');
%             elseif i < 19
%                 opened_figs{3} = wave_clus_aux2('Visible', 'off');
%             elseif i < 24
%                 opened_figs{4} = wave_clus_aux3('Visible', 'off');
%             elseif i < 29
%                 opened_figs{5} = wave_clus_aux4('Visible', 'off');
%             elseif i < 34
%                 opened_figs{6} = wave_clus_aux5('Visible', 'off');
%             %-------------------------------------------------------------------------
%             end
%         end
%     end
% end
% 
% 
% draw_histograms(handles, 0:min(nclusters,3),USER_DATA); %% plot ISI
% 
% if ~isempty(USER_DATA{5})
%     mark_clusters_temperature_diagram(handles,USER_DATA{5},clustering_results) %% plot temp graph
% end
% set(handles.file_name,'string', par.file_name_to_show);
% 
% set(allchild(handles.projections),'Visible','on')
% 
% %Resize axis
% if ~isempty(ylimit)
%     ymin = min(ylimit(:,1));
%     ymax = max(ylimit(:,2));
%     ylim(handles.spikes0,[ymin ymax]);
%     if get(handles.spike_shapes_button,'value') ==1
%         ylim(handles.projections,[ymin ymax]);
%     end
%     linkaxes(ax_v,'xy'); %drawnow inside
%     ylim(ax_v(1),[ymin ymax]);
% else
%     drawnow
% end
% 
% for i =1:figs_num
%     if ~isempty(opened_figs{i})  
%         set(opened_figs{i},'units','normalized','outerposition',[0 0 1 1])
%         set(opened_figs{i},'Visible', 'on'); 
%     end
% end
% 
% if exist('groot','builtin')
%     set(groot,'defaultfiguregraphicssmoothing','remove')
%     set(groot,'DefaultAxesFontSize','remove')
% end
