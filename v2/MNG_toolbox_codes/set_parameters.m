function par = set_parameters(fs) 
 
% LOAD PARAMS 
par.segments_length = 5;        % length (in minutes) of segments in which the data is cutted (default 5min). 
 
% PLOTTING PARAMETERS 
par.cont_segment = true; 
par.max_spikes_plot = 1000;     % max. number of spikes to be plotted 
par.print2file = true;          % If is not true, print the figure (only for batch scripts). 
par.cont_plot_samples = 100000; % number of samples used in the one-minute (maximum) sample of continuous data to plot. 
par.to_plot_std = 1;            % # of std from mean to plot 
par.all_classes_ax = 'mean';    % 'mean'/'all'. If it's 'mean' only the mean waveforms will be ploted in the axes with all the classes 
par.plot_feature_stats = false; 

% SPC PARAMETERS 
par.mintemp = 0.00;             % minimum temperature for SPC 
par.maxtemp = 0.251;            % maximum temperature for SPC 
par.tempstep = 0.01;            % temperature steps 
par.SWCycles = 100;             % SPC iterations for each temperature (default 100) 
par.KNearNeighb = 11;           % number of nearest neighbors for SPC 
par.min_clus = 20;              % minimum size of a cluster (default 60) 
par.randomseed = 0;             % if 0, random seed is taken as the clock value (default 0) 
par.temp_plot = 'log';          % temperature plot in log scale ('lin', 'log')
 
par.c_ov = 0.7;                 % Overlapping coefficient to use for the inclusion criterion. 
par.elbow_min  = 0.4;           % Thr_border parameter for regime border detection. 
 
% DETECTION PARAMETERS 
par.tmax = 'all';               % maximum time to load 
%par.tmax= 180;                 % maximum time to load (in sec) 
par.tmin= 0;                    % starting time for loading (in sec) 
par.w_pre = ceil(0.002*fs);%par.w_pre = 20;                 % number of pre-event data points stored (default 20) 
par.w_post = ceil(0.002*fs);%par.w_post = 20;                % number of post-event data points stored (default 44)) 
par.alignment_window = 10;      % number of points around the sample expected to be the maximum 
par.stdmin =4;                  % minimum threshold for detection %3.5 
par.stdmax =10;                 % maximum threshold for detection       %10 
par.detect_fmin = 300;          % high pass filter for detection 
par.detect_fmax = 3000;         % low pass filter for detection (default 1000) 
par.detect_order = 4;           % filter order for detection. 0 to disable the filter. 
par.sort_fmin = 300;            % high pass filter for sorting 
par.sort_fmax = 3000;           % low pass filter for sorting (default 3000) 
par.sort_order = 2;             % filter order for sorting. 0 to disable the filter. 
par.ref_ms = 10;                % detector dead time, minimum refractory period (in ms) 
par.detection = 'pos';          % type of threshold ('pos','neg','both') 
 
% INTERPOLATION PARAMETERS 
par.int_factor = 5;             % interpolation factor 
par.interpolation = 'y';        % interpolation with cubic splines ('y','n') (default) 
 
% FEATURES PARAMETERS 
par.min_inputs = 10;            % number of inputs to the clustering 
par.max_inputs = 0.75;          % number of inputs to the clustering. if < 1 it will the that proportion of the maximum. 
par.scales = 4;                 % number of scales for the wavelet decomposition 
par.features = 'wav';           % type of feature ('wav' or 'pca') 

% FORCE MEMBERSHIP PARAMETERS 
par.template_sdnum = 3;         % max radius of cluster in std devs. 
par.template_k = 10;            % # of nearest neighbors 
par.template_k_min = 10;        % min # of nn for vote 
par.template_type = 'center';   % 'nn', 'center', 'ml', 'mahal' 
par.force_feature = 'spk';      % feature use for forcing (whole spike shape) 
%par.force_feature = 'wav';     % feature use for forcing (wavelet coefficients). 
par.force_auto = false;          % automatically force membership (only for batch scripts). 
 
% TEMPLATE MATCHING 
par.match = 'y';                % for template matching 
%par.match = 'n';               % for no template matching 
par.max_spk = 40000;            % max. # of spikes before starting templ. match. 
par.permut = 'y';               % for selection of random 'par.max_spk' spikes before starting templ. match. 
% par.permut = 'n';             % for selection of the first 'par.max_spk' spikes before starting templ. match. 
 
% HISTOGRAM PARAMETERS 
par.nbins = 100;                % # of bins for the ISI histograms 
par.bin_step = 1;               % percentage number of bins to plot 
 