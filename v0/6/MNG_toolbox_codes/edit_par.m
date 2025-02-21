function par = edit_par(par_in)

fields = { 'max_spikes_plot',  'to_plot_std', 'all_classes_ax',...
          'mintemp', 'maxtemp', 'tempstep',... 
          'SWCycles', 'KNearNeighb', 'min_clus',...
          'randomseed', 'temp_plot', 'c_ov',...
          'elbow_min', 'int_factor', 'interpolation',...
          'min_inputs', 'max_inputs', 'scales',...
          'features', 'template_sdnum', 'template_k',... 
          'template_k_min', 'template_type', 'force_feature',...
          'match', 'max_spk', 'permut',...
          'nbins', 'bin_step'};
value = {1000,  1, 'mean',...
         0.00, 0.251, 0.01,...
         100, 11, 20,...
         0, 'log', 0.7,...
         0.4, 5, 'y',...
         10, 0.75, 4,... 
         'wav', 3, 10,...
         10,'center','spk',...
         'y', 40000, 'y',... 
         100, 1};
type = {'numeric', 'numeric', 'text',... 
        'numeric', 'numeric', 'numeric',...
        'numeric', 'numeric', 'numeric',...
        'numeric', 'text', 'numeric',...
        'numeric', 'numeric', 'text',...
        'numeric', 'numeric', 'numeric',...
        'text', 'numeric', 'numeric',...
        'numeric', 'text', 'text',...
        'text', 'numeric', 'text',...
        'numeric', 'numeric'};
lbls = {'Max # spikes to plot:', '# of std from mean to plot:', 'Spikes plot:',...
        'Minimal temperatre:', 'Maximal temperature:', 'Temperatur steps:',...
        'Iterations for each temp.:', '# nearest neighbors for SPC:', 'Minimal clustersize:',...
        'Randomseed:', 'Scale for temp.:', 'Overlapping coefficient for inclusion:',...
        'Param. for regime border detection:', 'Interpolation factor:',  'Cubic splines for interpol.:'...
        '# inputs to clustering:', 'Proportion of inputs', '# scales for wavelet recomp.:',...
        'Type of features:', 'Max. radius cluster (# std)', '# nearest neighbors:',...
        'Min # nn for vote', 'Template type:', 'Feature for forcing:'...
        'Template matching:', '# spikes to templ.match.:', 'Random spikes for templ.match:',...
        '# of bins for histograms:', 'Proportion of bins to plot:'};
defs = {'1000', '1', '"mean" ("mean" / "all")',...
        '0.00', '0.251', '0.01', '100', '11', '20',...
        '0 (0 = clock value)', '"log" ("lin" / "log")', '0.7',...
        '0.4', '5', '"y" ("y" / "n")',...
        '10', '0.75 (0 - 1)', '4',...
        '"wav" ("wav" / "pca")', '3', '10',...
        '10', '"center" ("nn" / "center" / "ml" / "mahal")', '"spk" ("spk" / "wav")',...
        '"y" ("y" / "n")', '40000', '"y" ("y" / "n")',...
        '100','1'};

pos = {[20, 424,300,22],[320, 424,100,22],[420, 424,80,22]};
     
       
ypos = [425, 401, 387, 339, 315, 291, 267, 243, 219, 195, 171, 147, 123, 75, 51,...
        425, 401, 387, 363, 315, 291, 267, 243, 219, 171, 147, 123, 75, 51];
xpos = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,...
        500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500];

par = struct([]);

for i = 1 : length(fields) 

    if isfield(par_in,fields{i})
        eval(['par.' fields{i} '= par_in.' fields{i} ';'])
    else
        eval([ 'par(1).' fields{i} ' = value{i};']) %= value{i};              
    end
end

h = uifigure('Position',[100,100,1000,500]) ;

for i = 1 : length(fields)
    eval(['lbl_' fields{i} '_n = uilabel(h,"Text",' lbls{i}  ',"Position", [' num2str(xpos(i)) ',' num2str(ypos(i)) ',300,22 ]);'])
    eval(['edt_' fields{i} '= uieditfield(h,' type{i} ',"Position", [' num2str(xpos(i)+300) ',' num2str(ypos(i)-1) ',100,24 ],"Value", ' value{i} ');'])
    eval(['lbl_' fields{i} '_d = uilabel(h,"Text",' defs{i}  ',"Position", [' num2str(xpos(i)+400) ',' num2str(ypos(i)) ',80,22 ]);'])
    eval(['lbl_' fields{i} '_d.Value']) = defs{i};
end




