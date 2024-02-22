function results =  correlation_analysis(app)


if ~isempty(app.burst_ints)
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = nan(length(tmp)+1,2);
    int_names = {'full interval'};

    for i = 1: length(tmp)        
        borders(i+1,:) = app.burst_ints(tmp(i)).borders;
        
        int_names{i+1} = app.burst_ints(tmp(i)).name;
    end
else
    borders = [nan, nan];
    int_names = {'full interval'};
end


lags = [app.edt_corre_lag1.Value, app.edt_corre_lag2.Value, app.edt_corre_lag3.Value];
ds = app.edt_downsampling.Value;

ch1_idx = find(strcmp(app.popup_coorelation_signal1.Value, app.popup_coorelation_signal1.Items));
[data_1,ts_1,name_1, unit_1] = current_signal(app, ch1_idx);
data_1(:,2) = ts_1(1):ts_1(1):ts_1(2);
name_1 = char(string(name_1));

ch2_idx = find(strcmp(app.popup_coorelation_signal2.Value, app.popup_coorelation_signal2.Items));
[data_2,ts_2,name_2, unit_2] = current_signal(app, ch2_idx);
data_2(:,2) = ts_2(1):ts_2(1):ts_2(2);
name_2 = char(string(name_2));

if ts_1(1) > ts_2(1)
    idx = find(data_2(:,2) == data_1(1,2));
    data_2(1:idx-1,:) = [];
else
    idx = find(data_1(:,2) == data_2(1,2));
    data_1(1:idx-1,:) = [];
end

data_1_ds = [downsample(data_1(:,1), .01/ts_1(1)),downsample(data_1(:,2), .01/ts_1(1))];

data_2_ds = [downsample(data_2(:,1), .01/ts_2(1)),downsample(data_2(:,2), .01/ts_2(1))];



data_1_dsds = [downsample(data_1_ds(:,1), ds), downsample(data_1_ds(:,2), ds)];
data_2_dsds = [downsample(data_2_ds(:,1), ds), downsample(data_2_ds(:,2), ds)];

borders(borders(:,1) < 1 ,1) = 1;
borders(borders(:,2) > floor(min([data_1_dsds(end,2), data_2_dsds(end,2)])),2) = floor(min([data_1_dsds(end,2), data_2_dsds(end,2)]));

if app.chkbx_movmean.Value
    data_1_dsds(:,1) = movmean(data_1_dsds(:,1),app.edt_movmean.Value /(.01*ds));
    data_2_dsds(:,1) = movmean(data_2_dsds(:,1),app.edt_movmean.Value /(.01*ds));
end


seriesname={name_1, name_2};

results = struct('lag1', struct('ax11',struct('Data',[],'Values',[],'NumBins',[],'BinEdges',[],'BinWidth',[],'BinLimits',[],'Normalization',[],'FaceColor',[],'EdgeColor',[],'xlim',[],'Ylim',[]),...
                                'ax12',struct('text',struct('String',[],'Position',[]),'line1',struct('XData',[],'YData',[]),'line2',struct('XData',[],'YData',[]),'XLim',[],'YLim',[]),...
                                'ax22',struct('Data',[],'Values',[],'NumBins',[],'BinEdges',[],'BinWidth',[],'BinLimits',[],'Normalization',[],'FaceColor',[],'EdgeColor',[],'xlim',[],'Ylim',[]), ...
                                'ax21',struct('text',struct('String',[],'Position',[]),'line1',struct('XData',[],'YData',[]),'line2',struct('XData',[],'YData',[]),'XLim',[],'YLim',[]),...
                                'r',[],'p',[]),...
                 'lag2', struct('ax11',struct('Data',[],'Values',[],'NumBins',[],'BinEdges',[],'BinWidth',[],'BinLimits',[],'Normalization',[],'FaceColor',[],'EdgeColor',[],'xlim',[],'Ylim',[]),...
                                'ax12',struct('text',struct('String',[],'Position',[]),'line1',struct('XData',[],'YData',[]),'line2',struct('XData',[],'YData',[]),'XLim',[],'YLim',[]),...
                                'ax22',struct('Data',[],'Values',[],'NumBins',[],'BinEdges',[],'BinWidth',[],'BinLimits',[],'Normalization',[],'FaceColor',[],'EdgeColor',[],'xlim',[],'Ylim',[]), ...
                                'ax21',struct('text',struct('String',[],'Position',[]),'line1',struct('XData',[],'yData',[]),'line2',struct('XData',[],'YData',[]),'XLim',[],'YLim',[]),...
                                'r',[],'p',[]),...
                 'lag3', struct('ax11',struct('Data',[],'Values',[],'NumBins',[],'BinEdges',[],'BinWidth',[],'BinLimits',[],'Normalization',[],'FaceColor',[],'EdgeColor',[],'xlim',[],'Ylim',[]),...
                                'ax12',struct('text',struct('String',[],'Position',[]),'line1',struct('XData',[],'YData',[]),'line2',struct('XData',[],'YData',[]),'XLim',[],'YLim',[]),...
                                'ax22',struct('Data',[],'Values',[],'NumBins',[],'BinEdges',[],'BinWidth',[],'BinLimits',[],'Normalization',[],'FaceColor',[],'EdgeColor',[],'xlim',[],'Ylim',[]), ...
                                'ax21',struct('text',struct('String',[],'Position',[]),'line1',struct('XData',[],'YData',[]),'line2',struct('XData',[],'YData',[]),'XLim',[],'YLim',[]),...
                                'r',[],'p',[]));
% results = struct('ax1', [],'ax2', [],'ax3', [],'ax4', []);
for i = 1: size(borders,1)
    
    if i == 1
        data_1_dsds_int = data_1_dsds;
        data_2_dsds_int = data_2_dsds;
    else
        data_1_dsds_int = data_1_dsds(int32(borders(i,1)*1/(.01*ds)):int32(borders(i,2)*1/(.01*ds)),:);
        data_2_dsds_int = data_2_dsds(int32(borders(i,1)*1/(.01*ds)):int32(borders(i,2)*1/(.01*ds)),:);
    end

    lag01=lags(1)*1/(.01*ds);
  

    clear h1
    h1 = figure('Visible','off');
    ca1 = gca;
    [r,p] = corrplot1(ca1,[data_1_dsds_int(1:end-lag01,1),data_2_dsds_int(lag01+1:end,1)],'varNames',{name_1,name_2} ,'type','Pearson','testR','on','alpha',0.05);
    
    results(i).lag1.r = r;
    results(i).lag1.p = p;
    results(i).lag1.ax11.Data = h1.Children(1).Children.Data;
    results(i).lag1.ax11.BinEdges = h1.Children(1).Children.BinEdges;
    results(i).lag1.ax11.XLim = h1.Children(1).Children.Parent.XLim;
    results(i).lag1.ax11.YLim = h1.Children(1).Children.Parent.YLim;

    results(i).lag1.ax12.text.String = h1.Children(4).Children(1).String;
    results(i).lag1.ax12.text.Position = h1.Children(4).Children(1).Position;
    results(i).lag1.ax12.line1.XData = h1.Children(4).Children(2).XData;
    results(i).lag1.ax12.line1.YData = h1.Children(4).Children(2).YData;
    results(i).lag1.ax12.line2.XData = h1.Children(4).Children(3).XData;
    results(i).lag1.ax12.line2.YData = h1.Children(4).Children(3).YData;
    results(i).lag1.ax12.XLim = h1.Children(4).Children(3).Parent.XLim;
    results(i).lag1.ax12.YLim = h1.Children(4).Children(3).Parent.YLim;

    results(i).lag1.ax21.text.String = h1.Children(5).Children(1).String;
    results(i).lag1.ax21.text.Position = h1.Children(5).Children(1).Position;
    results(i).lag1.ax21.line1.XData = h1.Children(5).Children(2).XData;
    results(i).lag1.ax21.line1.YData = h1.Children(5).Children(2).YData;
    results(i).lag1.ax21.line2.XData = h1.Children(5).Children(3).XData;
    results(i).lag1.ax21.line2.YData = h1.Children(5).Children(3).YData;
    results(i).lag1.ax21.XLim = h1.Children(5).Children(3).Parent.XLim;
    results(i).lag1.ax21.YLim = h1.Children(5).Children(3).Parent.YLim;

    results(i).lag1.ax22.Data = h1.Children(2).Children.Data;
    results(i).lag1.ax22.BinEdges = h1.Children(2).Children.BinEdges;
    results(i).lag1.ax22.XLim = h1.Children(2).Children.Parent.XLim;
    results(i).lag1.ax22.YLim = h1.Children(2).Children.Parent.YLim;

    close(h1)

     lag02=lags(2)*1/(.01*ds);
  

    clear h2
    h2 = figure('Visible','on');
    ca2 = gca;
    [r,p] = corrplot1(ca2,[data_1_dsds_int(1:end-lag02,1),data_2_dsds_int(lag02+1:end,1)],'varNames',{name_1,name_2} ,'type','Pearson','testR','on','alpha',0.05);
    
 
    results(i).lag2.r = r;
    results(i).lag2.p = p;
    results(i).lag2.ax11.Data = h2.Children(1).Children.Data;
    results(i).lag2.ax11.BinEdges = h2.Children(1).Children.BinEdges;
    results(i).lag2.ax11.XLim = h2.Children(1).Children.Parent.XLim;
    results(i).lag2.ax11.YLim = h2.Children(1).Children.Parent.YLim;

    results(i).lag2.ax12.text.String = h2.Children(4).Children(1).String;
    results(i).lag2.ax12.text.Position = h2.Children(4).Children(1).Position;
    results(i).lag2.ax12.line1.XData = h2.Children(4).Children(2).XData;
    results(i).lag2.ax12.line1.YData = h2.Children(4).Children(2).YData;
    results(i).lag2.ax12.line2.XData = h2.Children(4).Children(3).XData;
    results(i).lag2.ax12.line2.YData = h2.Children(4).Children(3).YData;
    results(i).lag2.ax12.XLim = h2.Children(4).Children(3).Parent.XLim;
    results(i).lag2.ax12.YLim = h2.Children(4).Children(3).Parent.YLim;

    results(i).lag2.ax21.text.String = h2.Children(5).Children(1).String;
    results(i).lag2.ax21.text.Position = h2.Children(5).Children(1).Position;
    results(i).lag2.ax21.line1.XData = h2.Children(5).Children(2).XData;
    results(i).lag2.ax21.line1.YData = h2.Children(5).Children(2).YData;
    results(i).lag2.ax21.line2.XData = h2.Children(5).Children(3).XData;
    results(i).lag2.ax21.line2.YData = h2.Children(5).Children(3).YData;
    results(i).lag2.ax21.XLim = h2.Children(5).Children(3).Parent.XLim;
    results(i).lag2.ax21.YLim = h2.Children(5).Children(3).Parent.YLim;

    results(i).lag2.ax22.Data = h2.Children(2).Children.Data;
    results(i).lag2.ax22.BinEdges = h2.Children(2).Children.BinEdges;
    results(i).lag2.ax22.XLim = h2.Children(2).Children.Parent.XLim;
    results(i).lag2.ax22.YLim = h2.Children(2).Children.Parent.YLim;

    close(h2)

    lag03=lags(3)*1/(.01*ds);
  
    clear h3
    h3 = figure('Visible','on');
    ca3 = gca;
    [r,p] = corrplot1(ca3,[data_1_dsds_int(1:end-lag03,1),data_2_dsds_int(lag03+1:end,1)],'varNames',{name_1,name_2} ,'type','Pearson','testR','on','alpha',0.05);
    
    results(i).lag3.r = r;
    results(i).lag3.p = p;
    results(i).lag3.ax11.Data = h3.Children(1).Children.Data;
    results(i).lag3.ax11.BinEdges = h3.Children(1).Children.BinEdges;
    results(i).lag3.ax11.XLim = h3.Children(1).Children.Parent.XLim;
    results(i).lag3.ax11.YLim = h3.Children(1).Children.Parent.YLim;

    results(i).lag3.ax12.text.String = h3.Children(4).Children(1).String;
    results(i).lag3.ax12.text.Position = h3.Children(4).Children(1).Position;
    results(i).lag3.ax12.line1.XData = h3.Children(4).Children(2).XData;
    results(i).lag3.ax12.line1.YData = h3.Children(4).Children(2).YData;
    results(i).lag3.ax12.line2.XData = h3.Children(4).Children(3).XData;
    results(i).lag3.ax12.line2.YData = h3.Children(4).Children(3).YData;
    results(i).lag3.ax12.XLim = h3.Children(4).Children(3).Parent.XLim;
    results(i).lag3.ax12.YLim = h3.Children(4).Children(3).Parent.YLim;

    results(i).lag3.ax21.text.String = h3.Children(5).Children(1).String;
    results(i).lag3.ax21.text.Position = h3.Children(5).Children(1).Position;
    results(i).lag3.ax21.line1.XData = h3.Children(5).Children(2).XData;
    results(i).lag3.ax21.line1.YData = h3.Children(5).Children(2).YData;
    results(i).lag3.ax21.line2.XData = h3.Children(5).Children(3).XData;
    results(i).lag3.ax21.line2.YData = h3.Children(5).Children(3).YData;
    results(i).lag3.ax21.XLim = h3.Children(5).Children(3).Parent.XLim;
    results(i).lag3.ax21.YLim = h3.Children(5).Children(3).Parent.YLim;

    results(i).lag3.ax22.Data = h3.Children(2).Children.Data;
    results(i).lag3.ax22.BinEdges = h3.Children(2).Children.BinEdges;
    results(i).lag3.ax22.XLim = h3.Children(2).Children.Parent.XLim;
    results(i).lag3.ax22.YLim = h3.Children(2).Children.Parent.YLim;

    close(h3)

    % .Toolbar = []
    % .PlotBoxAspectRatio = [1 0.7895 0.7895]
    
end

end
