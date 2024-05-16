function covar_out = covardlg(covar_in,use_rec,anocovas)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

uif = uifigure('Position',[100,100,800,600],'Color',[0.80,0.95,0.90]);
% uif.CloseRequestFcn = @(src,event)my_closereq(src);
col_ed = [false,false];
col_name = {'record', 'interval'};

for i = 1: length(covar_in)
    col_ed(i+2) = true;
    col_name{i+2} = ['Covars gr. ' num2str(anocovas(i))];
end
data = {'Covariate', 'name:'};
for i = 1:length(use_rec)
    data{i+1,1} = use_rec(i).recording;
    data{i+1,2} = use_rec(i).interval;
end
for i =1: length(covar_in)
    for j = 1: length(covar_in(i).checkarray)
        if ~covar_in(i).checkarray(j)
            data{1+j,2+i} = nan;
        else
            data{1+j,2+i} = [];
        end
    end
end

% covar_out = [];
uit = uitable(uif,'Position',[20,50,760,500],'ColumnEditable',col_ed, 'ColumnName', col_name,'RowName', [],'Data',data);
uit.CellEditCallback = @(src,event) cell_value_changed(src,event,covar_in);
uib = uibutton(uif, 'Text', 'Save & continue','ButtonPushedFcn', @(src,event) close_return(uif), 'Position',[300,10,200,30]);

uiwait(uif)
data = uit.Data;
covar_out = covar_in;
data(:,1:2) =[];
for i= 1: length(covar_in)
    covar_out(i).name = data{1,i};
    covar_out(i).values = vertcat(data{2:end,i});
end

close(uif)
end

function cell_value_changed(src,event,covar_in)
if event.Indices(1) == 1
    if ~isnan(covar_in(event.Indices(2)-2).asprev)
        if isequal(event.NewData, src.Data{1,covar_in(event.Indices(2)-2).asprev+2})
            src.Data(:,event.Indices(2)) = src.Data(:,covar_in(event.Indices(2)-2).asprev+2);
        end
    end
else
    [data,flag] = str2num(event.NewData);
    if flag
        if length(data) == 1
            src.Data{event.Indices(1),event.Indices(2)} = data;
        else
            src.Data{event.Indices(1),event.Indices(2)} = event.PreviousData;
        end
    else
        src.Data{event.Indices(1),event.Indices(2)} = event.PreviousData;
    end
end
end

function covar_out = close_return(uif)
data = uif.Children(2).Data;
data(:,1:2) = [];
done_flag = false;
for i = 1:size(data,1)
    for j = 1:size(data,2)
        done_flag = done_flag | isempty(data{i,j});
    end
end
if done_flag
    warndlg('please fill out all fields')
else
    uiresume(uif)
end

end

