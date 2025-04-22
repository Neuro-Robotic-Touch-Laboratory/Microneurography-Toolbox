ridx = nan(size(app.hb_res.t_events));
fs = 10000;
tmp_data = [x,zeros(1,3*fs)];
dur = 3*fs-1;

summ = zeros(fs*3,1);
for i = 1:length(app.hb_res.t_events)
    ridx(i) = find(t_msna>=app.hb_res.t_events(i),1);
    summ= summ+tmp_data(ridx(i):ridx(i)+dur)';
end
[b,a] = butter(4,5/5000,"low");
sumf = filtfilt(b,a,summ);
[~,mxi] = max(movmean(sumf,500));
[~,lcs] = findpeaks(-sumf-mean(-sumf),'MinPeakHeight',max(-sumf-mean(-sumf))/10);%findpeaks(-sumf,'MinPeakHeight',min(sumf)/10,"NPeaks",2);
[~,idx] = sort(abs(lcs-mxi),"ascend");
idx(3:end) = [];


sumts = (1:length(summ))/10000;
figure
plot (sumts,sumf)
hold on 
plot (sumts(mxi),sumf(mxi),'*r')
%plot (sumts(lcs),sumf(lcs),'*g')
plot (sumts(lcs(idx)),sumf(lcs(idx)),'*g')

pk_del = sumts(mxi);
% win = lcs;
win = sort(lcs(idx),'ascend');
% win = [win(1)-1000, win(2)+1000];
[b,a] = butter(2,[0.05,4]/5000,"bandpass");

x_f = filtfilt(b,a,tmp_data);
dx_f = diff(x_f).*circshift(diff(x_f),1);


brsts= nan(length(ridx),7);

for i = 1:length(ridx)
    [~,tmp_pk] = max(x_f(ridx(i)+win(1):ridx(i)+win(2)));
    if (tmp_pk ~= 1) & (tmp_pk ~= (abs(diff(win))+1))
        tmp_idx = find(dx_f(ridx(i)+win(1):ridx(i)+win(2))<0);
        strt = find(tmp_idx<tmp_pk,1,'last');
        stp =  find(tmp_idx>tmp_pk,1,'first');
        strt= tmp_idx(strt);
        if isempty(strt)
            strt =1;
        end
        stp = tmp_idx(stp);
        if isempty(stp)
            stp = win(2)-win(1)+1;
        end
        brsts(i,1:3) = [strt-1+ridx(i)+win(1), tmp_pk-1+ridx(i)+win(1), stp-1+ridx(i)+win(1)];
        brsts(i,4) =  max(tmp_data(brsts(i,2))-tmp_data(brsts(i,1)),tmp_data(brsts(i,2))-tmp_data(brsts(i,3)));
        brsts(i,5) = mean(diff(tmp_data(brsts(i,1):brsts(i,2))));
        brsts(i,6) = mean(diff(tmp_data(brsts(i,2):brsts(i,3))));
        brsts(i,7) = (brsts(i,2)-ridx(i)+1)/fs;

    end
end

brsts(isnan(brsts(:,3)),:) = [];

brsts((brsts(:,1)>= length(x)),:) = [];

brsts((brsts(:,2)>= length(x)),:) = [];

brsts((brsts(:,3)>= length(x)),:) = [];


%% threshold calc 1 
bins= min(brsts(:,7)):0.05:  max(brsts(:,7)+0.05);
pcntls =zeros(length(bins)-1,1);
for i = 1 :length(bins)-1
binvals = brsts(brsts(:,7)>bins(i)& brsts(:,7)<bins(i+1),4);
binvals = sort(binvals,"ascend");
len = length(binvals);
pcntls(i) = binvals(ceil(len*0.6));%[binvals(ceil(len*0.5)),binvals(ceil(len*0.6)),binvals(ceil(len*0.7)),binvals(ceil(len*0.8)),binvals(ceil(len*0.9))];
end
%bla = mean(pcntls,1);

%% threshold calc 2

%%
% Define your data points
xp = brsts(:,4); % Independent variable
yp = brsts(:,7); % Dependent variable

% Fit a linear model
coefficients = polyfit(xp, yp, 1); % 1st-degree polynomial (linear regression)

% Generate fitted values
x_fit = linspace(min(xp), max(xp), 100); % Generate more points for smooth plotting
y_fit = polyval(coefficients, x_fit);

% Compute distances from the regression line
y_pred = polyval(coefficients, xp);
distances = abs(yp - y_pred);

dist = 0.1; % Define distance threshold
selected_points = (distances <= dist);
selected_x = xp(selected_points);
selected_y = yp(selected_points);

% Compute selection area borders
y_upper = y_fit + dist;
y_lower = y_fit - dist;

% Plot original data and fitted line
figure;
scatter(xp, yp, 'bo', 'filled'); % Plot data points
hold on;

scatter(selected_x, selected_y, 'go', 'filled'); % Highlight selected points
plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Plot regression line
plot(x_fit, y_upper, 'k--', 'LineWidth', 1); % Upper boundary
plot(x_fit, y_lower, 'k--', 'LineWidth', 1); % Lower boundary
hold off;
title('Linear Regression Fit with Selection Area');
xlabel('X');
ylabel('Y');
grid on;
legend('Data Points', 'Fitted Line', 'Selected Points', 'Selection Borders');

tmp = xp(~selected_points);
mn = mean (tmp);
st =std(tmp);
%%

thrsh = std(pcntls);
thrsh = mn+st;

brsts(:,8)= brsts(:,4) >= thrsh;

st = std([brsts(:,5);-brsts(:,6)]);
mn = mean([brsts(:,5);-brsts(:,6)]);
brsts_idx = find (brsts(:,8) == 1);
rem_idx_up = find(brsts(brsts_idx,5)> mn+4.5*st);
rem_idx_un = find(-brsts(brsts_idx,6)> mn+4.5*st); 
rem_idx_lp = find(brsts(brsts_idx,5)< mn-st);
rem_idx_ln = find(-brsts(brsts_idx,6)< mn-st); 

rem_idx = unique([rem_idx_up;rem_idx_un;rem_idx_lp;rem_idx_ln]);
brsts(brsts_idx(rem_idx),8) = 0;

mn = mean(brsts(:,7));
st =std(brsts(:,7));
brsts_idx = find (brsts(:,8) == 1);
rem_idx_p = find(brsts(brsts_idx,7)> mn+2.5*st);
rem_idx_n = find(brsts(brsts_idx,7)< mn-2.5*st);

rem_idx = unique([rem_idx_p;rem_idx_n]);
brsts(brsts_idx(rem_idx),8) = 0;
%%
%intrbins =(bins(1:end-1)+bins(2:end))./2;
figure 
plot(brsts(logical(brsts(:,8)),4),brsts(logical(brsts(:,8)),7),'ro');
hold on
plot(brsts(~logical(brsts(:,8)),4),brsts(~logical(brsts(:,8)),7),'ko');
xlabel amplitude
ylabel latency

line([thrsh thrsh],[bins(1) bins(end)])

figure 
plot(brsts(logical(brsts(:,8)),5),brsts(logical(brsts(:,8)),6),'ro')
hold on
plot(brsts(~logical(brsts(:,8)),5),brsts(~logical(brsts(:,8)),6),'ko')
xlabel 'rise slope'
ylabel 'fall slope'

figure
tmp1 = diff([x(int64(brsts(:,1))); x(int64(brsts(:,2)))],1,1);
tmp2 = abs(diff([x(int64(brsts(:,2))); x(int64(brsts(:,3)))],1,1));
plot(tmp1(logical(brsts(:,8))),tmp2(logical(brsts(:,8))),'ro')
hold on
plot(tmp1(~logical(brsts(:,8))),tmp2(~logical(brsts(:,8))),'ko')
xlabel 'rise amplitude'
ylabel 'fall amplitude'





figure

hb_idx = app.hb_res.use_beats(:,1) ;

plot (t_msna,x)
hold on 

hb = nan(length(app.hb_res.t_events(hb_idx))*3,2);
tmp_idx = 1:3:length(app.hb_res.t_events(hb_idx))*3;
hb(tmp_idx,1)= app.hb_res.t_events(hb_idx);
hb(tmp_idx+1,1)= app.hb_res.t_events(hb_idx);
yl = ylim(app.ax_msna_int);
hb(tmp_idx,2)= yl(1);
hb(tmp_idx+1,2)= yl(2);
plot (hb(:,1),hb(:,2),'LineWidth',1.5,'Color',[0 0.5 0.2],'HitTest','off')

        
colors = [1.0000, 0, 0;...
          0, 1.0000, 0;...
          0, 0, 0.1724;...
          1.0000, 0.1034, 0.7241;...
          1.0000, 0.8276, 0];

burst_plot1 =[nan,nan]; burst_plot2 =[nan,nan]; burst_plot3 =[nan,nan]; burst_plot4 =[nan,nan]; burst_plot5 =[nan,nan]; burst_plot_rem = [nan,nan];
burst_amp_plot1 = [nan,nan]; burst_amp_plot2 = [nan,nan]; burst_amp_plot3 = [nan,nan]; burst_amp_plot4 = [nan,nan]; burst_amp_plot5 = [nan,nan]; burst_amp_plot_rem = [nan,nan]; 
tmp = find(brsts(:,8) == 1);
tmp_brsts = brsts(tmp,:);

for j = 1:5: ceil(size(tmp_brsts,1)/5)*5 %%%% floor
    if j <= size(tmp_brsts,1)
        burst_plot1 = [burst_plot1;[(t_msna(tmp_brsts(j,1):tmp_brsts(j,3)))',...
                       (x(tmp_brsts(j,1):tmp_brsts(j,3)))'];...
                       [t_msna(tmp_brsts(j,1)), x(tmp_brsts(j,1))];...
                       [nan,nan]];
    end
    if j+1 <= size(tmp_brsts,1)  
        burst_plot2 = [burst_plot2;[(t_msna(tmp_brsts(j+1,1):tmp_brsts(j+1,3)))',...
                       (x(tmp_brsts(j+1,1):tmp_brsts(j+1,3)))'];...
                       [t_msna(tmp_brsts(j+1,1)), x(tmp_brsts(j+1,1))];...
                       [nan,nan]];
    end
    if j+2 <= size(tmp_brsts,1)
        burst_plot3 = [burst_plot3;[(t_msna(tmp_brsts(j+2,1):tmp_brsts(j+2,3)))',...
                       (x(tmp_brsts(j+2,1):tmp_brsts(j+2,3)))'];...
                       [t_msna(tmp_brsts(j+2,1)), x(tmp_brsts(j+2,1))];...
                       [nan,nan]];
    end
    if j+3 <= size(tmp_brsts,1)
        burst_plot4 = [burst_plot4;[(t_msna(tmp_brsts(j+3,1):tmp_brsts(j+3,3)))',...
                       (x(tmp_brsts(j+3,1):tmp_brsts(j+3,3)))'];...
                       [t_msna(tmp_brsts(j+3,1)), x(tmp_brsts(j+3,1))];...
                       [nan,nan]];
    end
    if j+4 <= size(tmp_brsts,1) 
        burst_plot5 = [burst_plot5;[(t_msna(tmp_brsts(j+4,1):tmp_brsts(j+4,3)))',...
                       (x(tmp_brsts(j+4,1):tmp_brsts(j+4,3)))'];...
                       [t_msna(tmp_brsts(j+4,1)), x(tmp_brsts(j+4,1))];...
                       [nan,nan]];
    end
end
tmp = find(brsts(:,8)==0);
tmp_brsts = brsts(tmp,:);

for j = 1: length(tmp)
    burst_plot_rem = [burst_plot_rem;[(t_msna(tmp_brsts(j,1):tmp_brsts(j,3)))',...
                       (x(tmp_brsts(j,1):tmp_brsts(j,3)))'];...
                       [t_msna(tmp_brsts(j,1)), x(tmp_brsts(j,1))];...
                       [nan,nan]];
end
hold on 
plot(burst_plot1(:,1),burst_plot1(:,2),'Color',colors(1,:),'LineWidth',1.5)
plot(burst_plot2(:,1),burst_plot2(:,2),'Color',colors(2,:),'LineWidth',1.5)
plot(burst_plot3(:,1),burst_plot3(:,2),'Color',colors(3,:),'LineWidth',1.5)
plot(burst_plot4(:,1),burst_plot4(:,2),'Color',colors(4,:),'LineWidth',1.5)
plot(burst_plot5(:,1),burst_plot5(:,2),'Color',colors(5,:),'LineWidth',1.5)
plot(burst_plot_rem(:,1),burst_plot_rem(:,2),'Color',[.7,.7,.7],'LineWidth',1.5)


%% calc burst measures

for i = 1 : size(brsts,1)
    tmp = [t_msna(brsts(i,1):brsts(i,3))',x(brsts(i,1):brsts(i,3))'*1000];
    tmp(:,2) = tmp(:,2) - min (tmp(:,2));

    brsts(i,9)  = max(tmp(:,2)); 
    brsts(i,10) = trapz(tmp(:,1),tmp(:,2));
end

data = struct();

tmp = nan(1,length(app.data));
idx=0;
for i = 1: length(app.data)
    tmp(i) =  app.data(i).tic_multipl;
    if strcmp (app.data(i).name,'BURST INTEGRAL')
        idx = i;
    end
end

min_ts = app.data(find(tmp ==1,1)).ts(1);


tmp = 0.01:0.01:diff(app.settings.interval(1,:));
tmp_brsts = brsts(brsts(:,8)==1,:);



data(1).name = 'BURST INTEGRAL';
data(1).unit = 'µV*ms';
data(1).tic_multipl = 0.01/min_ts;
data(1).data = interp1(t_msna(tmp_brsts(:,2)),tmp_brsts(:,10),tmp)';
tmp_idx = [find(~isnan(data(1).data),1) ,  find(~isnan(data(1).data),1,'last')];
data(1).data(1:tmp_idx(1)-1) = data(1).data(tmp_idx(1));
data(1).data(tmp_idx(2)+1:end) = data(1).data(tmp_idx(2));
data(1).ts = [tmp(1),tmp(end)];
data(1).derived = true;

data(2).name = 'BURST DURATION';
data(2).unit = 's';
data(2).tic_multipl = 0.01/min_ts;
data(2).data = interp1(t_msna(tmp_brsts(:,2)),(tmp_brsts(:,3)-tmp_brsts(:,1))/fs,tmp)';
tmp_idx = [find(~isnan(data(2).data),1) ,  find(~isnan(data(2).data),1,'last')];
data(2).data(1:tmp_idx(1)-1) = data(2).data(tmp_idx(1));
data(2).data(tmp_idx(2)+1:end) = data(2).data(tmp_idx(2));
data(2).ts = [tmp(1),tmp(end)];
data(2).derived = true;

data(3).name = 'BURST amplitude';
data(3).unit = 's';
data(3).tic_multipl = 0.01/min_ts;
data(3).data = interp1(t_msna(tmp_brsts(:,2)),tmp_brsts(:,9),tmp)';
tmp_idx = [find(~isnan(data(3).data),1) ,  find(~isnan(data(3).data),1,'last')];
data(3).data(1:tmp_idx(1)-1) = data(3).data(tmp_idx(1));
data(3).data(tmp_idx(2)+1:end) = data(3).data(tmp_idx(2));
data(3).ts = [tmp(1),tmp(end)];
data(3).derived = true;

data(4).name = 'BURST INTERBURST INTERVAL';
data(4).unit = 's';
data(4).tic_multipl = 0.01/min_ts;
data(4).data = interp1(t_msna(tmp_brsts(1:end-1,2)),diff(tmp_brsts(:,2)),tmp)';
tmp_idx = [find(~isnan(data(4).data),1) ,  find(~isnan(data(4).data),1,'last')];
data(4).data(1:tmp_idx(1)-1) = data(4).data(tmp_idx(1));
data(4).data(tmp_idx(2)+1:end) = data(4).data(tmp_idx(2));
data(4).ts = [tmp(1),tmp(end)];
data(4).derived = true;

data(5).name = 'BURST delay';
data(5).unit = 's';
data(5).tic_multipl = 0.01/min_ts;
data(5).data = interp1(t_msna(tmp_brsts(:,2)),tmp_brsts(:,7),tmp)';
tmp_idx = [find(~isnan(data(3).data),1) ,  find(~isnan(data(3).data),1,'last')];
data(5).data(1:tmp_idx(1)-1) = data(3).data(tmp_idx(1));
data(5).data(tmp_idx(2)+1:end) = data(3).data(tmp_idx(2));
data(5).ts = [tmp(1),tmp(end)];
data(5).derived = true;


bla = true;