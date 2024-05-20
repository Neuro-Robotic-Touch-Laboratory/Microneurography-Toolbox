function save_trans_results(app)

[int_idxs,~] = listdlg('PromptString',{'Please select intervals ',...
    'to be plotted/saved.',''},...
    'SelectionMode','multiple','ListString',app.popup_int_spect.Items);

[form_idxs,~] = listdlg('PromptString',{'Please select fileformat ',...
    ''},...
    'SelectionMode','multiple','ListString',{'.fig','.jpg','.eps'});
%path = uigetdir;
path = app.settings.output_dir;

tmp(1) = find(app.settings.file_path == '\',1,'last');
tmp(2) = find(app.settings.file_path == '.',1,'last');
file = app.settings.file_path(tmp(1)+1:tmp(2)-1);
file(file == '.') = '-';

yinf=-1;
ysup=1;
x = [-0.5 0.5 0.5 -0.5];
y = [-1 -1 1 1];

for i = 1:length(int_idxs)
    h = figure('Position', get(0, 'Screensize'),'Visible','off');
    set(h, 'NumberTitle', 'off', ...
    'Name', app.popup_int_spect.Items{int_idxs(i)});

    subplot(2,8,[1:3,9:11])
    hold on
    for j = 1 : length(app.spike_res.transduction(int_idxs(i)).cluster)
        plot (app.spike_res.transduction(int_idxs(i)).cluster(j).mean_shape,'LineWidth',3)
        n_sp1{j}=num2str(app.spike_res.transduction(int_idxs(i)).cluster(j).n_spikes);
    end
    
    hold off
    grid on, grid minor
    legend(n_sp1','FontSize',9,'Location','southeast')
    title(['Spikes Sorted Shapes'  app.popup_int_transduction.Items{int_idxs(i)}])

    for j = 1 : length(app.spike_res.transduction(int_idxs(i)).cluster)

        subplot(2,8,8)
        hold on
        plot(app.spike_res.transduction(int_idxs(i)).cluster(j).phase_bp(:,1), app.spike_res.transduction(int_idxs(i)).cluster(j).phase_bp(:,2),'-o','linewidth',1.5)
        ylabel('Pearson First Spike PHASE - BP'), xlabel('cardiac cycle'),title('Pearson First Spike PHASE - BP')
        hold off
        grid on, grid minor
        ylim([yinf ysup])
          
        subplot(2,8,16)
        hold on
        plot(app.spike_res.transduction(int_idxs(i)).cluster(j).phase_rr(:,1), app.spike_res.transduction(int_idxs(i)).cluster(j).phase_rr(:,2),'-o','linewidth',1.5)
        ylabel('Pearson First Spike PHASE - RR'), xlabel('cardiac cycle'),title('Pearson First Spike PHASE - RR')
        hold off
        grid on, grid minor
        ylim([yinf ysup])
        
        subplot(2,8,6:7)
        hold on
        plot(app.spike_res.transduction(int_idxs(i)).cluster(j).latency_bp(:,1), app.spike_res.transduction(int_idxs(i)).cluster(j).latency_bp(:,2),'-o','linewidth',1.5)
        ylabel('Pearson First Spike Lat - BP'), xlabel('cardiac cycle'),title('Pearson First Spike Lat - BP')
        patch(x,y,[0.7 0.7 0.7])
        hold off
        grid on, grid minor
        ylim([yinf ysup])

        subplot(2,8,14:15)
        hold on
        plot(app.spike_res.transduction(int_idxs(i)).cluster(j).latency_rr(:,1), app.spike_res.transduction(int_idxs(i)).cluster(j).latency_rr(:,2),'-o','linewidth',1.5)
        ylabel('Pearson First Spike Lat - RR'), xlabel('cardiac cycle'),title('Pearson First Spike Lat - RR')
        patch(x,y,[0.7 0.7 0.7])
        hold off
        grid on, grid minor
        ylim([yinf ysup])

        subplot(2,8,4:5)
        hold on
        plot(app.spike_res.transduction(int_idxs(i)).cluster(j).fr_bp(:,1), app.spike_res.transduction(int_idxs(i)).cluster(j).fr_bp(:,2),'-o','linewidth',1.5),
        ylabel('Pearson FR - BP'), xlabel('cardiac cycle'),title('Pearson FR - BP')
        patch(x,y,[0.7 0.7 0.7])
        hold off
        grid on, grid minor
        ylim([yinf ysup])

        subplot(2,8,12:13)
        hold on
        plot(app.spike_res.transduction(int_idxs(i)).cluster(j).fr_rr(:,1), app.spike_res.transduction(int_idxs(i)).cluster(j).fr_rr(:,2),'-o','linewidth',1.5)
        patch(x,y,[0.7 0.7 0.7])
        hold off
        ylabel('Pearson FR - RR'), xlabel('cardiac cycle'),title('Pearson FR - RR')
        grid on, grid minor
        ylim([yinf ysup])

    end

    for j = 1 : length(form_idxs)
        switch form_idxs(j)
            case 1
                h.Visible = 'on';
                savefig(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_transduction_' simple_name(app.popup_int_transduction.Items{int_idxs(i)}) '.fig'],'compact')
                %switch_vis([path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_transduction_' simple_name(app.popup_int_transduction.Items{int_idxs(i)}) '.fig'])
            case 2
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_transduction_' simple_name(app.popup_int_transduction.Items{int_idxs(i)}) '.jpeg'])
            case 3
                saveas(h,[path '\' file '_INT_' num2str(app.settings.interval(1,1)) '-' num2str(app.settings.interval(1,2)) '_transduction_' simple_name(app.popup_int_transduction.Items{int_idxs(i)}) '.epsc'])
        end
    end
    close(h)
end
end 