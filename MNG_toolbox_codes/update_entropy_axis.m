function update_entropy_axis(app)
%%update_entropy_axis plots the result of the entropy analysis
cla(app.ax_entropy1)
cla(app.ax_entropy2)
cla(app.ax_entropy3)
cla(app.ax_entropy4)

int_idx = find(strcmp(app.popup_int_entropy.Value,app.popup_int_entropy.Items));
int_name = app.popup_int_entropy.Value;

channel_idx = find(strcmp(app.popup_entropy_signal.Value,app.popup_signal_spect.Items));
[data,ts,~, ~] = current_signal(app, channel_idx);
data(:,2) = ts(1):ts(1):ts(2);

if int_idx == 1
    borders = [1,size(data,1)];
else
    tmp =find(vertcat(app.burst_ints.type) == 1);
    borders = app.burst_ints(tmp(int_idx-1)).borders;
    borders = [ceil(borders(1)/ts(1)), floor(borders(2)/ts(1))];
end    
data_int = data(borders(1):borders(2),:);

plot(app.ax_entropy1, data_int(:,2), data_int(:,1), 'k', 'LineWidth', .2 );
app.lbl_entropy1.Text = ['Original time series ' int_name];


if ~isempty(app.entropy_res)
    switch app.entropy_res(int_idx).cfg.method
        case 'PE'
            plot(app.ax_entropy2, data_int( end - length( app.entropy_res(int_idx).outdata) + 1:end, 2 ), ...
                  app.entropy_res(int_idx).outdata, 'k', 'LineWidth', .2 ); 
            app.lbl_entropy2.Text = ['Values of permutation entropy ' int_name];
            app.lbl_entropy2.Visible = 'on';
            app.ax_entropy2.Visible = 'on';
            app.ax_entropy3.Visible = 'off';
            app.ax_entropy4.Visible = 'off';
            app.lbl_entropy3.Visible = 'off';
            app.lbl_entropy4.Visible = 'off';
            app.lbl_entropy_unit1.Visible = 'on';
            app.lbl_entropy_unit2.Visible = 'off';
        case 'PEeq'
            plot(app.ax_entropy2, data_int( end - length( app.entropy_res(int_idx).outdata) + 1:end, 2 ), ...
                  app.entropy_res(int_idx).outdata, 'k', 'LineWidth', .2 ); 
            app.lbl_entropy2.Text = ['Values of permutation entropy for ordinal patterns with tied ranks ' int_name];
            app.lbl_entropy2.Visible = 'on';
            app.ax_entropy2.Visible = 'on';
            app.ax_entropy3.Visible = 'off';
            app.ax_entropy4.Visible = 'off';
            app.lbl_entropy3.Visible = 'off';
            app.lbl_entropy4.Visible = 'off';
            app.lbl_entropy_unit1.Visible = 'on';
            app.lbl_entropy_unit2.Visible = 'off';
        case 'CE'
            plot(app.ax_entropy2, data_int( end - length( app.entropy_res(int_idx).outdata) + 1:end, 2 ), ...
                  app.entropy_res(int_idx).outdata, 'k', 'LineWidth', .2 ); 
            app.lbl_entropy2.Text = ['Values of conditional entropy of ordinal patterns ' int_name];
            app.lbl_entropy2.Visible = 'on';
            app.ax_entropy2.Visible = 'on';
            app.ax_entropy3.Visible = 'off';
            app.ax_entropy4.Visible = 'off';
            app.lbl_entropy3.Visible = 'off';
            app.lbl_entropy4.Visible = 'off';
            app.lbl_entropy_unit1.Visible = 'on';
            app.lbl_entropy_unit2.Visible = 'off';
        case 'rePE'
            plot(app.ax_entropy2, data_int( end - length( app.entropy_res(int_idx).outdata) + 1:end, 2 ), ...
                  app.entropy_res(int_idx).outdata, 'k', 'LineWidth', .2 ); 
            app.lbl_entropy2.Text = ['Values of robust permutation entropy ' int_name];
            app.lbl_entropy2.Visible = 'on';
            app.ax_entropy2.Visible = 'on';
            app.ax_entropy3.Visible = 'off';
            app.ax_entropy4.Visible = 'off';
            app.lbl_entropy3.Visible = 'off';
            app.lbl_entropy4.Visible = 'off';
            app.lbl_entropy_unit1.Visible = 'on';
            app.lbl_entropy_unit2.Visible = 'off';
        case 'opdPE'
            plot(app.ax_entropy2, data_int( end - length( app.entropy_res(int_idx).outdata.ePe) + 1:end, 2 ), ...
                  app.entropy_res(int_idx).outdata.ePe, 'k', 'LineWidth', .2 ); 
            app.lbl_entropy2.Text = ['Values of ordinary pattern distribution permutation entropy ' int_name];
            app.lbl_entropy2.Visible = 'on';
            app.ax_entropy2.Visible = 'on';
            app.ax_entropy3.Visible = 'off';
            app.ax_entropy4.Visible = 'off';
            app.lbl_entropy3.Visible = 'off';
            app.lbl_entropy4.Visible = 'off';
            app.lbl_entropy_unit1.Visible = 'on';
            app.lbl_entropy_unit2.Visible = 'off';
        case 'all'
            plot(app.ax_entropy1, data_int( end - length( app.entropy_res(int_idx).outdata.PE ) + 1:end, 2), data_int( end - length( app.entropy_res(int_idx).outdata.PE ) + 1:end, 1), 'k', 'LineWidth', 0.2 )
            app.lbl_entropy1.Text = ['Values of ordinary pattern distribution permutation entropy ' int_name];

            plot(app.ax_entropy2, data_int( end - length( app.entropy_res(int_idx).outdata.PE ) + 1:end, 2), app.entropy_res(int_idx).outdata.PE, 'k', 'LineWidth', 0.2 )
            app.lbl_entropy2.Text = ['Values of ordinary pattern distribution permutation entropy ' int_name];

            plot(app.ax_entropy3, data_int( end - length( app.entropy_res(int_idx).outdata.CE ) + 1:end, 2), app.entropy_res(int_idx).outdata.CE, 'k', 'LineWidth', 0.2 )
            app.lbl_entropy3.Text = ['Values of ordinary pattern distribution permutation entropy ' int_name];

            plot(app.ax_entropy4,  data_int( end - length( app.entropy_res(int_idx).outdata.RE ) + 1:end, 2), app.entropy_res(int_idx).outdata.RE, 'k', 'LineWidth', 0.2 )
            app.lbl_entropy4.Text = ['Values of ordinary pattern distribution permutation entropy ' int_name];
            app.lbl_entropy2.Visible = 'on';
            app.ax_entropy2.Visible = 'on';
            app.ax_entropy3.Visible = 'on';
            app.ax_entropy4.Visible = 'on';
            app.lbl_entropy3.Visible = 'on';
            app.lbl_entropy4.Visible = 'on';
            app.lbl_entropy_unit1.Visible = 'off';
            app.lbl_entropy_unit2.Visible = 'on';
    end
else
    app.lbl_entropy2.Visible = 'off';
    app.ax_entropy2.Visible = 'off';
    app.ax_entropy3.Visible = 'off';
    app.ax_entropy4.Visible = 'off';
    app.lbl_entropy3.Visible = 'off';
    app.lbl_entropy4.Visible = 'off';
    app.lbl_entropy_unit1.Visible = 'off';
    app.lbl_entropy_unit2.Visible = 'off';
end
end

      