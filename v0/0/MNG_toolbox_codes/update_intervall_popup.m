function update_intervall_popup(app)
%update_intervall_popup updates the selection options for all popup menues
%for the selection of specific intervals. to be called after intervalls
%were edited
%   Detailed explanation goes here
if ~isempty(app.burst_ints)
    names = {'full'};
    for i = 1 : length(app.burst_ints)
        names{i+1} = app.burst_ints(i).name;
    end   
    app.popup_int.Items = names;
    app.popup_int.Value = app.popup_int.Items{1} ;
    app.popup_int_spike.Items = names;
    app.popup_int_spike.Value = app.popup_int_spike.Items{1};
    
    names = {'full'};
    tmp = find(vertcat(app.burst_ints.type) == 1);
    for i = 1 : length(tmp)
        names{i+1} = app.burst_ints(tmp(i)).name;
    end 
    app.popup_int_spect.Items = names;
    app.popup_int_spect.Value = app.popup_int_spect.Items{1} ;
    app.popup_int_wavelet.Items = names;
    app.popup_int_wavelet.Value = app.popup_int_spect.Items{1} ;
    app.popup_int_entropy.Items = names;
    app.popup_int_entropy.Value = app.popup_int_spect.Items{1} ;
    app.popup_int_stepwise.Items = names;
    app.popup_int_stepwise.Value = app.popup_int_stepwise.Items{1};
    app.popup_int_annotate.Items = names;
    app.popup_int_annotate.Value = app.popup_int_annotate.Items{1};
    app.popup_int_corre.Items = names;
    app.popup_int_corre.Value = app.popup_int_corre.Items{1};
    
    names = {'full'};
    for i = 1 : length(app.burst_ints)
        if app.burst_ints(i).type == 1 && sum(diff(app.burst_ints(i).borders,1,2)) >= 60
            names{end+1} = app.burst_ints(i).name;
        end
    end
    app.popup_int_transduction.Items = names;
    app.popup_int_transduction.Value = app.popup_int_transduction.Items{1};

else
    names = {'full'};
    app.popup_int.Items = names;
    app.popup_int.Value = app.popup_int.Items{1} ;
    app.popup_int_spike.Items = names;
    app.popup_int_spike.Value = app.popup_int_spike.Items{1};
    app.popup_int_spect.Items = names;
    app.popup_int_spect.Value = app.popup_int_spect.Items{1} ;
    app.popup_int_wavelet.Items = names;
    app.popup_int_wavelet.Value = app.popup_int_spect.Items{1} ;
    app.popup_int_entropy.Items = names;
    app.popup_int_entropy.Value = app.popup_int_spect.Items{1} ;
    app.popup_int_stepwise.Items = names;
    app.popup_int_stepwise.Value = app.popup_int_stepwise.Items{1};
    app.popup_int_annotate.Items = names;
    app.popup_int_annotate.Value = app.popup_int_annotate.Items{1};
    app.popup_int_corre.Items = names;
    app.popup_int_corre.Value = app.popup_int_corre.Items{1};
    app.popup_int_transduction.Items = names;
    app.popup_int_transduction.Value = app.popup_int_transduction.Items{1};
end

end
