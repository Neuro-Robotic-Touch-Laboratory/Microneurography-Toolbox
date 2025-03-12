function update_signal_popup(app)
%UPATE_SIGNAL_POPUP updates the selection options for all popup menues
%for the selection of specific signals. to be called after calulation of
%derived signals
%   Detailed explanation goes here
if ~isempty(app.data)
    name = vertcat(app.data.name);
    if length(app.data) < 2
        i =1;
    else
        i = 2;
    end
else
    name = {'none'};
    i = 1;
end
app.popup_signal_spect.Items = name;
app.popup_signal_spect.Value = name{1};

app.popup_signal1.Items = name;
app.popup_signal1.Value = name{1};

app.popup_signal2.Items = name;
app.popup_signal2.Value = name{i};

% app.popup_entropy_signal.Items = name;
% app.popup_entropy_signal.Value = name{1};

app.popup_coorelation_signal1.Items = name;
app.popup_coorelation_signal1.Value = name{1};

app.popup_coorelation_signal2.Items = name;
app.popup_coorelation_signal2.Value = name{i};

app.popup_signal_to_predict.Items = name;
app.popup_signal_to_predict.Value = name{1};

app.popup_regressor.Items = name;
app.popup_regressor.Value = name{1};

app.popup_lag_sig1.Items = name;
app.popup_lag_sig1.Value = name{1};

app.popup_lag_sig2.Items = name;
app.popup_lag_sig2.Value = name{i};

end

