function update_cluster_popup(app)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

str_cell = {[]};
for i = 1: max(app.spike_res.cluster)
    str_cell{i} = num2str(i);
end
app.popup_cluster_transduction.Items = str_cell;
app.popup_cluster_transduction.Value = app.popup_cluster_transduction.Items{1};
app.popup_cluster.Items = str_cell;
app.popup_cluster.Value = app.popup_cluster.Items{1};
app.popup_cluster_transduction.Items = str_cell;
app.popup_cluster_transduction.Value = app.popup_cluster_transduction.Items{1};
end

