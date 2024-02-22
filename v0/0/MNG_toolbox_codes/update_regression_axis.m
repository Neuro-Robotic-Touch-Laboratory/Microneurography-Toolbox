function update_regression_axis(app)
if ~isempty(app.stepwise_res)
    int_items = find(strcmp(app.popup_int_stepwise.Value , app.popup_int_stepwise.Items));
    cla(app.ax_regression1)
    cla(app.ax_regression2)
    BlandAltman_ax(app.ax_regression1, app.ax_regression2, ...
                   app.stepwise_res(int_items).data.par, ...
                   app.stepwise_res(int_items).data.y_estimates)
end