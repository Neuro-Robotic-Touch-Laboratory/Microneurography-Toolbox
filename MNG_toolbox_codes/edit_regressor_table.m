function edit_regressor_table(app,action)

table_data = app.tbl_regressor.Data;

switch action
    case 'add'
        idx = find(strcmp(app.popup_regressor.Value, table_data));
        if isempty(idx)
            table_data{end+1,1} = app.popup_regressor.Value;
            table_data{end,2} = app.edt_lagtime.Value;
        else
            table_data{idx,2} = app.edt_lagtime.Value;
        end
    case 'rem'
        if size(table_data,1) >0
            table_data(app.tbl_regressor.Selection,:) =[]; 
        end
end
app.tbl_regressor.Data = table_data;