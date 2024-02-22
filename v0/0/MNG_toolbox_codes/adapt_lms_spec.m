function adapt_lms_spec(app,event)
%ADAPT_LMS_SPEC Summary of this function goes here
%   Detailed explanation goes here

xl = event.AffectedObject.XLim;
yl = [0 app.edt_amp.Value];

switch event.AffectedObject.Tag
    case app.ax_spect_1_1.Tag
        if ~isequal(app.ax_spect_2_1.XLim, xl)
            app.ax_spect_2_1.XLim = xl;
        end
    case app.ax_spect_2_1.Tag
        if ~isequal(app.ax_spect_3_1.XLim, xl)
            app.ax_spect_3_1.XLim = xl;
        end
    case app.ax_spect_3_1.Tag
        if ~isequal(app.ax_spect_4_1.XLim, xl)
            app.ax_spect_4_1.XLim = xl;
        end
    case app.ax_spect_4_1.Tag
        if ~isequal(app.ax_spect_1_2.XLim, xl)
            app.ax_spect_1_2.XLim = xl;
            app.ax_spect_1_2.YLim = yl;
        end
    case app.ax_spect_1_2.Tag
        if ~isequal(app.ax_spect_2_2.XLim, xl)
            app.ax_spect_2_2.XLim = xl;
            app.ax_spect_2_2.YLim = yl;
        end
    case app.ax_spect_2_2.Tag
        if ~isequal(app.ax_spect_3_2.XLim, xl)
            app.ax_spect_3_2.XLim = xl;
            app.ax_spect_3_2.YLim = yl;
        end
    case app.ax_spect_3_2.Tag
        if ~isequal(app.ax_spect_4_2.XLim, xl)
            app.ax_spect_4_2.XLim = xl;
            app.ax_spect_4_2.YLim = yl;
        end
    case app.ax_spect_4_2.Tag
        if ~isequal(app.ax_spect_1_1.XLim, xl)
            app.ax_spect_1_1.XLim = xl;
        end
end


