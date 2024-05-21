% % Button down function: EntropyTab
% function EntropyTabButtonDown(app, event)
%     working(app,'on')
%     auto_derived(app)
%     working(app,'off')
% end
 
 
% lbl_entropy_unit1
%     Text ='seconds'
%     Position =  [502,340,54,22]
 
% lbl_entropy_unit2
%     Text ='seconds'
%     Position =  [496,10,54,22]
 
% lbl_entropy1
%     Position =  [60,666,931,22]
%     FontSize = 15
% lbl_entropy2
%     Position =  [60,494,931,22]
%     FontSize = 15
% lbl_entropy3
%     Position =  [60,334,931,22]
%     FontSize = 15
% lbl_entropy4
%     Position =  [54,174,931,22]
%     FontSize = 15
 
% popup_int_entropy
%     Position =  [1020,250,240,30]
%         Position =  [1020,610,240,30]   
%     FontSize = 15
%     % Value changed function: edt_window_size
%     function edt_window_sizeValueChanged(app, event)
%         working(app,'on')
%         app.entropy_res = [];
%         update_entropy_axis(app)
%         app.btn_entropy.BackgroundColor = [.96, .96, .96];
%         working(app,'off')
%     end
%     >>>  update_intervall_popup(app)
 
% parametersPanel
%     BackgroundColor = [0.80,0.95,0.91]
%     FontSize = 15
%     Position =  [1020,320,240,240]
 
%     edt_order
%         Value = 3
%         Limits = -Inf,Inf
%         ValueDisplayFormat = %11.4g 
%         FontSize = 15
%         Position =  [120,170,100,30]
%             Position =  [20,170,100,30]
%         % Value changed function: edt_order
%         function edt_orderValueChanged(app, event)
%             working(app,'on')
%             app.entropy_res = [];
% %             cla(app.ax_entropy2)
% %             cla(app.ax_entropy3)
% %             cla(app.ax_entropy4)
% %             app.ax_entropy3.Visible = 'off';
% %             app.ax_entropy4.Visible = 'off';
% %             app.lbl_entropy3.Visible = 'off';
% %             app.lbl_entropy4.Visible = 'off';
% %             app.lbl_entropy_unit3.Visible = 'off';
% %             app.lbl_entropy_unit2.Visible = 'off';
%             update_entropy_axis(app)
%             app.btn_entropy.BackgroundColor = [.96, .96, .96];
%             working(app,'off')
%         end
 
%     edt_delay
%         Value = 2
%         Limits = -Inf,Inf
%         ValueDisplayFormat = %11.4g 
%         FontSize = 15
%         Position =  [120,70,100,30]
%             Position =  [20,120,100,30]
%         % Value changed function: edt_delay
%         function edt_delayValueChanged(app, event)
%             working(app,'on')
%             app.entropy_res = [];
%             update_entropy_axis(app)
%             app.btn_entropy.BackgroundColor = [.96, .96, .96];
%             working(app,'off')
%         end
 
%     edt_order_seq
%         Value = 6
%         Limits = 1,Inf
%         ValueDisplayFormat = %11.4g 
%         FontSize = 15
%         Position =  [120,120,100,30]
%             Position =  [20,70,100,30]
%         % Value changed function: edt_order_seq
%         function edt_order_seqValueChanged(app, event)
%             working(app,'on')
%             app.entropy_res = [];
%             update_entropy_axis(app)
%             app.btn_entropy.BackgroundColor = [.96, .96, .96];
%             working(app,'off')
%         end
 
%     edt_window_size
%         Value = 512
%         Limits = 0,Inf
%         ValueDisplayFormat = %11.4g 
%         FontSize = 15
%         Position =  [120,20,100,30]
%             Position =  [20,20,100,30]
%         % Value changed function: edt_window_size
%         function edt_window_sizeValueChanged(app, event)
%             working(app,'on')
%             app.entropy_res = [];
%             update_entropy_axis(app)
%             app.btn_entropy.BackgroundColor = [.96, .96, .96];
%             working(app,'off')
%         end
 
% btn_save_entropy_results
%     Text = 'save results'
%     FontSize = 15
%     Position =  [1020,90,240,30]
%     % Button pushed function: btn_save_enrtopy_results
%     function btn_save_enrtopy_resultsButtonPushed(app, event)
%         if strcmp(app.settings.output_dir , 'none')
%             warndlg('Please select output directory first!')
%         else    
%             working(app,'on')
%             save_entropy_res(app)
%             working(app,'off')
%         end
%     end
 
% popup_method
%     FontSize = 15
%     Value = {'Permutation entropy'}
%     Items = {'Permutation entropy', 'permutation entropy and ordinal distributions', 'all implemented ordinal-patterns-based measures', 'conditional entropy', 'robust permutation entropy', 'permutation entropy with tied ranks'}
%     Position =  [1020,90,240,30]
%         Position =  [1020,210,240,30]
 
% popup_entropy_signal
%     FontSize = 15
%     Position =  [1020,250,240,30]
%         Position =  [1020,280,240,30]
%     % Value changed function: popup_entropy_signal
%     function popup_entropy_signalValueChanged(app, event)
%         working(app,'on')
%         app.entropy_res = [];
% %             app.ax_entropy3.Visible = 'off';
% %             app.ax_entropy4.Visible = 'off';
% %             app.lbl_entropy3.Visible = 'off';
% %             app.lbl_entropy4.Visible = 'off';
% %             app.lbl_entropy_unit3.Visible = 'off';
% %             app.lbl_entropy_unit2.Visible = 'off';
%         update_entropy_axis(app)
%         app.btn_entropy.BackgroundColor = [.96, .96, .96];
%         working(app,'off')
%     end
%     >> update_signal_popup(app)
 
% btn_entropy
%     Text = 'entropy'
%     FontSize = 15
%     Position =  [1020,130,240,30]
%     % Button pushed function: btn_entropy
%     function btn_entropyButtonPushed(app, event)
%         working(app,'on')
%         app.entropy_res = entropy_analysis(app);
%         update_entropy_axis(app)
%         app.btn_entropy.BackgroundColor = 'g';
%         working(app,'off')
%     end
 
% ax_entropy1
%     Position = [20,515,980,155]
%     >>> clear_axis(app)

% ax_entropy2
%     Position = [20,185,980,155]
%     >>> clear_axis(app)
 
% ax_entropy3
%     Position = [20,515,980,155]
%     >>> clear_axis(app)

% ax_entropy4
%     Position = [20,20,980,155]
%     >>> clear_axis(app)