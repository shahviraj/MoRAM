function construct_subplots(err_data,pr,savename,plot_method,save_data)

% prepare the data for plotting according to the plot_method:
switch plot_method
    case 'mean-error'
       Y1 = squeeze(mean(err_data,3));
    case 'median-error'
       Y1 = squeeze(median(err_data,3));
    case 'max-error'
        Y1 = squeeze(max(err_data,3));
    case 'confidence'
        binary_err = (err_data < 0.0001);
        Y1 = squeeze(sum(binary_err,3))/pr.num_trials;
end



% Create figure
figure;
% 
% % Create axes
% axes1 = subplot(1,2,1);
% hold(axes1,'on');

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create plot
for i = 1:size(Y1,2)
    plot(pr.mspan,Y1(1:length(pr.mspan),i),'-s','DisplayName',['s=',num2str(pr.s_span(i))],'LineWidth',2);
end
% Create xlabel
xlabel({'\textbf{Number of measurements} $\mathbf{(m)}$'},...
    'LineWidth',1,...
    'Interpreter','latex');

% Create title
% switch method
%     case 'cosamp'
 title(['\textbf{Relative reconstruction }','\textbf{',pr.plot_method,' vs no. of measurements; for }','\textbf{',pr.method,' with} $\mathbf{||x^*||=',num2str(pr.amp),' ,R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
    'Interpreter','latex');
%     case 'robust-cosamp'
%     title(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
%     'Interpreter','latex');
% end

% Create ylabel
ylabel({'\textbf{Reconstruction Error}; $\mathbf{\frac{||x^*-x||}{||x^*||}}$'},...
    'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',18);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',20,'FontName','MS Sans Serif');
grid on;
grid minor;
ylim([0,2]);
%yticks(0:0.1:1.5)
pbaspect([1 1 1]);
hold(axes1,'off')

% axes2 = subplot(1,2,2);
% hold(axes2,'on');
% 
% % Create plot
% for i = 1:size(Y1,2)
%     plot(pr.mspan2,Y1(length(pr.mspan1)+1:end,i),'-s','DisplayName',['s=',num2str(pr.s_span(i))],'LineWidth',2);
% end
% % Create xlabel
% xlabel({'\textbf{Number of measurements} $\mathbf{(m)}$'},...
%     'LineWidth',1,...
%     'Interpreter','latex');
% 
% % Create title
% % switch method
% %     case 'cosamp'
% % title(['\textbf{Relative reconstruction error vs number of measurements; for CoSaMP with} $\mathbf{R=',num2str(pr.R), ',n=',num2str(pr.n),'}$'],...
% %     'Interpreter','latex');
% %     case 'robust-cosamp'
% %     title(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
% %     'Interpreter','latex');
% % end
% 
% 
% 
% 
% % Create ylabel
% %ylabel({'\textbf{Reconstruction Error}; $\mathbf{\frac{||x^*-x||}{||x^*||}}$'},...
% %    'Interpreter','latex');
% 
% box(axes2,'on');
% % Set the remaining axes properties
% set(axes2,'FontSize',16);
% % Create legend
% legend1 = legend(axes2,'show');
% set(legend1,'FontSize',20,'FontName','MS Sans Serif');
% grid on;
% grid minor;
% ylim([0,1]);
% %yticks(0:0.1:1.5)
% pbaspect([1 1 1]);
% hold(axes2,'off');

% p = mtit(['\textbf{Relative reconstruction }','\textbf{',pr.plot_method,' vs no. of measurements; for }','\textbf{',pr.method,' with} $\mathbf{||x^*||=',num2str(pr.amp),' ,R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
%     'Interpreter','latex','FontSize',22,'xoff',0.0,'yoff',-0.1);
% switch pr.method
%     case 'cosamp'
%  p = mtit(['\textbf{Relative reconstruction error vs number of measurements; for}','\textbf{',pr.method,' with} $\mathbf{||x^*||=',num2str(pr.amp),'}$'],...
%     'Interpreter','latex','FontSize',22,'xoff',0.0,'yoff',-0.1);
%     case 'robust-cosamp'
%  p =  mtit(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{||x^*||=',num2str(pr.amp),'}$'],...
%     'Interpreter','latex','FontSize',22,'xoff',0.0,'yoff',-0.1);
% end
savepath = ['./results/',savename,num2str(datenum(datetime('now')))];
if save_data == 1
    
    save([savepath,'.mat'],'err_data');
end

savefig([savepath,'_',pr.plot_method,'_single.fig']);


