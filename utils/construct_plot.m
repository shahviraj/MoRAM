function construct_plot(Y1,pr,savename, method)

% Create figure
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');

% Create plot
for i = 1:size(Y1,2)
    plot(pr.mspan,Y1(:,i),'DisplayName',['s=',num2str(pr.s_span(i))],'LineWidth',1.5);
end
% Create xlabel
xlabel({'\textbf{Number of measurements} $\mathbf{(m)}$'},...
    'LineWidth',1,...
    'Interpreter','latex');

% Create title
switch method
    case 'cosamp'
title(['\textbf{Relative reconstruction error vs number of measurements; for CoSaMP with} $\mathbf{R=',num2str(pr.R), ',n=',num2str(pr.n),'}$'],...
    'Interpreter','latex');
    case 'robust-cosamp'
    title(['\textbf{Relative reconstruction error vs number of measurements; for robust CoSaMP with} $\mathbf{R=',num2str(pr.R),',n=',num2str(pr.n),'}$'],...
    'Interpreter','latex');
end

% Create ylabel
ylabel({'\textbf{Reconstruction Error}; $\mathbf{\frac{||x^*-x||}{||x^*||}}$'},...
    'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',16,'FontName','MS Sans Serif');
grid on;
grid minor;
%yticks(0:0.1:1.5)
savefig(['./results/mod_recovery_results/',savename,'.fig'])

