function general_plot(X,Y,num_plots,legends,xlabeling,ylabeling,titling)

% Create figure
figure;

% Create axes
axes1 = axes;
hold(axes1,'on');
num_plots
for i = 1:num_plots
    plot(X(i,:),Y(i,:),'DisplayName',legends{i},'LineWidth',1.5);
end

% Create xlabel
xlabel(xlabeling,...
    'LineWidth',1,...
    'Interpreter','latex');

% Create title

title(titling,'Interpreter','latex');
  

% Create ylabel
ylabel(ylabeling,...
    'Interpreter','latex');

box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',18);
% Create legend
legend1 = legend(axes1,'show');
set(legend1,'FontSize',16,'FontName','MS Sans Serif');
grid on;
grid minor;


