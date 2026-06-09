%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function plot_simulated_annealing(result,codingNeurons)

figure('Color','w');

tiledlayout(1,3,'TileSpacing','compact');

%% ==========================================
%% PERFORMANCE
%% ==========================================

nexttile

plot(result.history.P,...
     'k',...
     'LineWidth',1.5);

hold on

plot(result.history.bestP,...
     'r',...
     'LineWidth',2);

xlabel('Iteration');
ylabel('P');

legend({'Current','Best'},...
       'Location','best');

title('Discrimination performance');

grid on


%% ==========================================
%% POPULATION SIZE
%% ==========================================

nexttile

plot(result.history.size,...
     'b',...
     'LineWidth',1.5);

hold on

yline(length(codingNeurons),...
      '--r',...
      'Ground truth');

xlabel('Iteration');
ylabel('Population size');

title('Population size');

grid on


%% ==========================================
%% TEMPERATURE
%% ==========================================

nexttile

semilogy(result.history.temperature,...
         'm',...
         'LineWidth',1.5);

xlabel('Iteration');
ylabel('Temperature');

title('Cooling schedule');

grid on


%% ==========================================
%% GLOBAL TITLE
%% ==========================================

sgtitle(sprintf( ...
    'Simulated Annealing | Best P = %.3f | Best Population = [%s]',...
    result.bestP,...
    num2str(sort(result.bestPopulation))));

end