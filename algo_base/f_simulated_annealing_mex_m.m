%% Simulated Annealing algorithm script (for SP hypothesis) 
% Date: June 2026
% Author : Laure WOLFF

function [nb_iterations] = f_simulated_annealing_mex_m(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, plotting )

[best_mask_overall, best_perf_overall, nb_iterations, Matrix_Grid, history_perf, ...
 hist_iter_P, hist_iter_bestP, hist_iter_size, hist_iter_temp] = ...
    f_simulated_annealing_core_mex(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing);

%% Final Wrap-up
best_subpop = find(best_mask_overall == 1)';
if showing
    fprintf('Optimal subpopulation found: [%s]\n', num2str(best_subpop));
    fprintf('Max performance P = %.4f\n', best_perf_overall);
    fprintf('number of iteration %.4f\n', nb_iterations);
end

%% Plotting 
if plotting == true && ~isempty(Matrix_Grid)
    num_paliers_reals = size(Matrix_Grid, 1);
    
    figure('Name', 'Results - Simulated Annealing','Color', [1 1 1]);
    
    % Matrix showing masks
    subplot(1, 4, 1:3);
    imagesc(1:num_neurons, 1:num_paliers_reals, Matrix_Grid);
    mymap = [0.2 0.4 0.8; 0.9 0.2 0.2]; 
    colormap(gca, mymap); 
    clim([0, 1]);
    set(gca, 'YTick', 1:num_paliers_reals)
    set(gca, 'YDir', 'normal'); 
    hold on;
    
    for i = 1:length(best_subpop)
        plot(best_subpop(i), num_paliers_reals, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    box on;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 10);
    set(gca, 'XTick', 1:num_neurons);
    xlabel('Neuron ID', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Temperature Steps (Cooling)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Simulated Annealing - Selected Neurons History', 'FontSize', 13, 'FontWeight', 'bold');
    
    cb = colorbar;
    set(cb, 'Ticks', [0.25, 0.75], 'TickLabels', {'Desactivated (0)', 'Activated (1)'}, 'FontSize', 10);
    
    subplot(1, 4, 4);
    plot(history_perf, 1:num_paliers_reals, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
    hold on;
    
    plot(best_perf_overall, num_paliers_reals, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
    plot(best_perf_overall, num_paliers_reals, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    grid on; box on;
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 1.2, 'YTickLabel', [], 'FontSize', 10);
    ylim([0.5, num_paliers_reals + 0.5]);
    xlabel('Performance P', 'FontSize', 12, 'FontWeight', 'bold');
    title('P(temp)', 'FontSize', 13, 'FontWeight', 'bold');
    hold off; 
    
    % Diagnostics de convergence
    figure('Name', 'Simulated Annealing - Convergence Diagnostics', 'Color', 'w', 'Position', [150, 150, 1000, 400]);
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    
    nexttile
    plot(hist_iter_P, 'k', 'LineWidth', 1.2);
    hold on;
    plot(hist_iter_bestP, 'r', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Performance P');
    legend({'Current', 'Best'}, 'Location', 'best');
    title('Discrimination performance');
    grid on; box on;
    
    nexttile
    plot(hist_iter_size, 'b', 'LineWidth', 1.2);
    hold on;
    yline(length(best_subpop), '--r', 'Optimal found', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Population size');
    title('Subpopulation Size Track');
    grid on; box on;
    
    nexttile
    semilogy(hist_iter_temp, 'm', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Temperature');
    title('Cooling Schedule (Log Scale)');
    grid on; box on;
    
    sgtitle(sprintf('Simulated Annealing Dynamics | Best P = %.4f | Best Population = [%s]',...
        best_perf_overall, num2str(sort(best_subpop))));
    shg;
end
end