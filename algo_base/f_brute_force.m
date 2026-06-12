%% Brute force wrapper calling the native C MEX function
% Date: June 2026
% Author : Laure WOLFF

function  f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, plotting)
    
    total_combinations = (2^num_neurons) - 1;
    
    % Safety cap to prevent hardware lock
    if num_neurons > 22
        error('Brute Force aborted: N is too large (%d). Reduce N between 1 and 21.', num_neurons);
    end
    
    if showing
        fprintf('-> Launching Brute Force via C-MEX Binary Engine (%d masks to evaluate...)\n', total_combinations);
    end

    %% CALL THE MEX ENGINE
    % [best_subpop, best_perf_overall, history_perf_brute]
    [best_subpop, best_perf_overall, history_perf_brute] = f_brute_force_mex(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);

    if showing
        fprintf('\n================ BRUTE FORCE CONVERGED ================\n');
        fprintf('Best binary combination found: [%s]\n', num2str(best_subpop));
        fprintf('Absolute maximum performance P = %.4f\n', best_perf_overall);
        fprintf('=======================================================\n');
    end

    %% Plotting the performance evolution
    if plotting == true && ~isempty(history_perf_brute)
        figure('Name', 'Brute Force - Combinatorial Search History', 'Color', [1 1 1]);
        
        plot(1:total_combinations, history_perf_brute, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        hold on;
        
        best_so_far = cummax(history_perf_brute);
        plot(1:total_combinations, best_so_far, 'b-', 'LineWidth', 2);
        
        idx_max = find(history_perf_brute == best_perf_overall, 1, 'first');
        plot(idx_max, best_perf_overall, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', [1 0.2 0.2]);
        
        box on; grid on;
        xlim([1, total_combinations]);
        min_p = max(0, min(history_perf_brute));
        max_p = max(history_perf_brute);
        ylim([min_p, max(max_p * 1.1, 0.1)]);
        
        xlabel('Binary Counter Iterations (Search Space)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Performance P', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Brute Force Search Tree Exploration (N = %d Neurons)', num_neurons), 'FontSize', 13, 'FontWeight', 'bold');
        
        legend('Evaluated Mask Performance', 'Global Maximum Progress', 'Absolute Best Solution', ...
               'Location', 'SouthEast');
        hold off;
        shg;
    end
end



