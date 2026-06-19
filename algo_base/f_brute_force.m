%% Brute Force (Exhaustive Search) Algorithm by Binary Incrementation
% Date: June 2026
% Author : Laure WOLFF
function [best_subpop, best_perf_overall] = f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, other_figs)
    
    % Total number of possible combinations (2^N - 1, ignoring the all-zero mask)
    total_combinations = (2^num_neurons) - 1;
    
    % Safety check to prevent computer freezing if N is set to 125
    if num_neurons > 20
        error('Brute Force aborted: N is too large (%d). Reduce N between 10 and 20 in your main script.', num_neurons);
    end
    
    if showing
        fprintf('-> Launching Brute Force by binary incrementation (%d masks to evaluate...)\n', total_combinations);
    end

    % Variables initialization
    best_perf_overall = -Inf;
    best_mask_overall = zeros(num_neurons, 1);
    
    % Array to store the evaluation history for plotting
    history_perf_brute = zeros(1, total_combinations);

    % Main execution loop
    for i = 1:total_combinations
        
        % Generate the binary string corresponding to the current iteration counter (+1 addition)
        % Example for i = 3 and N = 4: '0011'
        binary_string = dec2bin(i, num_neurons);
        
        % Convert the character string into a numeric vector of 0s and 1s
        current_mask = (binary_string - '0')'; % Column vector [0; 0; 1; 1]
        
        % Compute classification performance P for this specific subpopulation mask
        [perf, ~] = calculate_integrated_P_optimized(CellMatrix, current_mask, ...
            num_stimuli, num_repetitions, t1, t2, metric_choice);
        
        % Store the performance in the history array
        history_perf_brute(i) = perf;
        
        % Check if this binary mask yields the best performance found so far
        if perf > best_perf_overall
            best_perf_overall = perf;
            best_mask_overall = current_mask;
        end
    end

    % Convert the winning binary mask into neuron IDs (e.g., [1, 2, 5])
    best_subpop = find(best_mask_overall == 1)';

    if showing
        fprintf('\n================ BRUTE FORCE CONVERGED ================\n');
        fprintf('Best binary combination found: [%s]\n', num2str(best_subpop));
        fprintf('Absolute maximum performance P = %.4f\n', best_perf_overall);
        fprintf('=======================================================\n');
    end

    %% Plotting the performance evolution
    if other_figs == true && ~isempty(history_perf_brute)
        figure('Name', 'Brute Force - Combinatorial Search History', 'Color', [1 1 1]);
        
        % Plot every tested combination's performance
        plot(1:total_combinations, history_perf_brute, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.8);
        hold on;
        
        % Reconstruct and plot the step-by-step maximum progress line
        best_so_far = cummax(history_perf_brute);
        plot(1:total_combinations, best_so_far, 'b-', 'LineWidth', 2);
        
        % Highlight the global maximum point
        idx_max = find(history_perf_brute == best_perf_overall, 1, 'first');
        plot(idx_max, best_perf_overall, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', [1 0.2 0.2]);
        
        box on; grid on;
        xlim([1, total_combinations]);
        % Adjust Y limits based on data
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