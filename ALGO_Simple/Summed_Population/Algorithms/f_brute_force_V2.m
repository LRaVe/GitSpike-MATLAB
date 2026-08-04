%% Brute Force (Exhaustive Search) Algorithm by Binary Incrementation
% Date: June 2026
% Author : Laure WOLFF

function [best_subpop, best_perf_overall] = f_brute_force_V2(spikes,Tmax,Distances,threshold, showing,other_figs, useMex)
% F_BRUTE_FORCE_V2 Evaluates all :math:`2^N - 1` subpopulation combinations to find optimal neuronal subsets.
%
%   Performs an exhaustive search (Brute Force) over the entire combinatorial space 
%   of neuron subpopulations using binary incrementation. For each candidate subset, 
%   it calculates the classification performance :math:`P` using spike train distances.
%
%   The total number of tested combinations is defined as:
%
%   .. math::
%
%      N_{comb} = 2^{N_{neurons}} - 1
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Standard execution without MEX acceleration or figures
%      [best_subpop, best_perf] = f_brute_force_V2(spikes, Tmax, Distances, 0, 1, 0, false);
%
%      % Fast execution using C/MEX compilation with plotting enabled
%      [best_subpop, best_perf] = f_brute_force_V2(spikes, Tmax, Distances, 0.05, 1, 1, true);
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike times.
%   :type spikes: cell or double
%   :param Tmax: Upper temporal bound of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric type (e.g., `'SPIKE'`, `'RI-SPIKE'`, `'SPIKE_adaptive'`, `'RI-SPIKE adaptative'`).
%   :type Distances: 0 or 1 array  (e.g [1 0 0 0] for the SPIKE distance)
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%   :param showing: Flag to enable/disable console output messages (`1` = active, `0` = quiet).
%   :type showing: logical or integer
%   :param other_figs: Flag to control generation of the combinatorial search history figure (`true`/`false`).
%   :type other_figs: logical
%   :param useMex: Flag to delegate processing to C/MEX compiled function (`f_brute_force_mex`) for high speed.
%   :type useMex: logical
%
%   :returns: 
%             * **best_subpop** (*vector*) -- Indices of the neurons forming the optimal subpopulation.
%             * **best_perf_overall** (*double*) -- Highest classification performance :math:`P` achieved across all combinations.
%
%   .. note::
%      A safety check aborts execution if :math:`N_{neurons} > 20` to prevent memory allocation failure 
%      and exponential execution slowdowns (:math:`2^{20} = 1\,048\,575` iterations).
%
%   .. note::
%      This algorithms has a MEX version
%
%   :Author: Laure WOLFF
%   :Date: June 2026
       
    [num_neurons,num_stimuli, num_repetitions] = size(spikes);

    % Total number of possible combinations (2^N - 1, ignoring the all-zero mask)
    total_combinations = (2^num_neurons) - 1;
    
    % Safety check to prevent computer freezing if N is set too high
    if num_neurons > 20
        error('Brute Force aborted: N is too large (%d). Reduce N between 10 and 20 in your main script.', num_neurons);
    end
    
    if showing
        fprintf('-> Launching Brute Force by binary incrementation (%d masks to evaluate...)\n', total_combinations);
    end
    
    % Variables initialization
    best_perf_overall = -Inf;
    best_subpop = [];
    
    % Array to store the evaluation history for plotting
    history_perf_brute = zeros(1, total_combinations);
    % =====================================================================
    %  Mex version if mex_variable is true
    % =====================================================================
    if useMex 

        t1 = 0; 
        t2 = Tmax;

        [best_subpop, best_perf_overall, history_perf_brute] = f_brute_force_mex(spikes, num_neurons, num_stimuli, num_repetitions, t1, t2, Distances);

        if showing
            fprintf('\n================ BRUTE FORCE (MEX) CONVERGED ================\n');
            fprintf('Best binary combination found: [%s]\n', num2str(best_subpop));
            fprintf('Absolute maximum performance P = %.4f\n', best_perf_overall);
            fprintf('=============================================================\n');
        end
    else
        % =====================================================================
        % STANDARD SEQUENTIAL LOOP WITH FAST BINARY EXTRACTION
        % =====================================================================
        for i = 1:total_combinations
            % FAST BITWISE EXTRACTION: Direct numeric vector creation via bitget
            population = find(bitget(i,1:num_neurons));
    
            % if showing
            %     fprintf('Population : ');
            %     fprintf('%5d',population);
            %     fprintf('\n');
            % end
            
            % Compute classification performance P for this specific subpopulation mask
            perf = evaluate_population( ...
                        spikes,...
                        population,...
                        Tmax,...
                        Distances,...
                        threshold);
    
            % Store the performance in the history array
            history_perf_brute(i) = perf;
            
            % Check if this binary mask yields the best performance found so far
            if perf > best_perf_overall
                best_perf_overall = perf;
                best_subpop = population;
            end
        end
        % =====================================================================
        
        if showing
            fprintf('\n================ BRUTE FORCE CONVERGED ================\n');
            fprintf('Best binary combination found: [%s]\n', num2str(best_subpop));
            fprintf('Absolute maximum performance P = %.4f\n', best_perf_overall);
            fprintf('=======================================================\n');
        end
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
        
        % Adjust Y limits based on data dynamically
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