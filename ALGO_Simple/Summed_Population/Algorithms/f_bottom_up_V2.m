%% Bottom-up algorithm script 
% Date: May-June 2026
% Author : Laure WOLFF

function f_bottom_up_V2(spikes,Tmax,Distances,threshold,showing,plotting,other_figs)
% F_BOTTOM_UP_V2 Sequential greedy inclusion algorithm for optimal neuronal subpopulation selection.
%
%   Performs a iterative Bottom-Up (forward selection) search to identify 
%   the subpopulation of neurons that maximizes the classification performance :math:P. 
%   Starting from an empty set, at each step :math:`k`, the algorithm evaluates all available 
%   candidate neurons and permanently includes the one that yields the highest incremental performance.
%
%   The computational complexity of this forward selection process is polynomial:
%
%   .. math::
%
%      \mathcal{O}\left(\frac{N(N+1)}{2}\right)
%
%   which drastically reduces the evaluation space compared to the :math:`2^N - 1` exhaustive search.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Run Bottom-Up selection with console display and paper figure generation
%      f_bottom_up_V2(spikes, Tmax, Distances, 0, true, true, false);
%
%      % Quiet execution without figure generation
%      f_bottom_up_V2(spikes, Tmax, Distances, 'auto', false, false, false);
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric type (e.g., `'SPIKE'`, `'RI-SPIKE'`, `'SPIKE_adaptive'`, `'RI-SPIKE adaptative'`).
%   :type Distances: 0 or 1 array  (e.g [1 0 0 0] for the SPIKE distance)
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%   :param showing: Flag to enable/disable detailed console output and progression logs.
%   :type showing: logical
%   :param plotting: Master flag to enable figure generation.
%   :type plotting: logical
%   :param other_figs: Secondary flag to plot intermediate step-by-step exploration figures.
%   :type other_figs: logical
%
%   .. note::
%      The algorithm generates a dual-panel publication figure combining the **Bottom-Up Selection Matrix** 
%      (showing sequential neuron inclusions with ticks and color-coded performances) 
%      and the **Performance Function** curve indicating the global optimum size :math:`k_{opt}`.
%
%   :Author: Laure WOLFF
%   :Date: May-June 2026

    % Initialization variableq
    [num_neurons,~,~] = size(spikes);
    best_order = zeros(1, num_neurons);    % Any neuron in the list
    neurons_dispo = 1:num_neurons;        % All neurons before the beginning of the algorithm 
    history_perf = zeros(1, num_neurons);  % The performance of each neuron at the step k
    Matrix_Grid = NaN(num_neurons, num_neurons);
    currentPop = [];
    
    %if showing, t_start = tic; end

    for k = 1:num_neurons
        num_dispo = length(neurons_dispo);
        current_step_perf = -Inf(1, num_dispo);
        
        for i = 1:num_dispo

            candidate = [currentPop neurons_dispo(i)];

            P_candidate = evaluate_population(spikes,candidate,Tmax,Distances,threshold);

            current_step_perf(i) = P_candidate;

            Matrix_Grid(k,neurons_dispo(i)) = P_candidate;

        end
        
        [best_perf_step, best_idx] = max(current_step_perf);
        best_neurone_step = neurons_dispo(best_idx);
        
        for i = 1:num_dispo
            Matrix_Grid(k, neurons_dispo(i)) = current_step_perf(i);
        end

        currentPop = [currentPop best_neurone_step];
        
        best_order(1:k) = currentPop;
        neurons_dispo(best_idx) = [];
        history_perf(k) = best_perf_step;
        
        if showing
            fprintf('Step k = %d | Adding neuron : %d | Performance P = %.4f\n', k, best_neurone_step, best_perf_step);  
        end
    end
    
    % if showing
    %     fprintf('--- Algorithm completed in %.2f seconds ---\n', toc(t_start));
    % end

    % To find the best subpopulation
    [~, idx_max_absolu] = max(history_perf);
    best_subpop = best_order(1:idx_max_absolu); 
    
    %% Block to improve the lisibility in the command windows
    if showing == true 
        for i = 1:10:length(history_perf)
            last_idx = min(i+9, length(history_perf));
            fprintf('   [%d-%d] : %s\n', i, last_idx, num2str(history_perf(i:last_idx), ' %.4f'));
        end
        
        fprintf('\nOptimal neuron inclusion order:\n');
        for i = 1:10:length(best_order)
            last_idx = min(i+9, length(best_order));
            fprintf('   %s\n', num2str(best_order(i:last_idx)));
        end
        
        fprintf('\nThe best subpopulation found contains %d neurons:\n', length(best_subpop));
        for i = 1:10:length(best_subpop)
            last_idx = min(i+9, length(best_subpop));
            fprintf('   %s\n', num2str(best_subpop(i:last_idx)));
        end
    end
    
    if plotting == true 
        if num_neurons <= 15
            tick_step = 1;      
        elseif num_neurons <= 40
            tick_step = 5;      
        else
            tick_step = 10;     
        end
        if other_figs == true
            %% 1. The plot
            [max_P, idx_max] = max(history_perf);
            
            neuron_labels = cell(1, num_neurons);
            for k = 1:num_neurons
                neuron_labels{k} = sprintf('N%d', best_order(k));
            end
            
            figure('Name', 'Bottom-Up optimization results');
            plot(1:num_neurons, history_perf, '-o', 'LineWidth', 2.5, 'Color', [0.30 0.58 0.20], ...
                 'MarkerEdgeColor', [0.30 0.58 0.20], 'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerSize', 8);
            hold on;
            
            line([idx_max, idx_max], [min(history_perf) - 0.05, max_P + 0.05], ...
                 'Color', [0.85 0.33 0.1], 'LineStyle', '--', 'LineWidth', 1.5);
             
            grid on; box on;
            xlim([0.5, num_neurons + 0.5]);
            ylim([min(history_perf) - 0.02, max(history_perf) + 0.04]);
            set(gca, 'XTick', 1:tick_step:num_neurons, 'YTick', 1:tick_step:num_neurons);
             
            xlabel('Neurons integrated sequentially (Step k)', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Global performance P', 'FontSize', 11, 'FontWeight', 'bold');
            title('Evolution of performance using Bottom-Up selection', 'FontSize', 12, 'FontWeight', 'bold');
            
            text(idx_max + 0.15, max_P, sprintf('Optimal subpopulation:\nNeurons: [%s]\nMax P = %.4f', ...
                 num2str(best_subpop), max_P), 'FontSize', 9, 'FontWeight', 'bold', ...
                 'BackgroundColor', [0.96 0.96 0.96], 'EdgeColor', [0.7 0.7 0.7]);
            legend({'Performance P(k)', 'Optimal size threshold'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
            hold off;
            shg;
            
            %% 2. Matrix
            figure('Name', 'Bottom-Up selection matrix');
            chosen_background_color = [1 1 1]; 
            set(gca, 'Color', chosen_background_color); 
            imagesc(Matrix_Grid, 'AlphaData', ~isnan(Matrix_Grid)); 
            colormap(jet); 
            colorbar;
            hold on;
            
            for step = 1:num_neurons
                chosen_neuron = best_order(step);
                plot(chosen_neuron, step, 'kx', 'MarkerSize', 12, 'LineWidth', 2.5);
            end
            
            min_P = min(history_perf); max_P = max(history_perf);
            scaled_perf = 1 + (num_neurons - 1) * (history_perf - min_P) / (max_P - min_P);
            plot(scaled_perf, 1:num_neurons, '-r', 'LineWidth', 2.5);
            
            box on;
            set(gca, 'XAxisLocation', 'bottom', 'YDir', 'reverse'); 
            set(gca, 'XTick', 1:tick_step:num_neurons, 'YTick', 1:tick_step:num_neurons);
              
            xlabel('# Neuron', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Number of neurons (Step k)', 'FontSize', 11, 'FontWeight', 'bold');
            title('Bottom-Up selection matrix', 'FontSize', 12, 'FontWeight', 'bold');
            
            legend({'Selected neuron (\times)', 'Max performance P'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
            hold off;
            shg;
        end

        %% 3. The paper figure
        figure('Name', 'Bottom-Up selection figure');

        opt_size = length(best_subpop); 
        min_perf_val = min(history_perf) - 0.02; 

        Matrix_Paper = Matrix_Grid;
        for k = 1:num_neurons
            past_neurons = best_order(1:k-1);
            Matrix_Paper(k, past_neurons) = min_perf_val;
            Matrix_Paper(k, best_order(k)) = history_perf(k);
        end

        subplot(1, 5, 1:3); 
        imagesc(1:num_neurons, 1:num_neurons, Matrix_Paper);
        colormap(jet);
        clim([min_perf_val, max(history_perf)+0.02]); 
        set(gca, 'YDir', 'normal'); 
        hold on;

        for i = 0.5 : 1 : num_neurons+0.5
            line([0.5, num_neurons+0.5], [i, i], 'Color', [1 1 1 0.2], 'LineWidth', 0.5);
            line([i, i], [0.5, num_neurons+0.5], 'Color', [1 1 1 0.2], 'LineWidth', 0.5);
        end

        for k = 1:num_neurons
            n_id = best_order(k);
            text(n_id, k, char(10003), 'FontSize', 11, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'FontWeight', 'bold', 'Color', [0 0 0]);

            if k < num_neurons
                plot(repmat(n_id, 1, num_neurons-k), (k+1):num_neurons, '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 5);
            end
        end

        for k_sub = 1:opt_size
            curr_n = best_order(k_sub);
            plot(curr_n, opt_size, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
        end

        rectangle('Position', [0.53, opt_size-0.45, num_neurons-0.06, 0.9], ...
                  'EdgeColor', [0.15 0.62 0.15], 'LineWidth', 2);

        box on;
        set(gca, 'TickDir', 'out', 'LineWidth', 1.1);
        set(gca, 'XTick', 1:tick_step:num_neurons, 'YTick', 1:tick_step:num_neurons);
        set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        xlabel('Neuron ID', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Size of population (k)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Bottom-Up algorithm matrix', 'FontSize', 12, 'FontWeight', 'bold');

        pos_mat = get(gca, 'Position');
        cb = colorbar('Position', [pos_mat(1) + pos_mat(3) + 0.02, pos_mat(2), 0.015, pos_mat(4)]);
        set(cb, 'LineWidth', 1.1, 'FontWeight', 'bold');
        ylabel(cb, 'Global Performance P', 'FontSize', 11, 'FontWeight', 'bold');

        subplot(1, 5, 4:5); 
        plot(history_perf, 1:num_neurons, '-ko', 'LineWidth', 2, ...
             'MarkerFaceColor', [0 0 0], 'MarkerSize', 5);
        hold on;

        plot(history_perf(opt_size), opt_size, 'ro', 'MarkerSize', 11, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
        plot(history_perf(opt_size), opt_size, 'rx', 'MarkerSize', 7, 'LineWidth', 1.5);

        box on; grid on;
        set(gca, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 1.1);
        set(gca, 'YTick', 1:tick_step:num_neurons, 'YTickLabel', []); 
        set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        ylim([0.5, num_neurons + 0.5]);
        xlim([min_perf_val+0.02, max(history_perf)+0.03]);
        xlabel('Best performance P', 'FontSize', 11, 'FontWeight', 'bold');
        title('Performance function', 'FontSize', 12, 'FontWeight', 'bold');

        hold off; 
        shg;

        % %% 3. The paper figure with the rectangle and the call of the BF also using with very high moderation (very high complexity)
        % figure('Name', 'Bottom-Up selection figure', 'Color', 'w', 'Position', [100, 100, 850, 500]);
        % 
        % opt_size = length(best_subpop);
        % min_perf_val = min(history_perf);
        % max_perf_val = max(history_perf);
        % 
        % % The matrix
        % Matrix_Paper = Matrix_Grid;
        % 
        % ax_matrix = subplot(1, 5, 1:3);
        % imagesc(1:num_neurons, 1:num_neurons, Matrix_Paper);
        % colormap(ax_matrix, jet(256));
        % clim(ax_matrix, [min_perf_val - 0.01, max_perf_val + 0.01]); 
        % set(gca, 'YDir', 'normal');
        % hold on;
        % 
        % for i = 0.5 : 1 : num_neurons+0.5
        %     line([0.5, num_neurons+0.5], [i, i], 'Color', [1 1 1 0.2], 'LineWidth', 0.5);
        %     line([i, i], [0.5, num_neurons+0.5], 'Color', [1 1 1 0.2], 'LineWidth', 0.5);
        % end
        % 
        % for k = 1:num_neurons
        %     n_id = best_order(k);
        %     text(n_id, k, char(10003), 'FontSize', 11, ...
        %          'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        %          'FontWeight', 'bold', 'Color', [0 0 0]);
        % 
        %     if k < num_neurons
        %         plot(repmat(n_id, 1, num_neurons-k), (k+1):num_neurons, '.', 'Color', [0.3 0.3 0.3], 'MarkerSize', 5);
        %     end
        % end
        % 
        % % Red cross for the subpopulations found by the bottom-up
        % % algorithms
        % for idx = 1:length(best_subpop)
        %     neuron_id = best_subpop(idx);
        %     plot(ax_matrix, neuron_id, opt_size, 'rx', 'MarkerSize', 12, 'LineWidth', 2.5);
        % end
        % 
        % % The rectagle which shox the real best subpopulation found by the
        % % BF algorithms
        % [idx_bf_neurons,~] = f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, false,false); 
        % opt_bf_size = length(idx_bf_neurons); 
        % 
        % x_start = min(idx_bf_neurons) - 0.45; 
        % y_start = opt_bf_size - 0.45;
        % width_rect = length(idx_bf_neurons) - 0.1;
        % height_rect = 0.9;
        % 
        % rectangle(ax_matrix, 'Position', [x_start, y_start, width_rect, height_rect], ...
        %           'EdgeColor', [0.15 0.65 0.15], 'LineWidth', 2.5, 'Curvature', [0 0]);
        % 
        % 
        % box on;
        % set(gca, 'TickDir', 'out', 'LineWidth', 1.1);
        % set(gca, 'XTick', 1:num_neurons, 'YTick', 1:num_neurons);
        % set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        % xlabel('Neuron ID', 'FontSize', 11, 'FontWeight', 'bold');
        % ylabel('Size of population (k)', 'FontSize', 11, 'FontWeight', 'bold');
        % title('Bottom-Up algorithm matrix', 'FontSize', 12, 'FontWeight', 'bold');
        % 
        % pos_mat = get(gca, 'Position');
        % cb = colorbar('Position', [pos_mat(1) + pos_mat(3) + 0.015, pos_mat(2), 0.015, pos_mat(4)]);
        % set(cb, 'LineWidth', 1.1, 'FontWeight', 'bold');
        % ylabel(cb, 'Global Performance P', 'FontSize', 11, 'FontWeight', 'bold');
        % 
        % % the perfomance function
        % subplot(1, 5, 4:5);
        % plot(history_perf, 1:num_neurons, '-ko', 'LineWidth', 2, ...
        %      'MarkerFaceColor', [0 0 0], 'MarkerSize', 5);
        % hold on;
        % 
        % plot(history_perf(opt_size), opt_size, 'ro', 'MarkerSize', 11, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
        % plot(history_perf(opt_size), opt_size, 'rx', 'MarkerSize', 7, 'LineWidth', 1.5);
        % 
        % box on; grid on;
        % set(gca, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 1.1);
        % set(gca, 'YTick', 1:num_neurons, 'YTickLabel', []);
        % set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        % ylim([0.5, num_neurons + 0.5]);
        % xlim([min_perf_val-0.01, max(history_perf)+0.02]);
        % xlabel('Best performance P', 'FontSize', 11, 'FontWeight', 'bold');
        % title('Performance function', 'FontSize', 12, 'FontWeight', 'bold');
        % 
        % hold off;
        % shg;
    end
end