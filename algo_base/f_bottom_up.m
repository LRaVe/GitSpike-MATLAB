%% Bottom-up algorithm script 
% Date: May-June 2026
% Author : Laure WOLFF
function f_bottom_up(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, plotting,other_figs)
    % Initialization variableq
    best_order = zeros(1, num_neurons);    % Any neuron in the list
    neurons_dispo = 1:num_neurons;        % All neurons before the beginning of the algorithm 
    history_perf = zeros(1, num_neurons);  % The performance of each neuron at the step k
    Matrix_Grid = NaN(num_neurons, num_neurons);
    
    %if showing, t_start = tic; end

    for k = 1:num_neurons
        num_dispo = length(neurons_dispo);
        current_step_perf = zeros(1, num_dispo);
        
        selection = zeros(num_neurons, 1);
        if k > 1
            selection(best_order(1:k-1)) = 1;
        end
        
        parfor i = 1:num_dispo
            neuron_test = neurons_dispo(i);
            local_selection = selection;
            local_selection(neuron_test) = 1;
            
            [perf, ~] = calculate_integrated_P_optimized(CellMatrix, local_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
            current_step_perf(i) = perf;
        end
        
        [best_perf_step, best_idx] = max(current_step_perf);
        best_neurone_step = neurons_dispo(best_idx);
        
        for i = 1:num_dispo
            Matrix_Grid(k, neurons_dispo(i)) = current_step_perf(i);
        end
        
        best_order(k) = best_neurone_step;
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
            set(gca, 'XTick', 1:num_neurons, 'XTickLabel', neuron_labels, 'FontSize', 10, 'FontWeight', 'bold');
            
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
            set(gca, 'XTick', 1:num_neurons, 'YTick', 1:num_neurons);
            
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
        set(gca, 'XTick', 1:num_neurons, 'YTick', 1:num_neurons);
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
        set(gca, 'YTick', 1:num_neurons, 'YTickLabel', []); 
        set(gca, 'FontSize', 10, 'FontWeight', 'bold');
        ylim([0.5, num_neurons + 0.5]);
        xlim([min_perf_val+0.02, max(history_perf)+0.03]);
        xlabel('Best performance P', 'FontSize', 11, 'FontWeight', 'bold');
        title('Performance function', 'FontSize', 12, 'FontWeight', 'bold');
        
        hold off; 
        shg;
    end
end