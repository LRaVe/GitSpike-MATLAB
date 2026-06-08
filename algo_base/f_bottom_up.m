%% Bottom-up algorithm script
% Date: May-June 2026
% Author : Laure WOLFF

function f_bottom_up(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, plotting)
    % Initialisation des variables
    best_order = zeros(1, num_neurons);    % Any neuron in the list
    neurons_dispo = 1:num_neurons;        % All neurons before the beginning of the algorithm 
    history_perf = zeros(1, num_neurons);  % The performance of each neuron at the step k
    Matrix_Grid = NaN(num_neurons, num_neurons);
    
    for k = 1:num_neurons
        best_perf_step = -Inf;
        best_neurone_step = -1;

        for i = 1:length(neurons_dispo)
            neuron_test = neurons_dispo(i);
            
            % Creation of a mask to evaluate the best perfomance
            selection = zeros(num_neurons, 1);
            if k > 1
                selection(best_order(1:k-1)) = 1;
            end
            selection(neuron_test) = 1; % adding the current testing neuron
            
            % Calcul of the perfomance with the Summed Population hypothesis
            [current_perf, ~] = calculate_integrated_P_optimized(CellMatrix, selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
            Matrix_Grid(k, neuron_test) = current_perf;

            if current_perf > best_perf_step
                best_perf_step = current_perf;
                best_neurone_step = neuron_test;
            end
        end
        best_order(k) = best_neurone_step;

        % if k < num_neurons
        %     Matrix_Grid(k+1:end, best_neurone_step) = best_perf_step;
        % end
        
        % Remove the neuron from the list of available neurons
        neurons_dispo(neurons_dispo == best_neurone_step) = [];
        
        history_perf(k) = best_perf_step;
        fprintf('Step k = %d | Adding neuron : %d | Performance P = %.4f\n', k, best_neurone_step, best_perf_step);  
    end

    % To find the best subpopulation
    stop_idx = 1;
    while stop_idx < num_neurons && history_perf(stop_idx + 1) > history_perf(stop_idx)
        stop_idx = stop_idx + 1;
    end
    best_subpop = best_order(1:stop_idx);

    if showing == true 
        disp('History of performances P :');
        disp(history_perf);
        disp('Optimal neuron inclusion order');
        disp(best_order);
        disp('The best subpopulation are the neurons :');
        disp(best_subpop);
    end
    if plotting == true 

        %% 1. the plot
        % Find the best performance
        [max_P, idx_max] = max(history_perf);
        
        % Creation of the labels
        neuron_labels = cell(1, num_neurons);
        for k = 1:num_neurons
            neuron_labels{k} = sprintf('N%d', best_order(k));
        end
        
        figure('Name', 'Bottom-Up optimization results', 'Position', [200, 200, 700, 480]);
        plot(1:num_neurons, history_perf, '-o', 'LineWidth', 2.5, 'Color', [0.30 0.58 0.20], ...
             'MarkerEdgeColor', [0.30 0.58 0.20], 'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerSize', 8);
        hold on;
        % Indicate the best 
        line([idx_max, idx_max], [min(history_perf) - 0.05, max_P + 0.05], ...
             'Color', [0.85 0.33 0.1], 'LineStyle', '--', 'LineWidth', 1.5);
         
        grid on; box on;
        xlim([0.5, num_neurons + 0.5]);
        ylim([min(history_perf) - 0.02, max(history_perf) + 0.04]);
        
        % --- FORCER L'AXE X A AFFICHER L'ORDRE DES NEURONES ---
        set(gca, 'XTick', 1:num_neurons, 'XTickLabel', neuron_labels, 'FontSize', 10, 'FontWeight', 'bold');
        
        xlabel('Neurons integrated sequentially (Step k)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Global performance P', 'FontSize', 11, 'FontWeight', 'bold');
        title('Evolution of performance using Bottom-Up selection', 'FontSize', 12, 'FontWeight', 'bold');
        
        % Text box to have all information
        text(idx_max + 0.15, max_P, sprintf('Optimal subpopulation:\nNeurons: [%s]\nMax P = %.4f', ...
             num2str(best_subpop), max_P), 'FontSize', 9, 'FontWeight', 'bold', ...
             'BackgroundColor', [0.96 0.96 0.96], 'EdgeColor', [0.7 0.7 0.7]);
        legend({'Performance P(k)', 'Optimal size threshold'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
        hold off;
        shg;

        %% 1. Matrix
        figure('Name', 'Bottom-Up selection matrix', 'Position', [250, 250, 600, 500]);

        % Plotting the perfomance matrix
        chosen_background_color = [1 1 1]; 
        set(gca, 'Color', chosen_background_color); 
        imagesc(Matrix_Grid, 'AlphaData', ~isnan(Matrix_Grid)); 
        colormap(jet); 
        colorbar;
        hold on;
        
        % Drawing the black cross to indicate the chosen neuron
        for step = 1:num_neurons
            chosen_neuron = best_order(step);
            plot(chosen_neuron, step, 'kx', 'MarkerSize', 12, 'LineWidth', 2.5);
        end
        
        % Adding the maximum discrimination performance curve (red line)
        % Scale history_perf values to fit the plot range from 1 to N
        min_P = min(history_perf); max_P = max(history_perf);
        scaled_perf = 1 + (num_neurons - 1) * (history_perf - min_P) / (max_P - min_P);
        plot(scaled_perf, 1:num_neurons, '-r', 'LineWidth', 2.5);
        
        box on;
        set(gca, 'XAxisLocation', 'bottom', 'YDir', 'reverse'); % Étape 1 en haut
        set(gca, 'XTick', 1:num_neurons, 'YTick', 1:num_neurons);
        
        xlabel('# Neuron', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Number of neurons (Step k)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Bottom-Up selection matrix', 'FontSize', 12, 'FontWeight', 'bold');
        
        legend({'Selected neuron (\times)', 'Max performance P'}, 'Location', 'southoutside', 'Orientation', 'horizontal');
        hold off;
        shg;
    end
end