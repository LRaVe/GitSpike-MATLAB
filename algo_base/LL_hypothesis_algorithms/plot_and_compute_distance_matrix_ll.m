%% Script to calculate the performance an plot the matrix in labeled line hypothesis 
% Date:June 2026
% Author : Laure WOLFF


function plot_and_compute_distance_matrix_ll(CellMatrix, num_indi, num_stimuli, num_repetitions, t1, t2, plotting)
    num_trials = num_stimuli * num_repetitions;
    
    %% 1. Allocation of the 3D Distance Array (T x T x N)
    % Stores one T x T matrix for each of the N neurons
    Dn_all_neurons = zeros(num_trials, num_trials, num_indi);
    Pn_all_neurons = zeros(num_indi, 1); % To store individual performances if needed
    
    %% 2. Loop over each individual neuron
    for n = 1:num_indi
        % Create a single-neuron mask (only neuron 'n' is active)
        single_neuron_selection = zeros(num_indi, 1);
        single_neuron_selection(n) = 1;
        
        % Compute SPIKE-distance matrix specifically for this single neuron
        % Note: Ensure your function returns the T x T distance matrix as the 2nd output
        [perf_n, Matrix_Dn] = calculate_labeled_line_P(CellMatrix, single_neuron_selection, num_stimuli, num_repetitions, t1, t2);
        
        % Store result in our 3D tensor
        Dn_all_neurons(:, :, n) = Matrix_Dn;
        Pn_all_neurons(n) = perf_n;
    end
    
    %% 4. Plotting the N individual distance matrices (Grid Layout)
    if plotting == true
            %% 1. Internal generation of trial labels (S1-R1, S1-R2...)
        trial_labels = cell(1, num_trials);
        counter = 1;
        for st = 1:num_stimuli
            for rp = 1:num_repetitions
                trial_labels{counter} = sprintf('S%d-R%d', st, rp);
                counter = counter + 1;
            end
        end
        
        figure('Name', 'LL Pairwise Distance Matrices Dn per Neuron', 'Color', 'w', 'Position', [50, 50, 1400, 900]);
        
        % Dynamically compute grid layout based on the number of neurons
        % e.g., if num_neurons = 20, layout will be 4 rows x 5 columns
        cols = ceil(sqrt(num_indi * 1.25)); 
        rows = ceil(num_indi / cols);
        
        for n = 1:num_indi
            subplot(rows, cols, n);
            
            % Display the T x T distance matrix for neuron 'n'
            imagesc(Dn_all_neurons(:, :, n)); 
            colormap('jet'); 
            axis square;
            
            % Colorbar can over-crowd small subplots, add it only if grid is small 
            if num_indi <= 6
                colorbar;
            end
            
            % Ticks management: only show full ticks if trials are low, otherwise keep it clean
            if num_trials <= 12
                set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
                set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'FontSize', 6, 'TickLabelInterpreter', 'none');
                xtickangle(45);
            else
                % For larger T (like T=20), just show numbers or keep it clear
                set(gca, 'XTick', [1, num_trials], 'YTick', [1, num_trials]);
                set(gca, 'XTickLabel', {'1', num2str(num_trials)}, 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 8);
            end
            
            % Title identification (Red text for your coding neurons, blue for non-coding)
            % Assuming your input CellMatrix has coding neurons first
            global_num_indi = evalin('caller', 'num_indi'); % Retrieves num_indi from caller script safely
            if n <= global_num_indi
                title(sprintf('Neuron %d (Coding)', n), 'Color', [0.85 0 0], 'FontSize', 9, 'FontWeight', 'bold');
            else
                title(sprintf('Neuron %d (Noise)', n), 'Color', [0 0.45 0.74], 'FontSize', 9, 'FontWeight', 'normal');
            end
            
            xlabel('Trials', 'FontSize', 7); 
            ylabel('Trials', 'FontSize', 7);
        end
        
        shg;
    end
end