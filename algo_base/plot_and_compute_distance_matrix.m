%% Spike distance matrix for coding, uncoding and all script 
% Date: May-June 2026
% Author : Laure WOLFF

function  plot_and_compute_distance_matrix(CellMatrix, num_neurons, num_coding_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice)
% PLOT_AND_COMPUTE_MATRICES Computes and plots the distance matrices for 
% Coding, Non-coding, and Full populations under the Summed Population hypothesis.
%
% Note: trial_labels are now automatically generated internally.

    num_trials = num_stimuli * num_repetitions;

    %% 1. Generation des labels en interne
    trial_labels = cell(1, num_trials);
    counter = 1;
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            trial_labels{counter} = sprintf('S%d-R%d', st, rp);
            counter = counter + 1;
        end
    end

    %% 2. Creation of the three selection masks
    coding_selection = [ones(num_coding_neurons, 1); zeros(num_neurons - num_coding_neurons, 1)];
    noise_selection  = [zeros(num_coding_neurons, 1); ones(num_neurons - num_coding_neurons, 1)];
    full_selection   = ones(num_neurons, 1);

    %% 3. Efficient calculation of the 3 Matrices
    [perf_A, Matrix_A] = calculate_integrated_P_optimized(CellMatrix, coding_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
    [perf_B, Matrix_B] = calculate_integrated_P_optimized(CellMatrix, noise_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
    [perf_C, Matrix_C] = calculate_integrated_P_optimized(CellMatrix, full_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);

    %% 4. Plotting the three distance matrices
    figure('Name', 'SP Distances matrix', 'Position', [25, 150, 1500, 450]);

    % Matrix A : Coding
    subplot(1, 3, 1);
    imagesc(Matrix_A); colormap('jet'); colorbar; axis square;
    set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
    set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none');
    xtickangle(45);
    title(sprintf('A. Coding subpopulation (C)\nP = %.4f', perf_A), 'Color', 'r', 'FontWeight', 'bold');
    xlabel('Trials'); ylabel('Trials');

    % Matrix B : Non-Coding 
    subplot(1, 3, 2);
    imagesc(Matrix_B); colormap('jet'); colorbar; axis square;
    set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
    set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none');
    xtickangle(45);
    title(sprintf('B. Non-coding subpopulation (NC)\nP = %.4f', perf_B), 'Color', 'b', 'FontWeight', 'bold');
    xlabel('Trials');

    % Matrix C : All
    subplot(1, 3, 3);
    imagesc(Matrix_C); colormap('jet'); colorbar; axis square;
    set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
    set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none');
    xtickangle(45);
    title(sprintf('C. Full population (All)\nP = %.4f', perf_C), 'Color', 'k', 'FontWeight', 'bold');
    xlabel('Trials');

    hold off; shg;
end