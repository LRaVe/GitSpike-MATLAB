%% Benchmark Script: Brute Force vs Greedy vs Simulated Annealing Evaluation Count
% Date: June 2026
% Author: Laure WOLFF
clear; clc; close all;

% --- 1. Simulation Parameters ---
neuron_counts = [5, 10, 15, 20, 25]; % Sizes to test (X-axis)
num_stimuli = 4;
num_repetitions = 5;
t1 = 0; t2 = 4; % Time window
metric_choice = 'SPIKE_DISTANCE';

% Pre-allocate evaluation counters
eval_greedy = zeros(size(neuron_counts));
eval_sa     = zeros(size(neuron_counts));
eval_sa_u     = zeros(size(neuron_counts));
eval_brute  = NaN(size(neuron_counts)); % NaN for missing points at N >= 20

fprintf('Starting Evaluation Count Benchmark...\n');

% --- 2. Benchmark Loop ---
for idx = 1:length(neuron_counts)
    N = neuron_counts(idx);
    fprintf('\n----------------------------------------\n');
    fprintf('Testing with %d neurons...\n', N);
    
    % --- Generate Artificial Dataset (Satuvuori 2018 Style) ---
    FakeCellMatrix = cell(N, num_stimuli, num_repetitions);
    for n = 1:N
        for s = 1:num_stimuli
            for r = 1:num_repetitions
                num_spikes = randi([5, 15]);
                FakeCellMatrix{n,s,r} = sort(t1 + (t2 - t1) * rand(1, num_spikes));
            end
        end
    end
    
    % --- Benchmark Algorithm 1: Bottom-Up ---
    fprintf('Evaluating Bottom-Up...\n');
    % Formule mathématique exacte du Bottom-Up codé dans f_bottom_up.m
    eval_greedy(idx) = (N * (N + 1)) / 2;
    
    % --- Benchmark Algorithm 2: Simulated Annealing ---
    fprintf('Evaluating Simulated Annealing...\n');
   nb_iterations_list = zeros(1, 10); 
    for i = 1:10
        [~, ~, nb_iter] = f_simulated_annealing(FakeCellMatrix, N, num_stimuli, num_repetitions, t1, t2, metric_choice, false, false);
        nb_iterations_list(i) = nb_iter;
    end
    eval_sa(idx) = mean(nb_iterations_list);

    % --- Benchmark Algorithm 2.5: Simulated Annealing unique ---
    fprintf('Evaluating Simulated Annealing Unique...\n');
    [nb_iter] = f_simulated_annealing(FakeCellMatrix, N, num_stimuli, num_repetitions, t1, t2, metric_choice, false, false);
    eval_sa_u(idx) = nb_iter;
    
    
    % --- Benchmark Algorithm 3: Brute Force) ---
    fprintf('Evaluating Brute Force...\n');
    % Nombre exact de masques binaires générés par dec2bin
    eval_brute(idx) = (2^N) - 1;

end
fprintf('\nBenchmark completed successfully!\n');

% --- 3. Plotting the Complexity Curves ---
figure('Name', 'Algorithmic Complexity Benchmark', 'Color', [1 1 1], 'Position', [200, 200, 800, 550]);

% Bottom-Up Curve (Blue Circle)
plot(neuron_counts, eval_greedy, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0 0.4470 0.7410], 'Color', [0 0.4470 0.7410], 'DisplayName', 'Bottom-Up/Top-Down (Polynomial: N(N+1)/2)');
hold on;

% Simulated Annealing Curve (Orange Square)
plot(neuron_counts, eval_sa, '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', [0.8500 0.3250 0.0980], 'DisplayName', 'Simulated Annealing (Heuristic)');

% Brute Force Curve (Purple Diamond - Stops at N=15)
plot(neuron_counts, eval_brute, '-d', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0.4940 0.1840 0.5560], 'Color', [0.4940 0.1840 0.5560], 'DisplayName', 'Brute Force (Exponential: 2^N-1)');

grid on; box on;
set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11);
set(gca, 'XTick', neuron_counts);

xlabel('Number of Neurons in Pool (N)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Number of Evaluated Subpopulations', 'FontSize', 12, 'FontWeight', 'bold');
title('Search Space Exploration Scale (Log Scale)', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'NorthWest', 'FontSize', 11);

% Keeping the log scale is essential to show the explosive behavior of 2^N
set(gca, 'YScale', 'log'); 

hold off;
shg;