%% Benchmark Script: Greedy vs Simulated Annealing Execution Time
% Date: June 2026
% Author: Laure WOLFF

clear; clc; close all;

% --- 1. Simulation Parameters ---
neuron_counts = [5, 10, 15, 20, 25]; % Sizes to test (X-axis)
num_stimuli = 4;
num_repetitions = 5;
t1 = 0; t2 = 4; % Time window
metric_choice = 'SPIKE_DISTANCE';

% Pre-allocate time vectors
time_greedy = zeros(size(neuron_counts));
time_sa = zeros(size(neuron_counts));

fprintf('Starting Benchmark...\n');

% --- 2. Benchmark Loop ---
for idx = 1:length(neuron_counts)
    N = neuron_counts(idx);
    fprintf('Testing with %d neurons...\n', N);
    
    % Generate dummy CellMatrix spike data for N neurons
    FakeCellMatrix = cell(N, num_stimuli, num_repetitions);
    for n = 1:N
        for s = 1:num_stimuli
            for r = 1:num_repetitions
                % Generate random spike times between t1 and t2
                num_spikes = randi([5, 15]);
                FakeCellMatrix{n,s,r} = sort(t1 + (t2 - t1) * rand(1, num_spikes));
            end
        end
    end
    
    % --- Benchmark Algorithm 1: Iterative / Greedy ---
    % Replace "f_greedy_algorithm" with your exact function name
    tic;
    f_bottom_up(FakeCellMatrix, N, num_stimuli, num_repetitions, t1, t2, metric_choice, false, false,false);
    % (Simulated placeholder for the test - remove the line below when un-commenting yours)
    pause(0.02 * N); 
    time_greedy(idx) = toc;
    
    % --- Benchmark Algorithm 2: Simulated Annealing ---
    % Showing and plotting set to false to avoid pop-ups during benchmark
    tic;
    f_simulated_annealing(FakeCellMatrix, N, num_stimuli, num_repetitions, t1, t2, metric_choice, false, false);
    time_sa(idx) = toc;
end

fprintf('Benchmark completed!\n');

% --- 3. Plotting the Performance Curves ---
figure('Name', 'Execution Time Benchmark', 'Color', [1 1 1], 'Position', [200, 200, 700, 500]);
plot(neuron_counts, time_greedy, '-o', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0 0.4470 0.7410], 'DisplayName', 'Iterative (Greedy) Algo');
hold on;
plot(neuron_counts, time_sa, '-s', 'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.8500 0.3250 0.0980], 'DisplayName', 'Simulated Annealing');

grid on; box on;
set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11);
xlabel('Number of Neurons in Pool', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Execution Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
title('Computational Complexity Analysis', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'NorthWest', 'FontSize', 11);

% Optional: Use log-scale if Simulated Annealing takes way more time
% set(gca, 'YScale', 'log'); 

hold off;