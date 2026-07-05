%% Benchmark Script: Brute Force vs Bottom-up/ Top-down vs Simulated Annealing Evaluation Count
% Date: June 2026
% Author: Laure WOLFF
clear; clc; close all;

% Make sure subfolders containing functions and order helpers are available
thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end

% 1. Simulation Parameters
neuron_counts = [5, 10, 15]; % Sizes to test (X-axis)

params.S = 4;
params.R = 5;

params.Tmax = 1;
params.rate = 10;

params.Distances = [0 0 1 0];
params.threshold = 'auto';

eval_greedy = zeros(size(neuron_counts));
eval_sa     = zeros(size(neuron_counts));
eval_sa_u   = zeros(size(neuron_counts));
eval_brute  = zeros(size(neuron_counts));

fprintf('Starting Evaluation Count Benchmark...\n');

% 2. Benchmark Loop 
for idx = 1:length(neuron_counts)
    
    params.N = neuron_counts(idx);

    params.c = floor(params.N/2);
    params.nIndi = 0;
    params.indiJitter = 0;

    fprintf('\n---------------------------------\n');
    fprintf('Testing N = %d\n',params.N);

    %% fake dataset

    spikes = generate_SP_dataset(params);
    
    % Benchmark Algorithm 1: Bottom-up/ Top-down 
    fprintf('Evaluating Bottom-up/ Top-down...\n');
    eval_greedy(idx) = params.N*(params.N+1)/2;
    
    % Benchmark Algorithm 2: Simulated Annealing
    fprintf('Evaluating Simulated Annealing...\n');

    loop = 3;

    paramsSA.N0 = 50;
    paramsSA.steps = 5*params.N;
    paramsSA.coolingFactor = 0.9;
    paramsSA.codingNeurons = 1:params.c;
    
    iterations = zeros(1,loop);
    uniquePop  = zeros(1,loop);
    
    fprintf('Running Simulated Annealing...\n')
    
    for k = 1:loop
    
        SA = simulated_annealing(spikes,params.Tmax,params.Distances,params.threshold,paramsSA);
    
        iterations(k) = SA.iterations;
        uniquePop(k)  = SA.cacheMisses;
    
    end
    
    eval_sa(idx)   = mean(iterations);
    eval_sa_u(idx) = mean(uniquePop);
    
    
    % Benchmark Algorithm 3: Brute Force
    fprintf('Evaluating Brute Force...\n');
    eval_brute(idx) = 2^params.N-1;

end
fprintf('\nBenchmark completed successfully!\n');

% 3. Plotting the Complexity Curves
figure('Name', 'Algorithmic Complexity Benchmark', 'Color', [1 1 1], ...
    'Position', [200, 200, 800, 550]);

% Bottom-Up Curve (Blue Circle)
plot(neuron_counts, eval_greedy, '-o', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0 0.4470 0.7410], 'Color', [0 0.4470 0.7410], ...
     'DisplayName', 'Bottom-Up/Top-Down (Polynomial: N(N+1)/2)');
hold on;

% Simulated Annealing Curve (Orange Square)
plot(neuron_counts, eval_sa, '-s', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0.8500 0.3250 0.0980], 'Color', ...
     [0.8500 0.3250 0.0980], 'DisplayName', 'Simulated Annealing (Heuristic)');

% Simulated Annealing unique Curve (Red triangle)
plot(neuron_counts, eval_sa_u, '-^', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'r', 'Color', 'r', ...
     'DisplayName', 'Simulated Annealing (Heuristic and unique)');


% Brute Force Curve (Purple Diamond)
plot(neuron_counts, eval_brute, '-d', 'LineWidth', 2, 'MarkerSize', 6, ...
     'MarkerFaceColor', [0.4940 0.1840 0.5560], 'Color', ...
     [0.4940 0.1840 0.5560], 'DisplayName', 'Brute Force (Exponential: 2^N-1)');

grid on; 
box on;
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