%% Script to test my functions in the LL hypothesis (using as a main script)
% Date: May-June 2026
% Author : Laure WOLFF
clear; clc; close all;


name_diary = 'simulation_results_LL.txt';
diary(name_diary); 
diary on; 

profile off;
profile('-historysize', 2000000000);

thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end

% Command to find the C++ compiler and link with MATLAB
setenv('MW_MINGW64_LOC', 'C:\mingw64') 

%% 1. Global parameters

num_stimuli = 4;        % S
num_repetitions = 5;    % R
num_neurons = 4;        % N
num_indi = 4;           % Individually coding neurons (they have some information each)
                       
t1 = 0; t2 = 1;         % Time window
refrac = 0.002;         % "an absolute refractory period of 2 ms" paper 2018
base_rate = 20;         % Frequency of the coding neurons (Hz)

% Metric selection
% metric_choice = 'ISI_ADAPTIVE';
metric_choice = 'SPIKE_DISTANCE';

showing = true;
plotting = true;         % Boolean to plot or not the main graphics
other_figs = true;       % Boolean to plot or not auxiliary figures
rng(12);                 % To reproduce the script without new random values

%% 2. Dataset creation (Labeled line hypothesis)

CellMatrix = generate_and_plot_raster_ll(num_stimuli, num_repetitions, ...
    num_indi, num_neurons, t1, t2, base_rate, refrac, showing, plotting, other_figs);

%% 3. Distance matrix visualization
plot_and_compute_distance_matrix_ll(CellMatrix, num_indi, num_stimuli, num_repetitions, t1, t2)


% %%  4. Algorithms and performance profiling
% profile on;
% 
% % 4.1 Brute Force Algorithm 
% if num_neurons < 21 
%     fprintf('\n--- Running Brute Force Optimization ---\n');
%     f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
%         t2, metric_choice, showing, other_figs);
% end
% 
% % 4.2 Bottom-Up Algorithm
% fprintf('\n--- Running Bottom-Up Optimization ---\n');
% f_bottom_up(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
%     t2, metric_choice, showing, plotting, other_figs);
% 
% % 4.3 Simulated Annealing Algorithm
% fprintf('\n--- Running Simulated Annealing Optimization ---\n');
% f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
%     t2, metric_choice, showing, plotting, other_figs);
% 
% profile off;
% diary off;
% profile viewer;