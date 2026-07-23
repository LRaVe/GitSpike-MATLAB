%% Script to test my functions in teh SP hypothesis (using as a main script)
% Date: May-June 2026
% Author : Laure WOLFF
clear; clc; close all;


name_diary = 'simulation_results.txt';
diary(name_diary); 
diary on; 

profile off;
profile('-historysize', 200000000000);

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
num_neurons = 7;        % N
num_coll = 3;           % Coding neurons initial (the information is in the sum of these neurons)
num_indi = 0;           % Individually coding neurons (they have some information each)
                          % --> crashes the algorithms written for SP hypothesis if there are any num_coll
t1 = 0; t2 = 1;         % Time window
refrac = 0.002;         % "an absolute refractory period of 2 ms" paper 2018
base_rate = 20;         % Frequency of the coding neurons (Hz)

% Metric selection
% metric_choice = 'ISI_ADAPTIVE';
%metric_choice = 'SPIKE_DISTANCE';

num_coding_neurons = num_indi + num_coll;
showing = true;
plotting = true;         % Boolean to plot or not the main graphics
other_figs = true;       % Boolean to plot or not auxiliary figures
rng(12);                 % To reproduce the script without new random values

%% 2. Dataset creation (Summed Population Hypothesis)

CellMatrix = generate_and_plot_raster(num_stimuli, num_repetitions, ...
    num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs);

% Created file to understand the structure of the data and use the txt
% format
f_export_simulation_to_txt(CellMatrix, 'Simulated_data.txt');

% % Created file ith the mat file
% save('Simulated_data.mat', 'CellMatrix');

% Alternative: Fail-case dataset for Bottom-Up testing
% To try to fail the BU algorithms (parameters : num_neurons = 10; num_coll = 4; num_indi= 3;)
% CellMatrix = generate_and_plot_raster_fail_BU(num_stimuli, num_repetitions, ...
%     num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs);

% % Test using a "simulated data" since a document file
% CellMatrix = f_import_data_secure('Simulated_data.txt', num_stimuli, num_repetitions, num_neurons);

% % Test using a "simulated data" since a MATLAB file
% load ("Simulated_data.mat");

% data = load("SP_python_data.mat");
% CellMatrix = data.SP_python_data;

%% 3. Distance matrix visualization
plot_and_compute_distance_matrix(CellMatrix, num_neurons, ...
    num_coding_neurons, num_stimuli, num_repetitions, t1, t2);

% 3.25 Individual performance matrix (Figure 7C)
P_individuelles = zeros(1, num_neurons);
for nc = 1:num_neurons
    selection_solo = zeros(1, num_neurons);
    selection_solo(nc) = 1;
    [P_solo, ~] = calculate_integrated_P_optimized(CellMatrix, selection_solo, ...
                    num_stimuli, num_repetitions, t1, t2);

    P_individuelles(nc) = P_solo;
end

figure('Name', 'Individual performances of the dataset', 'Color', 'w');
hBar = bar(P_individuelles, 'FaceColor', [0.30, 0.75, 0.93], 'EdgeColor', [0 0 0]);
hold on;

% Demarcation lines to separate the groups (Coll | Indi | NC)
 line([num_coll + 0.5, num_coll + 0.5], [0, max(P_individuelles)*1.2], 'Color', 'k', 'LineWidth', 1.5);
 line([num_coll + num_indi + 0.5, num_coll + num_indi + 0.5], [0, max(P_individuelles)*1.2], 'Color', 'k', 'LineWidth', 1.5);
 grid on; box on;
 xlim([0.5, num_neurons + 0.5]);
 ylim([0, max(P_individuelles)*1.2]);

% Adaptive tick spacing for the x-axis
 if num_neurons <= 15
     tick_step = 1;      
 elseif num_neurons <= 40
     tick_step = 5;     
 else
     tick_step = 10;    
 end
 set(gca, 'XTick', 1:tick_step:num_neurons);

% Subpopulation Labels
 xlabel('Neuron Index', 'FontSize', 11, 'FontWeight', 'bold');
 ylabel('Individual Performance', 'FontSize', 11, 'FontWeight', 'bold');
 title('Individual Performance Profile (Fig 7C)', 'FontSize', 12, 'FontWeight', 'bold');

 text(num_coll/2, max(P_individuelles)*1.1, 'Coll', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
 text(num_coll + num_indi/2, max(P_individuelles)*1.1, 'Indi', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
 text(num_coll + num_indi + (num_neurons - num_coll - num_indi)/2, max(P_individuelles)*1.1, 'NC', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
 hold off;

%%  4. Algorithms and performance profiling
profile on;

% 4.1 Brute Force Algorithm 
if num_neurons < 21 
    fprintf('\n--- Running Brute Force Optimization ---\n');
    f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
        t2, showing, other_figs);
end

% 4.2 Bottom-Up Algorithm
fprintf('\n--- Running Bottom-Up Optimization ---\n');
f_bottom_up(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2,showing, plotting, other_figs);

% 4.3 Simulated Annealing Algorithm
fprintf('\n--- Running Simulated Annealing Optimization ---\n');
f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2,showing, plotting, other_figs);

profile off;
diary off;
%profile viewer;
p = profile('info');
profsave(p, 'profile_results_html');
disp('Profiling completed successfully. Results saved in ./profile_results_html');
