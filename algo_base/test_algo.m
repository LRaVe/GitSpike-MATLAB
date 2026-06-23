%% Script to test my functions (using as a main script)
% Date: May-June 2026
% Author : Laure WOLFF
clear; clc; close all;

thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end

setenv('MW_MINGW64_LOC', 'C:\mingw64') % Command to find the C++ compiler and link with MatLab

%% Global parameters
num_stimuli = 2;        % S
num_repetitions = 3;    % R
num_neurons = 100;        % N
num_coll = 40;           % coding neurons initial (the information is in the sum of these neurons)
num_indi=10;            % individually coding neurons (they have some information each)
                         % --> crash the algorithms written for SP
                         % hypothesis if there are any num_coll

t1 = 0; t2 = 1;          % Time window

refrac = 0.002;  % "an absolute refractory period of 2 ms" paper 2018
base_rate= 20;   % Frequency of the coding neurons (Hz)
%metric_choice = 'ISI_ADAPTIVE'
metric_choice = 'SPIKE_DISTANCE';

num_coding_neurons = num_indi + num_coll;

showing = true;
plotting = true; % Boolean to plotting or not the graphics
other_figs = true; % Boolean to plotting or not other figures

rng(12); % To reproduce the script witout new values

%% 2. Creation of the dataset with summed population hypothesis
CellMatrix = generate_and_plot_raster(num_stimuli, num_repetitions, ...
    num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs);

% To try to fail the BU algorithms (parameters : num_neurons = 10; num_coll = 4; num_indi= 3;)
% CellMatrix = generate_and_plot_raster_fail_BU(num_stimuli, num_repetitions, ...
%     num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs);


%% 3. Plooting the 3 matrix
plot_and_compute_distance_matrix(CellMatrix, num_neurons, ...
    num_coding_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);

%% 3.25 Plooting the indivitual performance for classic dataset(Panneau C de la Fig. 7)
P_individuelles = zeros(1, num_neurons);

for nc = 1:num_neurons
    selection_solo = zeros(1, num_neurons);
    selection_solo(nc) = 1;
    
    [P_solo, ~] = calculate_integrated_P_optimized(CellMatrix, selection_solo, ...
                    num_stimuli, num_repetitions, t1, t2, 'SPIKE_DISTANCE');
                
    P_individuelles(nc) = P_solo;
end

figure('Name', 'Individual performances of the dataset', 'Color', 'w');
hBar = bar(P_individuelles, 'FaceColor', [0.30, 0.75, 0.93], 'EdgeColor', [0 0 0]);
hold on;

% Lines to separate the several groups (Coll | Indi | NC)
line([num_coll + 0.5, num_coll + 0.5], [0, max(P_individuelles)*1.2], 'Color', 'k', 'LineWidth', 1.5);
line([num_coll + num_indi + 0.5, num_coll + num_indi + 0.5], [0, max(P_individuelles)*1.2], 'Color', 'k', 'LineWidth', 1.5);

grid on; box on;
xlim([0.5, num_neurons + 0.5]);
ylim([0, max(P_individuelles)*1.2]);
if num_neurons <= 15
    tick_step = 1;      
elseif num_neurons <= 40
    tick_step = 5;     
else
    tick_step = 10;    
end
 set(gca, 'XTick', 1:tick_step:num_neurons);

% Labels
xlabel('Neuron Index', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Individual Performance', 'FontSize', 11, 'FontWeight', 'bold');
title('Individual Performance Profile (Fig 7C)', 'FontSize', 12, 'FontWeight', 'bold');

text(num_coll/2, max(P_individuelles)*1.1, 'Coll', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(num_coll + num_indi/2, max(P_individuelles)*1.1, 'Indi', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(num_coll + num_indi + (num_neurons - num_coll - num_indi)/2, max(P_individuelles)*1.1, 'NC', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
hold off;

%% 3.5 Brute force algorithm 
if num_neurons < 21 
    t_start = tic;
    f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2, metric_choice, showing, other_figs);
    fprintf ("Spent time is : %.4f \n", toc(t_start));
end

%% 4. Bottom-up algorithm
t_start = tic;
f_bottom_up(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2, metric_choice, showing, plotting, other_figs);
fprintf ("Spent time is : %.4f \n", toc(t_start));

%% 5. Annealing
t_start = tic;
f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2, metric_choice, showing, plotting,other_figs);
fprintf ("Spent time is : %.4f", toc(t_start));
