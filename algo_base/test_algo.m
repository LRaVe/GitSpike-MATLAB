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
num_stimuli = 2;         % S
num_repetitions = 2;     % R
num_neurons = 20;        % N
num_coll = 10;            % coding neurons initial (the information is in the sum of these neurons)
num_indi= 0;            % individually coding neurons (they have some information each)
                         % --> crash the algorithms writtsen for SP
                         % hypothesis

t1 = 0; t2 = 1;          % Time window

refrac = 0.002;  % "an absolute refractory period of 2 ms" paper 2018
base_rate= 30;   % Frequency of the coding neurons (Hz)
%metric_choice = 'ISI_ADAPTIVE'
metric_choice = 'SPIKE_DISTANCE';

num_coding_neurons = num_indi + num_coll;

showing = true;
plotting = true; % Boolean to plotting or not the graphics
other_figs = true; % Boolean to plotting or not other figures

rng(50); % To reproduce the script witout new values

%% 2. Creation of the dataset with summed population hypothesis
CellMatrix = generate_and_plot_raster(num_stimuli, num_repetitions, ...
    num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs);

%% 3. Plooting the 3 matrix
plot_and_compute_distance_matrix(CellMatrix, num_neurons, ...
    num_coding_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);

%% 3.5 Brute force algorithm 
if num_neurons < 21 
    t_start = tic;
    f_brute_force(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2, metric_choice, showing, plotting);
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
    t2, metric_choice, showing, plotting);
fprintf ("Spent time is : %.4f", toc(t_start));
