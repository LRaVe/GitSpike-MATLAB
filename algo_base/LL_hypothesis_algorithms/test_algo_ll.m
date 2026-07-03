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
base_rate = 35;         % Frequency of the coding neurons (Hz) 

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
All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, plotting);

All_Matrices_M = calculate_plot_matrix_M(All_Matrix_D, num_neurons,num_stimuli, num_repetitions,plotting);

% profile off;
diary off;
% profile viewer;