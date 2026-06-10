%% Script to test my functions (using as a main script)
% Date: May-June 2026
% Author : Laure WOLFF
clear; clc; close all;

%% Global parameters
num_stimuli = 4;         % S
num_repetitions = 5;     % R
num_neurons = 125;         % N
num_coding_neurons = 50;  % c
t1 = 0; t2 = 1;          % Time window

refrac = 0.002;  % "an absolute refractory period of 2 ms" paper 2018
base_rate= 30;   % Frequency of the coding neurons (Hz)
%metric_choice = 'ISI_ADAPTIVE'
metric_choice = 'SPIKE_DISTANCE';

showing = true;
plotting = true; % Boolean to plotting or not the graphics
other_figs = true; % Boolean to plotting or not other figures

rng(50); % To reproduce the script witout new values

%% 2. Creation of the dataset with summed population hypothesis
CellMatrix = generate_and_plot_raster(num_stimuli, num_repetitions, ...
    num_coding_neurons, num_neurons, t1, t2, base_rate, refrac, false, false);

%% 3. Plooting the 3 matrix
plot_and_compute_distance_matrix(CellMatrix, num_neurons, ...
    num_coding_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice);

%% 4. Bottom-up algorithm
f_bottom_up(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2, metric_choice, showing, plotting, other_figs)

%% 5. Annealing
f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, ...
    t2, metric_choice, showing, plotting)

