%% ISI-distance tests
% Author: Laure WOLFF
% Date: May 2026
clear; 
clc; 
close all;

%% Dataset Configuration
num_trains = 3;
spikes_trains = cell(1, num_trains);
spikes_trains{1} = [1 2 4 7 10]; 
spikes_trains{2} = [0 3 4 6 ];
spikes_trains{3} = [2 5];

tmin = 0;
tmax = 10;
threshold = 1000; % Manual test value for MRTS

[spikes, ~, ~] = add_auxiliary_spikes_bis(spikes_trains, tmin, tmax);

showing = 15; 
plotting = 15;

%% Classic ISI-distance 
fprintf('==================================================\n');
fprintf('        RUNNING CLASSIC ISI DISTANCE FUNCTION     \n');
fprintf('==================================================\n');

f_ISI_distance(spikes, tmin, tmax, showing, plotting);

%% Adaptive/Classic ISI-distance
fprintf('\n==================================================\n');
fprintf('     RUNNING ADAPTIVE ISI DISTANCE FUNCTION       \n');
fprintf('==================================================\n');

fprintf('\n--- Running classic mode (Threshold = 0) ---\n');
f_ISI_distance_adaptive_v1(spikes, tmin, tmax, 0, showing, plotting);

fprintf('\n--- Running adaptive mode (Threshold = ''auto'') ---\n');
f_ISI_distance_adaptive_v1(spikes, tmin, tmax, 'auto', showing, plotting);

fprintf('\n--- Running adaptive mode (Threshold = %.1f) ---\n', threshold);
f_ISI_distance_adaptive_v1(spikes, tmin, tmax, threshold, showing, plotting);