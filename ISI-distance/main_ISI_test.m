%% ISI-distance tests
% Author: Laure WOLFF
% Date: May 2026
clear; 
clc; 
close all;

%% Dataset Configuration
num_trains = 3;
spikes_trains = cell(1, num_trains);
spikes_trains{1} = [0 1 2 4 7]; 
spikes_trains{2} = [3 4 6 10];
spikes_trains{3} = [2 5];

tmin = 0;
tmax = 10;
threshold = 1000; % Manual test value for MRTS

showing = 15; 
plotting = 15;

%% Classic ISI-distance 
fprintf('==================================================\n');
fprintf('        RUNNING CLASSIC ISI DISTANCE FUNCTION     \n');
fprintf('==================================================\n');

f_ISI_distance(spikes_trains, tmin, tmax, showing, plotting);


%% Adaptive/Classic ISI-distance
fprintf('\n==================================================\n');
fprintf('     RUNNING ADAPTIVE ISI DISTANCE FUNCTION       \n');
fprintf('==================================================\n');

% Classic ISI-Distance via the adaptive function (MRTS = 0)
fprintf('\n--- Running classic mode (Threshold = 0) ---\n');
f_ISI_distance_adaptive_v1(spikes_trains, tmin, tmax, 0, showing, plotting);

% Adaptive ISI-Distance (MRTS = 'auto')
fprintf('\n--- Running adaptive mode (Threshold = ''auto'') ---\n');
f_ISI_distance_adaptive_v1(spikes_trains, tmin, tmax, 'auto', showing, plotting);

% Adaptive ISI-Distance (Manual MRTS based on threshold)
fprintf('\n--- Running adaptive mode (Threshold = %.1f) ---\n', threshold);
f_ISI_distance_adaptive_v1(spikes_trains, tmin, tmax, threshold, showing, plotting);
