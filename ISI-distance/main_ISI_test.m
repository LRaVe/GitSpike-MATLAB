%% ISI-distance tests
% Author: Laure WOLFF
% Date: May 2026

clear; 
clc; 
close all;

%% Dataset
num_trains = 3;
spikes_trains=cell(1,num_trains);
spikes_trains{1} = [0 1 2 4 7]; 
spikes_trains{2} = [3 4 6 10];
spikes_trains{3} = [2 5];

tmin=0;
tmax=10;
threshold = 1000; 

showing=14;                % +2:Distance,+4:Profile,+8:Matrix
plotting=14;               % +4:Profile,+8:Matrix

%% Call of the function of f_ISI_distance
[maMatrice, I,moyenneDistance] = f_ISI_distance(spikes_trains, tmin, tmax, ...
    showing,plotting);

%% Call of the function of f_ISI_distance_adaptative

%ISI-Distance Classique (MRTS = 0)
fprintf('--- Running classic ISI-distance with the adaptative ISI-distance function---\n');
try
    [mat_classic, I_list_c, mean_c] = f_ISI_distance_adaptive_v1(spikes_trains, tmin, tmax, 0, showing, plotting);
catch ME
    fprintf('\n[INFO MAIN] Classic mode was successfully blocked by security:\n');
    fprintf('--> %s\n\n', ME.message);
    mat_classic = []; I_list_c = []; mean_c = NaN;
end

% ISI-Distance Adaptive (MRTS = 'auto')
fprintf('\n--- Running adaptive auto ISI-distance (auto) ---\n');
[mat_adapt, I_list_a, mean_a] = f_ISI_distance_adaptive_v1(spikes_trains, ...
    tmin, tmax, 'auto',showing, plotting);

% ISI-Distance Adaptive (MRTS = 1000)
fprintf('\n--- Running adaptive manual ISI-distance (1000) ---\n');
[mat_adapt_d, I_list_d, mean_d] = f_ISI_distance_adaptive_v1(spikes_trains, ...
    tmin, tmax, threshold,showing, plotting);

% % Comparaison of the resultats in the window
% fprintf('\nRESULTS COMPARISON:\n');
% fprintf('Classic nean:  %.4f\n', mean_c);
% fprintf('Adaptive auto mean: %.4f\n', mean_a);
% fprintf('Adaptive manual mean: %.4f\n', mean_d);


