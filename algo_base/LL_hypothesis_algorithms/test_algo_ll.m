%% Script to test my functions in the LL hypothesis (using as a main script)
% Date: May-June 2026
% Author : Laure WOLFF
clear; clc; close all;
% clear functions; 


% name_diary = 'simulation_results_LL.txt';
% diary(name_diary); 
% diary on; 

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
refrac = 0.02;         % "an absolute refractory period of 2 ms" paper 2018
base_rate = 35;         % Frequency of the coding neurons (Hz) 
jitter_std = 0.005;   % Force du mouvement (Écart-type: 0.5 ms)

% Metric selection
% metric_choice = 'ISI_ADAPTIVE';
% metric_choice = 'SPIKE_DISTANCE';

showing = true;
plotting = true;         % Boolean to plot or not the main graphics
other_figs = true;       % Boolean to plot or not auxiliary figures
rng(12);                 % To reproduce the script without new random values

%% 2. Dataset creation (Labeled line hypothesis)

CellMatrix = generate_and_plot_raster_ll(num_stimuli, num_repetitions, ...
    num_indi, num_neurons, t1, t2, base_rate, refrac, jitter_std,showing, plotting, other_figs);
save('LL_matlab_data.mat', 'CellMatrix', '-v7');

% load ("LL.mat")
% CellMatrix = CELL;
% CellMatrix = cellfun(@double, CellMatrix, 'UniformOutput', false);
% test_spike = CellMatrix{1, 1, 1}(1);
% disp(test_spike);

% load("LL.mat"); 
% [num_neurons, num_stimuli, num_repetitions] = size(CELL); 
% 
% % Cleaning of the dataset 
% CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
% for n = 1:num_neurons
%     for st = 1:num_stimuli
%         for rp = 1:num_repetitions
%             spikes = CELL{n, st, rp};
% 
%             if isstruct(spikes) || isobject(spikes)
%                 if isfield(spikes, 'times') || isprop(spikes, 'times')
%                     spikes = spikes.times;
%                 elseif isfield(spikes, 'spikes') || isprop(spikes, 'spikes')
%                     spikes = spikes.spikes;
%                 end
%             end
% 
%             % Using double type and not uint32
%             CellMatrix{n, st, rp} = double(spikes(:)');
%         end
%     end
% end

% data = load("LL_python_data.mat");
% CellMatrix = data.LL_python_data;

%% 3. Distance matrix visualization
All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, plotting);

% data = load('distances_python.mat');         
% All_Matrix_D = data.All_Matrix_D_py;

% %% Comparason between MATLAB and Python
% data_py = load('distances_python.mat');
% D_py = data_py.All_Matrix_D_py;
% 
% % 1. Calcul de la différence absolue
% diff_matrix = abs(All_Matrix_D - D_py);
% max_diff = max(diff_matrix(:));
% mean_diff = mean(diff_matrix(:));
% 
% fprintf('\n=== DIAGNOSTIC DE COMPARAISON ===\n');
% fprintf('Différence maximale sur une cellule : %f\n', max_diff);
% fprintf('Différence moyenne : %f\n', mean_diff);
% 
% % 2. Trouver où se trouve la plus grosse erreur
% [row, col] = find(diff_matrix == max_diff, 1);
% if ~isempty(row)
%     fprintf('Plus grand écart trouvé à l''indice global : (%d, %d)\n', row, col);
%     fprintf('Valeur MATLAB : %f | Valeur Python : %f\n', All_Matrix_D(row,col), D_py(row,col));
% end
% fprintf('=================================\n\n');
% 
% fprintf('Somme : %f\n', sum(All_Matrix_D, 'all'));
% fprintf('Moyenne : %f\n', mean(All_Matrix_D, 'all'));

All_Matrices_M = calculate_plot_matrix_M(All_Matrix_D, num_neurons,num_stimuli, num_repetitions,plotting);

[P_all_neurons, P_pop, M_max, opt_LL, PLL_total] = calculate_and_plot_performance_matrix(All_Matrices_M, All_Matrix_D, num_neurons, num_stimuli, num_repetitions, showing, plotting);

% profile off;
% diary off;
% profile viewer;