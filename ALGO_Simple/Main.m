%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026

clear;
close all;

% Make sure subfolders containing functions and order helpers are available
thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end

%% =====================================================
%% SEED FOR RANDOMNESS
%% =====================================================

ind_seed = 580633;  %good visual

change_seed = 1;
if change_seed
    ind_seed = randi([1 1000000],1,1);
end;
rng(ind_seed);

%% =====================================================
%% PARAMETERS
%% =====================================================

params.N = 7;
params.c = 3;

params.S = 4;
params.R = 5;

params.Tmax = 200;
params.rate = 0.03;

params.Distances = [0 0 1 0];  %Use the SPIKE-Distance Classic / RI / A / RIA according to the position of the 1
params.threshold = 'auto';

%% =====================================================
%% GENERATE DATASET
%% =====================================================

spikes = generate_SP_dataset(params);

%% =====================================================
%% FIGURE 1 PARAMETERS
%% =====================================================

plotParams.stimuli = [1 2];
plotParams.repetitions = [1 2];

plotParams.showPooling = true;

%% =====================================================
%% FIGURE 1
%% =====================================================

plot_SP_figure(spikes,params,plotParams);

%% =====================================================
%% FIGURE 2
%% =====================================================

codingNeurons = 1:params.c;
nonCoding = params.c+1:params.N;
allNeurons = 1:params.N;

%% -----------------------------
%% coding
%% -----------------------------

[trialsC,labels] = ...
    build_trials(spikes,codingNeurons);

DC = compute_population_distance_matrix( ...
    trialsC,params.Tmax,params.Distances,params.threshold);

%% -----------------------------
%% non coding
%% -----------------------------

[trialsNC,~] = ...
    build_trials(spikes,nonCoding);

DNC = compute_population_distance_matrix( ...
    trialsNC,params.Tmax,params.Distances,params.threshold);

%% -----------------------------
%% all neurons
%% -----------------------------

[trialsAll,~] = ...
    build_trials(spikes,allNeurons);

DALL = compute_population_distance_matrix( ...
    trialsAll,params.Tmax,params.Distances,params.threshold);

%% =====================================================
%% FIGURE 2 DISPLAY
%% =====================================================

figure(2);

subplot(3,1,1);
plot_distance_matrix(DC,labels,'C');

subplot(3,1,2);
plot_distance_matrix(DNC,labels,'NC');

subplot(3,1,3);
plot_distance_matrix(DALL,labels,'All');

