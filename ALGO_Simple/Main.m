%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026

clear;
close all;
tic;
% Make sure subfolders containing functions and order helpers are available
thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end

%% =====================================================
%% SEED FOR RANDOMNESS
%% =====================================================

ind_seed = 288277;  %good visual

change_seed = 0;
if change_seed
    ind_seed = randi([1 1000000],1,1);
end
rng(ind_seed);

%% =====================================================
%% PARAMETERS
%% =====================================================

params.N = 10;              % Number of neurons (codant/Indi/non-codant)

params.c = 4;               % Coll (codant neurons)

params.nIndi = 0;           % Indi 
params.indiJitter = 0.2;

params.S = 4;               %Number of stimuli
params.R = 5;               %Number of repetition

params.Tmax = 200;

params.rate = 0.03;

params.Distances = [0 0 1 0];  %Use the SPIKE-Distance Classic / RI / A / RIA according to the position of the 1
params.threshold = 'auto';

%% =====================================================
%% GENERATE DATASET
%% =====================================================

spikes = generate_SP_dataset(params);

% load('dataset_topdown_fail.mat');



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

[trialsC,labels] = build_trials(spikes,codingNeurons);

DC = compute_population_distance_matrix(trialsC,params.Tmax,params.Distances,params.threshold);

%% -----------------------------
%% non coding
%% -----------------------------

[trialsNC,~] = build_trials(spikes,nonCoding);

DNC = compute_population_distance_matrix(trialsNC,params.Tmax,params.Distances,params.threshold);

%% -----------------------------
%% all neurons
%% -----------------------------

[trialsAll,~] = build_trials(spikes,allNeurons);

DALL = compute_population_distance_matrix(trialsAll,params.Tmax,params.Distances,params.threshold);



%% =====================================================
%% Single Performance (relevant in case of Indi neurons)
%% =====================================================

Psingle = zeros(1,params.N);

for n = 1:params.N

    Psingle(n) = evaluate_population( ...
        spikes,...
        n,...
        params.Tmax,...
        params.Distances,...
        params.threshold);



end

fprintf('Psingle : %.3f\n',Psingle);

figure(10);
bar(Psingle)
xlabel('Neuron')
ylabel('P')
title('Single neuron performance')
grid on



%% =====================================================
%% FIGURE 2 DISPLAY
%% =====================================================

figure(2);

subplot(1,3,1);
plot_distance_matrix(DC,labels,'C');

subplot(1,3,2);
plot_distance_matrix(DNC,labels,'NC');

subplot(1,3,3);
plot_distance_matrix(DALL,labels,'All');



%% =========================================
%% TOP-DOWN
%% =========================================

disp('Neuron numbers :')
disp(params.N)

result = top_down_gradient(spikes,params.Tmax,params.Distances,params.threshold);

disp('Best population found :')

disp(result.bestPopulation)

disp(['Best P = ' num2str(result.bestP)])

plot_top_down_gradient(result,codingNeurons);

toc;



%% =========================================
%% PARAMETERS FOR SIMULATED ANNEALING
%% =========================================

paramsSA.N0 = 50;
paramsSA.steps = 5 * params.N; %number of tests per plateau
paramsSA.coolingFactor = 0.95;
paramsSA.codingNeurons = codingNeurons;


%% =========================================
%% SIMULATED ANNEALING
%% =========================================

SA = simulated_annealing(spikes,params.Tmax,params.Distances,params.threshold,paramsSA);

disp('SA best population :')
disp(SA.bestPopulation)

disp(['SA best P = ' num2str(SA.bestP)])

plot_simulated_annealing(SA, codingNeurons);








toc;
