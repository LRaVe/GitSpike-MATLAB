%% ALGO Computation
% Author: Maxime BELTOISE & Laure WOLFF
% Date: May-July 2026


clear; close all; clc;

% Make sure subfolders containing functions and order helpers are available
thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    repoRoot = fileparts(thisFile);
    addpath(genpath(repoRoot));
end

%% =====================================================
%% SEED FOR RANDOMNESS
%% =====================================================

ind_seed = 493480;  % to keep track of a seed's index

change_seed = 0;
if change_seed
    ind_seed = randi([1 1000000],1,1);
end
rng(ind_seed);


%% =====================================================
%% PARAMETERS
%% =====================================================

MODE = 'SP';               % Choose SP or LL dataset

global useMex
useMex = true;             % a must have if you have a mex compiler


%% =====================================================
%% SUMMED POPULATION PARAMETERS
%% =====================================================

params.N = 7;              % Number total of neurons (codant/Indi/non-codant)

params.c = 3;               % Coll (codant neurons)

params.nIndi = 0;           % Indi 
params.indiJitter = 0.01;

params.S = 4;               % Number of stimuli
params.R = 5;               % Number of repetition

params.Tmax = 1;
params.rate = 10;           % Number of spikes per train

params.Distances = [1 0 0 0];  % Use the SPIKE-Distance Classic / RI / A / RIA according to the position of the 1
params.threshold = 'auto';


%% =====================================================
%% LABELED LINE PARAMETERS
%% =====================================================

% Comment and uncomment the following sections to choose your parameters

%% -----------------------------
% Structured mode where you choose the response matrix
%% -----------------------------

%% Parameters to reproduce something similar to Fig.4

paramsLL.mode = 'structured';

paramsLL.N = 4;
paramsLL.S = 4;
paramsLL.R = 5;

paramsLL.Tmax = 1;
paramsLL.meanRate = 20;

paramsLL.Distances = [1 0 0 0];  % Use the SPIKE-Distance Classic / RI / A / RIA according to the position of the 1
paramsLL.threshold = 'auto';

paramsLL.jitter = 0.01;
paramsLL.jitterIntensity = [1 1 1 1];

paramsLL.responseMatrix = [

1 1 0 0
0 0 1 1
1 0 1 0
0 1 0 1

];

paramsLL.sameResponse = [1 1 1 1];  % neuron respond to the stimulus in the same way or not


%% Parameters to reproduce something similar to Fig.9

% paramsLL.mode = 'structured';
% 
% paramsLL.N = 10;
% paramsLL.S = 8;
% paramsLL.R = 5;
% 
% paramsLL.Tmax = 1;
% paramsLL.meanRate = 20;
% 
% paramsLL.Distances = [1 0 0 0];  % Use the SPIKE-Distance Classic / RI / A / RIA according to the position of the 1
% paramsLL.threshold = 'auto';
% 
% paramsLL.jitter = 0.01;
% paramsLL.jitterIntensity   = [1 1/3 1.5 1/3 2 2 2 1/3 1/2 1.5];
% 
% paramsLL.responseMatrix = [
% 
% 1 1 1 1 0 0 0 0
% 1 1 1 1 0 0 0 0
% 1 1 0 0 0 0 0 0
% 0 0 1 1 0 0 0 0
% 0 0 1 0 0 0 0 0
% 0 0 0 1 0 0 0 0
% 0 0 0 0 1 1 0 0
% 0 0 0 0 1 1 0 0
% 0 0 0 0 0 0 0 0
% 0 0 0 0 0 0 1 1
% 
% ];
% 
% paramsLL.sameResponse = [1 0 1 0 1 1 0 0 1 1];  % neuron respond to the stimulus in the same way or not


%% -----------------------------
% Random mode
%% -----------------------------

% paramsLL.mode = 'random';
% 
% paramsLL.N = 10;
% paramsLL.S = 8;
% paramsLL.R = 10;
% 
% paramsLL.Tmax = 1;
% paramsLL.meanRate = 20;
% 
% paramsLL.Distances = [1 0 0 0];  % Use the SPIKE-Distance Classic / RI / A / RIA according to the position of the 1
% paramsLL.threshold = 'auto';
% 
% paramsLL.jitter = 0.01;
% paramsLL.jitterIntensity   = [1 1 1 1 1 1 1 1 1 1];
% 
% paramsLL.connectionProbability = 0.4;
% 
% paramsLL.sameResponse = [1 1 1 1 1 1 1 1 1 1];  % neuron respond to the stimulus in the same way or not



%% =====================================================
%% GENERATE DATASET
%% =====================================================

% To create a Summed population (SP) dataset
if strcmp(MODE,'SP')
    spikes = generate_SP_dataset(params);
end

% To create a Labeled line (LL) dataset
if strcmp(MODE,'LL')
    [spikes,responseMatrix] = generate_LL_dataset(paramsLL);
end

% import_spikes_from_txt('filename',S,R,N);
% spikes = ans;

if strcmp(MODE,'SP')

    %% =====================================================
    %% FIGURE 1 (PART OF THE DATASET)
    %% =====================================================
    
    plotParams.stimuli = [1 2];         % Index of stimuli shown
    plotParams.repetitions = [1 2];     % Index of repetition shown
    
    plotParams.showPooling = true;
    
    %% =====================================================
    %% FIGURE 1 DISPLPAY
    %% =====================================================
    
    plot_SP_figure(spikes,params,plotParams);
    
    %% =====================================================
    %% FIGURE 2 (PERFORMANCE MATRICES)
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
        Psingle(n) = evaluate_population(spikes,n,params.Tmax,params.Distances,params.threshold);
    end
    
    disp('Neuron numbers :')
    disp(params.N)
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
    %% BRUTE FORCE
    %% =========================================
    
    if (params.N<8 && ~useMex) || (params.N<14 && useMex)  % variable security
        fprintf('Starting Brute Force...\n');
        [bestPop,bestP] = f_brute_force_V2(spikes,params.Tmax,params.Distances,params.threshold,true,true);
    end
    
    
    %% =========================================
    %% BOTTOM-UP
    %% =========================================
    
    fprintf('Starting Bottom-Up...\n');

    f_bottom_up_V2(spikes,params.Tmax,params.Distances,params.threshold,true,true,false);

    
    %% =========================================
    %% TOP-DOWN
    %% =========================================
    
    fprintf('Starting Top-Down...\n');

    result = top_down_gradient(spikes,params.Tmax,params.Distances,params.threshold);

    disp('Best population found :')

    disp(result.bestPopulation)

    disp(['Best P = ' num2str(result.bestP)])

    plot_top_down_gradient(result,codingNeurons);

    
    %% =========================================
    %% PARAMETERS FOR SIMULATED ANNEALING
    %% =========================================
    
    paramsSA.N0 = 50;
    paramsSA.steps = 5 * params.N;              % number of tests per plateau
    paramsSA.coolingFactor = 0.9;               
    paramsSA.codingNeurons = codingNeurons;
    
    
    %% =========================================
    %% SIMULATED ANNEALING
    %% =========================================
    
    fprintf('Starting Simulated Annealing...\n');
    
    SA = simulated_annealing(spikes,params.Tmax,params.Distances,params.threshold,paramsSA);
    
    disp('SA best population :')
    disp(SA.bestPopulation)
    
    disp(['SA best P = ' num2str(SA.bestP)])
    
    plot_simulated_annealing(SA, codingNeurons);

end


if strcmp(MODE,'LL');

    %% =====================================================
    %% FIGURE 4 AND 9 DISPLAY
    %% =====================================================
    
    result_LL = evaluate_LL_population(spikes,paramsLL.Tmax,paramsLL.Distances,paramsLL.threshold);
    
    if paramsLL.N ~= 4
        plot_LL_results(spikes,result_LL);
    end

    if paramsLL.N == 4
        plot_LL_figure4(spikes,result_LL,paramsLL,paramsLL.Distances,paramsLL.threshold);
    end

end


