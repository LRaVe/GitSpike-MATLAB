%% Script to test MEX files of the ALGO  project 
% Date: June 2026
% Author: Laure WOLFF
% Description: Compare the algorithms Brute Force and the simulated Annealing

clear; clc; close all;

setenv('MW_MINGW64_LOC', 'C:\mingw64') 

%% Dataset configuration
num_neurons = 20;          
num_stimuli = 3;           
num_repetitions = 4;       
num_coding_neurons = 16;   


t1 = 0;                    
t2 = 1;                    
base_rate = 30;            
refrac = 0.002;            

metric_choice = 'SPIKE_DISTANCE'; 

showing = true;            
plotting = true;           
other_fig = true; 

fprintf('==================================================\n');
fprintf('   Starting with %d neurons \n', num_neurons);
fprintf('==================================================\n\n');

%% Generate the dataset

CellMatrix = generate_and_plot_raster_mex_m(num_stimuli, num_repetitions, ...
    num_coding_neurons, num_neurons, t1, t2, base_rate, refrac, false, false);

fprintf(' Dataset was corretly generated \n');

%% Algorithm 1 : Force Brute (Brute Force)
fprintf('-> Force Brute algorithms (%d combinaisons)...\n', (2^num_neurons)-1);
tic;

% Call the MEX function
[best_subpop_bf, best_perf_bf, history_perf_brute] = f_brute_force_mex(...
    CellMatrix, ...
    double(num_neurons), ...
    double(num_stimuli), ...
    double(num_repetitions), ...
    double(t1), ...
    double(t2), ...
    metric_choice);
time_brute_force = toc;

if showing
    fprintf('\n================ BRUTE FORCE CONVERGED ================\n');
    fprintf('Best subpopulation : [%s]\n', num2str(best_subpop_bf));
    fprintf('Best performance P     = %.4f\n', best_perf_bf);
    fprintf('=======================================================\n');
end

%% Algorithm 2 : Simulated Annealing
fprintf('\n-> Simuated Annealing algorithm \n');
tic;

[nb_iterations] = f_simulated_annealing_mex_m(...
    CellMatrix, ...
    double(num_neurons), ...
    double(num_stimuli), ...
    double(num_repetitions), ...
    double(t1), ...
    double(t2), ...
    metric_choice, ...
    showing, ...
    plotting);
time_sa = toc;

[best_mask_sa, best_perf_sa, ~, ~, ~, ~, ~, ~, ~] = f_simulated_annealing_mex(...
    CellMatrix, ...
    double(num_neurons), ...
    double(num_stimuli), ...
    double(num_repetitions), ...
    double(t1), ...
    double(t2), ...
    metric_choice, ...
    false);

best_subpop_sa = find(best_mask_sa == 1);

%% Comparaison 
fprintf('\n==================================================\n');
fprintf('      Comparaison of the two algorithms             \n');
fprintf('==================================================\n');
fprintf('Paramètres : %d neurones | %d stimuli | %d rép\n\n', num_neurons, num_stimuli, num_repetitions);

fprintf('1. FORCE BRUTE :\n');
fprintf('   - Best subpopulation : [%s]\n', num2str(best_subpop_bf));
fprintf('   - Best performance P : %.4f\n', best_perf_bf);
fprintf('   - Running time       : %.4f secondes\n\n', time_brute_force);

fprintf('2. SIMULATED ANNEALING :\n');
fprintf('   - Best subpopulation : [%s]\n', num2str(best_subpop_sa));
fprintf('   - Best performance P : %.4f\n', best_perf_sa);
fprintf('   - Number of iterations : %d\n', nb_iterations);
fprintf('   - Running time       : %.4f secondes\n\n', time_sa);

%% Analysis and advantages of the simulated annealing algorithm
fprintf('3. ANALYSIS AND THE ADVANTAGE OF THE SIMUATED ANNEALING:\n');
if best_perf_bf > 0
    precision = (best_perf_sa / best_perf_bf) * 100;
    fprintf('   - Simulated Annealing Accuracy: %.2f%% of the absolute optimum\n', precision);
else
    precision = 100;
    fprintf('   - Simulated Annealing Accuracy: N/A (Brute Force Optimum = 0)\n');
end

gain_vitesse = time_brute_force / time_sa;
fprintf('   - Algorithm Speedup: x%.2f faster\n', gain_vitesse);

if precision >= 95
    fprintf('   - RESULT: [SUCCESS] Simulated annealing successfully converged to the optimum.\n');
else
    fprintf('   - RESULT: [WARNING] Discrepancy detected. Consider adjusting cooling_factor or iterations_per_temp.\n');
end
fprintf('==================================================\n');