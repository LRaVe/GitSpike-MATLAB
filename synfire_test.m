addpath(genpath('.'))
tmin = 0;
tmax = 100;
n_trains = 6; % Number of spike trains
% Choose the number of events for each type: correct, random, and inversed 
n_correct_events = 2; % Number of correctly ordered events (synfire spikes)
n_random_events = 3;  % Number of random events with jitter
n_inversed_events = 4; % Number of inverse-ordered events (reversed synfire spikes)
overlap = 0.1; % Overlap between events
reference_train_idx = 3; % Choose the train for the reference (won't be used in simulated annealing, but still needed for the function signature)
mode = {'row','first_diagonal','sim_ann'}; % Modes to test


trains_synfire = f_synfire(tmin, tmax, n_trains, n_correct_events, n_random_events, n_inversed_events, overlap); 

for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end



f_plot_trains_with_correction(trains_synfire,reference_train_idx,'sim_ann',tmin,tmax);