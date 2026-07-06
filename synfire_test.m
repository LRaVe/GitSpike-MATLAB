addpath(genpath('.'))
tmin = 0;
tmax = 100;
n_trains = 6;
n_correct_events = 3;
n_random_events = 0;
n_inversed_events = 0;
overlap = 0.1; % 10% overlap between events
reference_train_idx = 3; % Use the third train as the reference (not used in simulated annealing)
mode = {'row','first_diagonal','sim_ann'}; % Modes to test


trains_synfire = f_synfire(tmin, tmax, n_trains, n_correct_events, n_random_events, n_inversed_events, overlap); 

for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end



f_plot_trains_with_correction(trains_synfire,reference_train_idx,'sim_ann',tmin,tmax);