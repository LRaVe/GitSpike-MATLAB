addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

train_synfire = f_synfire(0, 100, 6, 9, 2, 4, 0.1, 1); 
disp('Generated Synfire Trains:');
for i = 1:length(train_synfire)
    fprintf('Train %d: %s\n', i, mat2str(train_synfire{i}));
end

% Plot the generated synfire trains
plot_synfire_trains(train_synfire, 'Generated Synfire Trains');
