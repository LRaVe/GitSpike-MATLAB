addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

trains_synfire = f_synfire(0, 100, 6, 5, 1, 1, 0.1, 1); 
disp('Generated Synfire Trains:');
for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end

[sortedOrders, sortedTimes] = order_spikes(0, 100, trains_synfire);
disp('Sorted Orders:');
disp(sortedOrders);

% Plot the generated synfire trains
plot_synfire_trains(trains_synfire, sortedOrders, sortedTimes, 'Generated Synfire Trains');
