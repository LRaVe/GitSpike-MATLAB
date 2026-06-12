addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

trains_synfire = f_synfire(0, 100, 6, 4, 3, 1, 0.1, 1); 
disp('Generated Synfire Trains:');
for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end

[sortedOrders, sortedTimes] = order_spikes(0, 100, trains_synfire);
disp('Sorted Orders:');
disp(sortedOrders);

% Plot the generated synfire trains
plot_synfire_trains(trains_synfire, sortedOrders, sortedTimes, 'Generated Synfire Trains');

td_matrix = f_TD_matrix(trains_synfire);
disp('Time Difference Matrix:');    
disp(td_matrix);
figure('Name', 'Time Difference Matrix', 'NumberTitle', 'off');
imagesc(td_matrix, 'CDataMapping', 'scaled');
colormap(gca, jet(256));
colorbar;
title('Time Difference Matrix');
[Cost_matrix,Cost_value] = f_Cost_matrix(trains_synfire);
disp('Cost Matrix:');
disp(Cost_matrix);
disp('Cost Value:');
disp(Cost_value);
figure('Name', 'Cost Matrix', 'NumberTitle', 'off');
imagesc(Cost_matrix, 'CDataMapping', 'scaled');
colormap(gca, jet(256));
colorbar;
title('Cost Matrix, Cost Value: ' + string(Cost_value));