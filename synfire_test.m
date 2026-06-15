addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

trains_synfire = f_synfire(0, 100, 6, 5, 5, 0, 0.1); 
%swap the order of the trains to test the robustness of the measures
for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end

[sortedOrders, sortedTimes] = order_spikes(0, 100, trains_synfire);
disp('Sorted Orders:');
disp(sortedOrders);

%% Plot the synfire trains

% Plot the generated synfire trains
plot_synfire_trains(trains_synfire, sortedOrders, sortedTimes, 'Generated Synfire Trains');

%% Time Difference Matrix 

td_matrix = f_TD_matrix(trains_synfire);
disp('Time Difference Matrix:');    
disp(td_matrix);
figure('Name', 'Time Difference Matrix', 'NumberTitle', 'off');
imagesc(td_matrix, 'CDataMapping', 'scaled');
colormap(gca, jet(256));
colorbar;
title('Time Difference Matrix');

%% Cost Matrix

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


shifts = f_first_diagonal(td_matrix, 2);
disp('Shifts:');
disp(shifts);

%% Shifts with row method

shifts_row=f_row(td_matrix,1);
plot_shifts_row(trains_synfire,sortedOrders,shifts_row,sortedTimes);

f_plot_trains_with_correction(trains_synfire,2);



