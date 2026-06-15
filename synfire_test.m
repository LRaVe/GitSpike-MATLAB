addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

trains_synfire = f_synfire(0, 100, 6, 4, 2, 0, 0.5);
%swap the order of the trains to test the robustness of the measures
for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end

[sortedOrders, sortedTimes] = order_spikes(0, 100, trains_synfire);
disp('Sorted Orders:');
disp(sortedOrders);


f_plot_trains_with_correction(trains_synfire,2);



