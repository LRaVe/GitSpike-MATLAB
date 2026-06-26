addpath(genpath('.'))

trains_synfire = f_synfire(0, 100, 6, 3, 0, 0, 0.3); 

for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end


f_plot_trains_with_correction(trains_synfire,3,'sim_ann');


% shifts = f_latency_correction_sim_ann(trains_synfire, 0, 100);

% corrected_trains = cell(1, length(trains_synfire));
% for i = 1:length(trains_synfire)
%     corrected_trains{i} = trains_synfire{i} - shifts(i);
% end
    

% fprintf('Shifts: %s\n', mat2str(shifts));

% fprintf('Corrected Trains:\n');
% for i = 1:length(corrected_trains)
%     fprintf('Train %d: %s\n', i, mat2str(corrected_trains{i}));
% end