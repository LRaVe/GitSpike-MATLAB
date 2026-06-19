addpath(genpath('.'))

trains_synfire = f_synfire(0, 100, 6, 1, 1, 0, 0.1); 
%swap the order of the trains to test the robustness of the measures
for i = 1:length(trains_synfire)
    fprintf('Train %d: %s\n', i, mat2str(trains_synfire{i}));
end


f_plot_trains_with_correction(trains_synfire,2,'row');



