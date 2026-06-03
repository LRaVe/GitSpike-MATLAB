train_synfire = f_synfire(10, 5, 0, 100, 1, 0.2, 0.1); % Example usage: 10 trains, 5 spikes per train, time window 0-100 ms, jitter of 5 ms, shuffle 20% of trains, and skip 10% of spikes
disp('Generated Synfire Trains:');
for i = 1:length(train_synfire)
    fprintf('Train %d: %s\n', i, mat2str(train_synfire{i}));
end

num_trains = length(train_synfire);
plot_spikes = cell(1, num_trains);
tmin = 0; % Minimum time for plotting
tmax = 100; % Maximum time for plotting
for i = 1:num_trains
    plot_spikes{i} = train_synfire{i}; % Use the generated synfire trains for plotting
end
% Plot the generated synfire trains and their jittered versions

figure('name', 'Synfire Trains', 'NumberTitle', 'off');
    hold on;
    set(gcf, 'Color','w');
    for i=1:num_trains
        y = num_trains - i + 1;
        for j=1:numel(plot_spikes{i})
            line([plot_spikes{i}(j) plot_spikes{i}(j)], [y-0.35 y+0.35], 'Color', 'k', 'LineWidth', 1.5);
        end
    end
    set(gca, 'Color', 'w');
    set(gca, 'XColor', 'k');
    set(gca, 'YColor', 'k');
    set(gca, 'YTick', []);
    xMargin = 0.001;
    xlim([tmin - xMargin, tmax + xMargin]);
    ylim([0.5 num_trains+0.5]);
    title(figure_title, 'FontSize', 14, 'FontWeight', 'bold');
    hold off;