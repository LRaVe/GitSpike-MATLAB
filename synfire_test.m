addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

train_synfire = f_synfire(10, 5, 0, 100, 0, 0, 0); % Example usage: 10 trains, 5 spikes per train, time window 0-100 ms, jitter of 5 ms, shuffle 20% of trains, and skip 10% of spikes
disp('Generated Synfire Trains:');
for i = 1:length(train_synfire)
    fprintf('Train %d: %s\n', i, mat2str(train_synfire{i}));
end

[sortedOrders, sortedTimes, SO_matrix] = order_spikes(0, 100, train_synfire); % Compute SPIKE-order matrix for the generated synfire trains

num_trains = length(train_synfire);
plot_spikes = cell(1, num_trains);
tmin = 0; % Minimum time for plotting
tmax = 100; % Maximum time for plotting
for i = 1:num_trains
    plot_spikes{i} = train_synfire{i}; % Use the generated synfire trains for plotting
end

% Plot the generated synfire trains with colors based on individual spike order
figure('name', 'Synfire Trains', 'NumberTitle', 'off');
    hold on;
    set(gcf, 'Color','w');
    colormap(jet);  % Use jet colormap
    
    for i=1:num_trains
        y = num_trains - i + 1;
        for j=1:numel(plot_spikes{i})
            spike_time = plot_spikes{i}(j);
            % Find this spike in sortedTimes and get its order value
            [~, idx] = min(abs(sortedTimes - spike_time));
            if abs(sortedTimes(idx) - spike_time) < 1e-6  % Check if we found a close match
                order_val = sortedOrders(idx);
            else
                order_val = 0; % Default to middle of colormap
            end
            % Convert order value to RGB color using jet colormap
            % +1 should map to RED (high index), -1 should map to BLUE (low index)
            color = jet(256);
            color_idx = max(1, min(256, round(((order_val + 1) / 2) * 255 + 1)));
            spike_color = color(color_idx, :);
            
            line([spike_time spike_time], [y-0.35 y+0.35], 'Color', spike_color, 'LineWidth', 1.5);
        end
    end
    set(gca, 'Color', 'w');
    set(gca, 'XColor', 'k');
    set(gca, 'YColor', 'k');
    set(gca, 'YTick', []);
    xMargin = 0.001;
    xlim([tmin - xMargin, tmax + xMargin]);
    ylim([0.5 num_trains+0.5]);
    title('Synfire Trains', 'FontSize', 14, 'FontWeight', 'bold');
    
    % Add colorbar to show order value scale
    cbar = colorbar;
    cbar.Label.String = 'Spike Order Value';
    hold off;

    figure('name', 'SPIKE-Order Matrix', 'NumberTitle', 'off');
    imagesc(SO_matrix);
    colormap(jet);
    colorbar;
    title('SPIKE-Order Matrix', 'FontSize', 14, 'FontWeight', 'bold');
