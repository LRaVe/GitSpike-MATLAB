addpath('.');
addpath('SPIKE_order');
addpath('SPIKE_synchro');
addpath('spike_common');

tmin = 0;
tmax = 10;
num_trains = 5;
spikes = cell(1, num_trains);

% Synfire-style events shared by all trains.
% Set jitter = 0 for exact coincidences, or a small value to keep a color gradient.
event_times = [1.0 2.5 4.0 6.0 8.5];
jitter = 0.1;

for i = 1:num_trains
	if jitter == 0
		spikes{i} = event_times;
	else
		spikes{i} = event_times + (i - 1) * jitter;
	end
end

% order_spikes already returns the spike values aligned with the sorted times
[sortedOrders, sortedTimes] = order_spikes(tmin, tmax, spikes);
disp(sortedOrders);

% Flatten spike times once and build the matching train rows vectorized
flatTimes = [spikes{:}];
flatRows = repelem(num_trains:-1:1, cellfun(@numel, spikes));
[~, sortIdx] = sort(flatTimes);
flatRows = flatRows(sortIdx);

values = sortedOrders(:);
if isempty(values)
	error('No spike values returned by order_spikes.');
end

vmin = min(values);
vmax = max(values);
if vmin == vmax
	vmax = vmin + 1;
end

figure('Name', 'Colored Spike Trains', 'NumberTitle', 'off');
hold on;
box on;
set(gcf, 'Color', 'w');

cmap = jet(256);
colormap(cmap);
caxis([vmin vmax]);

numColors = size(cmap, 1);
colorIdx = 1 + round((values - vmin) * (numColors - 1) / (vmax - vmin));
colorIdx = max(1, min(numColors, colorIdx));
spikeColors = cmap(colorIdx, :);

%%
for k = 1:numel(sortedTimes)
	line([sortedTimes(k) sortedTimes(k)], [flatRows(k)-0.35 flatRows(k)+0.35], ...
		'Color', spikeColors(k, :), 'LineWidth', 1.5);
end

set(gca, 'Color', 'w');
set(gca, 'XColor', 'k');
set(gca, 'YColor', 'k');
set(gca, 'YTick', 1:num_trains);
set(gca, 'YTickLabel', arrayfun(@(i) sprintf('Spike train %d', i), num_trains:-1:1, 'UniformOutput', false));
xlabel('Time');
xlim([tmin tmax]);
ylim([0.5 num_trains + 0.5]);
colorbar;
title('Spike trains colored by spike-order value');
hold off;

%%
plot_synfire_trains(spikes, sortedOrders, sortedTimes, 'Synfire Trains Colored by Spike-Order Value');