%% Map sortedOrders to colors for spike times in spikes
% Author: Agathe JULIEN
% Date: June 2026

function [spikeColors,flatRows] = values_to_colors(spikes, sortedOrders,num_trains)
    % This function maps sortedOrders to colors and returns the corresponding color for each spike time in spikes.
    %
    % It also returns the corresponding row indices for each spike time for plotting. 
    % 
    % Args:
    %     spikes (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %     sortedOrders (array): An array of all spike orders sorted according to the sorted spike times.
    %     num_trains (int): The total number of spike trains.
    %
    % Returns:
    %     spikeColors (array): An array where each row corresponds to the RGB color for the corresponding spike time in spikes. 
    %     flatRows (array): An array where each element corresponds to the train index of the corresponding spike time in spikes. 


    flatRows = repelem(num_trains:-1:1, cellfun(@numel, spikes));
    [~, sortIdx] = sort([spikes{:}]);
    flatRows = flatRows(sortIdx);

    cmap = jet(256);
    sortedOrders = sortedOrders(:);
    vmin = -1;
    vmax = 1;

    spikeColors = repmat([0.7 0.7 0.7], length(sortedOrders), 1); % default for NaN / missing values of sortedOrders

    valid = ~isnan(sortedOrders); % only consider valid values for color mapping
    if ~any(valid)
        return
    end

    numColors = size(cmap, 1);
    idx = 1 + round((sortedOrders(valid) - vmin) * (numColors - 1) / (vmax - vmin));
    idx = max(1, min(numColors, idx));
    spikeColors(valid, :) = cmap(idx, :);
end
