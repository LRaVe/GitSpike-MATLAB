%% Function to compute the overall spike-train-order value F
% Author: Agathe JULIEN
% Date: May 2026

function [F,sortedTimes,sortedOrders] = compute_spike_train_order_value(spikes, orders, number_spikes)
    % This function computes the overall SPIKE-train-order value F by concatenating all spike times and their corresponding order values, sorting them, and normlalizing by the total number of spikes. 
    %
    % Args:
    %     spikes (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %     orders (cell array): A cell array where each cell contains the SPIKE-order values for the corresponding spike times in the spikes cell array.
    %     number_spikes (int): The total number of spikes across all spike trains.
    %
    % Returns:
    %     F (float): The overall SPIKE-train-order value, normalized by the total number of spikes.
    %     sortedTimes (array): An array of all spike times sorted in ascending order.
    %     sortedOrders (array): An array of SPIKE-order values corresponding to the sorted spike times.

    
    time = horzcat(spikes{:});
    value = horzcat(orders{:});

    [sortedTimes,orderInd]=sort(time);
    sortedOrders=value(orderInd);
    F=sum(sortedOrders);

    % Normalize by number of spikes
    if number_spikes~=0
        F=F/number_spikes;
    end
end