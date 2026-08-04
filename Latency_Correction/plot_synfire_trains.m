%% Function to plot synfire spike trains with optional latency correction
% Author : Lucas Raveloarinoro
% Date : 2024-06-05

function plot_synfire_trains(trains, sortedOrders, sortedTimes, title_str, is_corrected, shifts, trains_corrected, sortedOrders_corrected, sortedTimes_corrected, ~)
    % This function plots the synfire spike trains. It can also plot the corrected spike trains if latency correction is applied. 
    %
    % Args:
    %     trains (cell array): A cell array where each cell contains the spike times for a specific spike train. 
    %     sortedOrders (array): An array of all spike orders sorted according to the sorted spike times.
    %     sortedTimes (array): An array of all spike times sorted in ascending order.
    %     title_str (str): The title for the plot. If not provided, a default title will be used.
    %     is_corrected (bool): A boolean indicating whether the spike trains have been corrected for latency. If true, the corrected spike trains will be plotted.
    %     shifts (array): An array where each element represents the shift for the corresponding spike train. 
    %     trains_corrected (cell array): A cell array where each cell contains the corrected spike times for a specific spike train. 
    %     sortedOrders_corrected (array): An array of all spike orders for the corrected spike trains sorted according to the sorted corrected spike times.
    %     sortedTimes_corrected (array): An array of all corrected spike times sorted in ascending order.

    
    hold on;
    box on;
    colormap(gca, jet(256));
    clim([-1 1]);
    
    num_trains = length(trains);
    
    % Color mapping for original trains
    [spikeColors, flatRows] = values_to_colors(trains, sortedOrders, num_trains);
    
    if nargin < 5 || ~is_corrected
        % Plotting original spikes
        for k = 1:numel(sortedTimes)
            line([sortedTimes(k) sortedTimes(k)], [flatRows(k)-0.5 flatRows(k)+0.5], 'Color', spikeColors(k,:), 'LineWidth', 1.5);
        end
    else
        % Color mapping for corrected trains
        [spikeColors_corrected, flatRows_corrected] = values_to_colors(trains_corrected, sortedOrders_corrected, num_trains);
        
        alignedRows = zeros(size(flatRows_corrected)); 
        alignedColors = zeros(numel(sortedTimes_corrected), 3); 
        alignedXStart = zeros(size(sortedTimes_corrected)); 
        already_associated_spikes = false(size(sortedTimes));
        
        % Align corrected spikes with original spikes
        for k = 1:numel(sortedTimes_corrected)
            current_row_corr = flatRows_corrected(k);
            current_time_corr = sortedTimes_corrected(k);
            
            candidate_indices = find(flatRows == current_row_corr & ~already_associated_spikes);
            
            if ~isempty(candidate_indices)
                [~, min_id] = min(abs(sortedTimes(candidate_indices) - current_time_corr));
                original_index = candidate_indices(min_id);
                
                already_associated_spikes(original_index) = true;
                alignedRows(k) = flatRows(original_index);
                alignedXStart(k) = sortedTimes(original_index);
                alignedColors(k, :) = spikeColors_corrected(k, :);
            else
                alignedRows(k) = flatRows_corrected(k);
                alignedColors(k, :) = spikeColors_corrected(k, :);
                alignedXStart(k) = current_time_corr;
            end
        end
        
        % Plot the corrected spikes
        for k = 1:numel(sortedTimes_corrected)
            line([sortedTimes_corrected(k) sortedTimes_corrected(k)], [alignedRows(k)-0.5 alignedRows(k)+0.5], 'Color', alignedColors(k,:), 'LineWidth', 1.5);
        end
        
        % Draw arrows from original spikes to corrected spikes
        if nargin > 5 && ~isempty(shifts)
            for k = 1:numel(sortedTimes_corrected)
                x_start = alignedXStart(k);         
                x_end = sortedTimes_corrected(k);   
                y_pos = alignedRows(k);             
                
                if abs(x_start - x_end) > 1e-5
                    line([x_start, x_end], [y_pos, y_pos], 'Color', 'k', 'LineWidth', 1.2);
                    orientation = sign(x_start - x_end);
                    distance = abs(x_start - x_end);
                    if distance > 1
                        size_x = 0.45 * orientation;
                    elseif 0.45 * distance < 0.15
                        size_x = 0.15 * orientation;
                    else
                        size_x = 0.45 * distance * orientation; 
                    end
                    size_y = 0.1;
                    line([x_end, x_end + size_x], [y_pos, y_pos + size_y], 'Color', 'k', 'LineWidth', 1.2);
                    line([x_end, x_end + size_x], [y_pos, y_pos - size_y], 'Color', 'k', 'LineWidth', 1.2);
                end
            end
        end
    end
    
    set(gca, 'YTick', 1:num_trains);
    set(gca, 'YTickLabel', arrayfun(@num2str, num_trains:-1:1, 'UniformOutput', false));
    xlabel('Time (ms)');
    ylabel('Train Index');
    ylim([0.5, num_trains + 0.5]);
    colorbar;
    
    if nargin >= 4 && ~isempty(title_str)
        title(title_str);
    else
        title('Generated Synfire Trains');
    end
    
    hold off;
end