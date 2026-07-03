% Author : Lucas Raveloarinoro
% Date : 2024-06-05

function plot_synfire_trains(trains, sortedOrders, sortedTimes, title_str, is_corrected, shifts, trains_corrected, sortedOrders_corrected, sortedTimes_corrected, cost_value)
    % Plot synfire trains with train 1 at top and last train at bottom
    
    hold on;
    box on;
    colormap(gca, jet(256));
    clim([-1 1]);
    
    num_trains = length(trains);
    cmap = jet(256);
    
    % Color mapping for original trains
    [~, flatRows] = values_to_colors(trains, sortedOrders, num_trains);
    spikeColors = zeros(numel(sortedTimes), 3);
    
    total_spikes = numel(sortedTimes);
    num_events = ceil(total_spikes / num_trains);
    
    for event = 1:num_events
        idx_start = (event - 1) * num_trains + 1;
        idx_end = min(event * num_trains, total_spikes);
        current_block_size = idx_end - idx_start + 1;
        
        if current_block_size > 1
            val_norm = -1 + 2 * (current_block_size - 1:-1:0) / (current_block_size - 1 + eps);
            id_color = 1 + round((val_norm - (-1)) * 255 / 2);
            id_color = max(1, min(256, id_color));
            spikeColors(idx_start:idx_end, :) = cmap(id_color, :);
        else
            spikeColors(idx_start:idx_end, :) = cmap(128, :); 
        end
    end
    
    if nargin < 5 || ~is_corrected
        % Plotting original spikes
        for k = 1:numel(sortedTimes)
            line([sortedTimes(k) sortedTimes(k)], [flatRows(k)-0.5 flatRows(k)+0.5], 'Color', spikeColors(k,:), 'LineWidth', 1.5);
        end
    else
        % Color mapping for corrected trains
        [~, flatRows_corrected] = values_to_colors(trains_corrected, sortedOrders_corrected, num_trains);
        spikeColors_corrected = zeros(numel(sortedTimes_corrected), 3);
        
        total_spikes_corr = numel(sortedTimes_corrected);
        num_events_corr = ceil(total_spikes_corr / num_trains);
        
        for event = 1:num_events_corr
            idx_start = (event - 1) * num_trains + 1;
            idx_end = min(event * num_trains, total_spikes_corr);
            current_block_size = idx_end - idx_start + 1;
            
            if current_block_size > 1
                val_norm_corr = -1 + 2 * (current_block_size - 1:-1:0) / (current_block_size - 1 + eps);
                id_color_corr = 1 + round((val_norm_corr - (-1)) * 255 / 2);
                id_color_corr = max(1, min(256, id_color_corr));
                spikeColors_corrected(idx_start:idx_end, :) = cmap(id_color_corr, :);
            else
                spikeColors_corrected(idx_start:idx_end, :) = cmap(128, :);
            end
        end
    
        if nargin < 10 || isempty(cost_value)
            aligned = false;
        else
            aligned = (cost_value == 0);
        end
        
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
            
                if aligned
                    alignedColors(k, :) = spikeColors(original_index, :);
                else
                    alignedColors(k, :) = spikeColors_corrected(k, :);
                end
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