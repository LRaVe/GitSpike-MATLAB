% Author : Lucas Raveloarinoro
% Date : 2024-06-05

function plot_synfire_trains(trains, sortedOrders,sortedTimes,title_str,is_corrected,shifts,trains_corrected,sortedOrders_corrected,sortedTimes_corrected)
    % Plot synfire trains with train 1 at top and last train at bottom
    % trains: cell array where trains{i} contains spike times for train i
    % title_str: (optional) title for the figure
    % is_corrected: (optional) boolean indicating if the trains are corrected, if so we add arrow from the previous spike to the current spike
    % shifts: (optional) vector of shifts applied to each train, used for plotting arrows if is_corrected is true
    hold on;
    box on;
    colormap(gca, jet(256));
    clim([-1 1]);
    
    % Plot each train

    if nargin < 5 || ~is_corrected
        [spikeColors,flatRows]=values_to_colors(trains, sortedOrders, length(trains));
        for k=1:numel(sortedTimes)
            line([sortedTimes(k) sortedTimes(k)],[flatRows(k)-0.5 flatRows(k)+0.5], 'Color', spikeColors(k,:),'LineWidth',1.5);
        end
    else
        [spikeColors_corrected,flatRows_corrected]=values_to_colors(trains_corrected, sortedOrders_corrected, length(trains));
        for l=1:numel(sortedTimes_corrected)
            for k=1:numel(sortedTimes_corrected)
                if is_corrected
                    line([sortedTimes_corrected(k)  sortedTimes_corrected(k)],[flatRows_corrected(k)-0.5 flatRows_corrected(k)+0.5],'Color',spikeColors_corrected(k,:),'LineWidth',1.5);
                else
                    line([sortedTimes_corrected(k) sortedTimes_corrected(k)],[flatRows_corrected(k)-0.5 flatRows_corrected(k)+0.5],'Color',spikeColors(k,:),'LineWidth',1.5);
                end
            end
        end
    end
    
    if nargin > 4 && is_corrected
        chosen_event = round(length(trains_corrected{1})/2); % choose the middle train as the reference for arrows
        for k=1:length(trains_corrected)
            y_start = length(trains_corrected) - k + 1; % y position for the current train
            x_start = trains{k}(chosen_event); % x position for the current train's chosen event
            dx = - shifts(k); % shift for the current train
            quiver(x_start, y_start, dx, 0, 0, 'MaxHeadSize', 0.5, 'Color', 'k', 'LineWidth', 1.2);
        end
    end
            
    
    set(gca, 'YTick', 1:length(trains));
    set(gca, 'YTickLabel', arrayfun(@num2str, length(trains):-1:1, 'UniformOutput', false));
    xlabel('Time (ms)');
    ylabel('Train Index');
    ylim([0.5 length(trains) + 0.5]);
    colorbar;
    
    % Set title if provided, otherwise use default
    if nargin > 1
        title(title_str);
    else
        title('Generated Synfire Trains');
    end
    
    hold off;
end
