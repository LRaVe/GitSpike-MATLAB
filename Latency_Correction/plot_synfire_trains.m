% Author : Lucas Raveloarinoro
% Date : 2024-06-05

function plot_synfire_trains(trains, sortedOrders,sortedTimes,title_str,is_corrected,shifts)
    % Plot synfire trains with train 1 at top and last train at bottom
    % trains: cell array where trains{i} contains spike times for train i
    % title_str: (optional) title for the figure
    % is_corrected: (optional) boolean indicating if the trains are corrected, if so we add arrow from the previous spike to the current spike
    % shifts: (optional) vector of shifts applied to each train, used for plotting arrows if is_corrected is true
    hold on;
    box on;
    [spikeColors,flatRows]=values_to_colors(trains, sortedOrders, length(trains));
    colormap(gca, jet(256));
    clim([-1 1]);
    
    % Plot each train
    for k=1:numel(sortedTimes)
        line([sortedTimes(k) sortedTimes(k)],[flatRows(k)-0.5 flatRows(k)+0.5], 'Color', spikeColors(k,:),'LineWidth',1.5);
    end
    
    % If the trains are corrected, plot arrows from previous spike to current spike based on shifts
    
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
