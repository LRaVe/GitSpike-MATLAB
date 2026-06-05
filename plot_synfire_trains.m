function plot_synfire_trains(trains, title_str)
    % Plot synfire trains with train 1 at top and last train at bottom
    % trains: cell array where trains{i} contains spike times for train i
    % title_str: (optional) title for the figure
    
    figure;
    hold on;
    
    % Plot each train
    for i = 1:length(trains)
        % Invert so train 1 is at top, last train at bottom
        y_pos = length(trains) - i + 1;
        
        % Plot vertical line for each spike
        for j = 1:length(trains{i})
            plot([trains{i}(j), trains{i}(j)], [y_pos-0.5, y_pos+0.5], 'k');
        end
    end
    
    xlabel('Time (ms)');
    ylabel('Train Index');
    
    % Set title if provided, otherwise use default
    if nargin > 1
        title(title_str);
    else
        title('Generated Synfire Trains');
    end
    
    hold off;
end
