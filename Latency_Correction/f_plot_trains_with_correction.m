function [figures] = f_plot_trains_with_correction(trains,row,mode)
    % Generate synfire trains and apply correction based on time difference matrix
    % trains: cell array where trains{i} contains spike times for train i
    % row: the row of the time difference matrix to use for correction (default is 1)
    % mode: correction mode (default is 'row')
    
    if nargin < 3
        mode = 'row'; % Default to row-based correction
    end
    if nargin < 2
        row = 1; % Default to using the first row for correction
    end

    % Compute sorted orders and times for original trains
    [sortedOrders, sortedTimes] = order_spikes(0, 100, trains);
    
    % Compute time difference matrix and cost matrix
    td_matrix = f_TD_matrix(trains, 0, 100);
    [Cost_matrix, Cost_value] = f_Cost_matrix(trains, 0, 100);
    
    % Compute shifts based on first diagonal of time difference matrix
    if strcmp(mode, 'row')
        shifts = f_row(td_matrix, row);
    elseif strcmp(mode, 'first_diagonal')
        shifts = f_first_diagonal(td_matrix, row);
    end
    
    disp('Shifts');
    disp(shifts);

    % Apply correction to the original trains
    trains_corrected = cell(1, length(trains));
    for i = 1:length(trains)
        trains_corrected{i} = trains{i} - shifts(i);
    end
    
    % Compute sorted orders and times for corrected trains
    [sortedOrders_corrected, sortedTimes_corrected] = order_spikes(0, 100, trains_corrected);
    %[~,~] = order_spikes(0, 100, trains_corrected);
    
    td_matrix_corrected = f_TD_matrix(trains_corrected, 0, 100);
    [Cost_matrix_corrected, Cost_value_corrected] = f_Cost_matrix(trains_corrected, 0, 100);



    figures = figure(1);
    set(gcf, 'Name', 'Synfire Trains and Corrections');
    tl = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>

    ax1 = nexttile(1);
    ax2 = nexttile(2);
    ax3 = nexttile(3);
    ax4 = nexttile(4);
    ax5 = nexttile(5);
    ax6 = nexttile(6);

    axes(ax1);
    plot_synfire_trains(trains, sortedOrders, sortedTimes, 'Original Synfire Trains');
    title('Original Synfire Trains');
    
    axes(ax2);
    imagesc(td_matrix, 'CDataMapping', 'scaled');
    colormap(gca, jet(256));
    colorbar;
    title('Time Difference Matrix');

    axes(ax3);
    imagesc(Cost_matrix, 'CDataMapping', 'scaled');
    colormap(gca, jet(256));
    colorbar;
    title(['Cost Matrix, Cost Value: ' num2str(Cost_value)]);

    axes(ax4);
    plot_synfire_trains(trains, sortedOrders, sortedTimes, 'Corrected Synfire Trains',true,shifts,trains_corrected,sortedOrders_corrected,sortedTimes_corrected);
    %plot_shifts_row(trains,sortedOrders,shifts,sortedTimes);
    title('Corrected Synfire Trains');

    axes(ax5);
    imagesc(td_matrix_corrected, 'CDataMapping', 'scaled');
    colormap(gca,jet(256));
    colorbar;
    title('New Time Difference Matrix');

    axes(ax6);
    imagesc(Cost_matrix_corrected, 'CDataMapping', 'scaled');
    colormap(gca,jet(256));
    colorbar;
    title(['New Cost Matrix, Cost Value: ' num2str(Cost_value_corrected)]);
    
end