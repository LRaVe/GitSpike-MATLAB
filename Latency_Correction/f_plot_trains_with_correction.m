function [figures] = f_plot_trains_with_correction(trains,row,mode,t_min,t_max)
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

    if nargin < 1
        error('Please provide a cell array of trains as input.');
    end

    if row < 1 || row > length(trains)
        error('Row index must be between 1 and the number of trains.');
    end

    if nargin == 3 && ~ismember(mode, {'row', 'first_diagonal','sim_ann'})
        error('Mode must be either "row" or "first_diagonal" or "sim_ann".');
    end

    % Compute sorted orders and times for original trains
    [sortedOrders, sortedTimes] = order_spikes(t_min, t_max, trains);

    fprintf('(Debugging) Initial Sorted Orders:\n');
    disp(sortedOrders);
    
    % Compute time difference matrix and cost matrix
    td_matrix = f_TD_matrix(trains, t_min, t_max);
    [Cost_matrix, Cost_value] = f_Cost_matrix(trains, t_min, t_max);
    
    % Compute shifts based on first diagonal of time difference matrix
    if strcmp(mode, 'row')
        shifts = f_row(td_matrix, row);
    elseif strcmp(mode, 'first_diagonal')
        shifts = f_first_diagonal(td_matrix, row);
    elseif strcmp(mode, 'sim_ann')
        [shifts, costs] = f_lc_simulated_annealing(trains, t_min, t_max);
    end
    
    fprintf('(Debugging) Shifts values:\n');
    fprintf('%s\n', mat2str(shifts));

    % Apply correction to the original trains
    trains_corrected = cell(1, length(trains));
    for i = 1:length(trains)
        trains_corrected{i} = trains{i} - shifts(i);
    end
    
    % Compute sorted orders and times for corrected trains
    [sortedOrders_corrected, sortedTimes_corrected] = order_spikes(t_min, t_max, trains_corrected);
    %[~,~] = order_spikes(t_min, t_max, trains_corrected);
    fprintf('(Debugging) Sorted Orders Corrected:\n');
    disp(sortedOrders_corrected);
    
    td_matrix_corrected = f_TD_matrix(trains_corrected, t_min, t_max);
    [Cost_matrix_corrected, Cost_value_corrected] = f_Cost_matrix(trains_corrected, t_min, t_max);

    td_common_limits = [min([td_matrix(:); td_matrix_corrected(:)]), max([td_matrix(:); td_matrix_corrected(:)])];
    if td_common_limits(1) == td_common_limits(2)
        td_common_limits(2) = td_common_limits(2) + eps(td_common_limits(2) + 1);
    end

    cost_common_limits = [min([Cost_matrix(:); Cost_matrix_corrected(:)]), max([Cost_matrix(:); Cost_matrix_corrected(:)])];
    if cost_common_limits(1) == cost_common_limits(2)
        cost_common_limits(2) = cost_common_limits(2) + eps(cost_common_limits(2) + 1);
    end


    % Plotting the original and corrected trains, time difference matrices, and cost matrices, and if needed, the cost over iterations for simulated annealing
    figures = figure(1);
    set(gcf, 'Name', 'Synfire Trains and Corrections');

    if strcmp(mode, 'sim_ann')
        tl= tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    else
        tl = tiledlayout(2,3,'TileSpacing','Compact','Padding','Compact'); %#ok<*NASGU>
    end

    ax1 = nexttile(1);
    ax2 = nexttile(2);
    ax3 = nexttile(3);
    ax4 = nexttile(4);
    ax5 = nexttile(5);
    ax6 = nexttile(6);

    if strcmp(mode, 'sim_ann')
        ax7 = nexttile(7,[1 3]);
    end


    axes(ax1);
    plot_synfire_trains(trains, sortedOrders, sortedTimes, 'Original Synfire Trains');
    title('Original Synfire Trains');
    
    axes(ax2);
    imagesc(td_matrix, 'CDataMapping', 'scaled');
    colormap(gca, jet(256));
    caxis(td_common_limits);
    colorbar;
    title('Time Difference Matrix');

    axes(ax3);
    imagesc(Cost_matrix, 'CDataMapping', 'scaled');
    colormap(gca, jet(256));
    caxis(cost_common_limits);
    colorbar;
    title(['Cost Matrix, Cost Value: ' num2str(Cost_value)]);

    axes(ax4);
    plot_synfire_trains(trains, sortedOrders, sortedTimes, 'Corrected Synfire Trains',true,shifts,trains_corrected,sortedOrders_corrected,sortedTimes_corrected,Cost_value_corrected);
    %plot_shifts_row(trains,sortedOrders,shifts,sortedTimes);
    title('Corrected Synfire Trains');

    axes(ax5);
    imagesc(td_matrix_corrected, 'CDataMapping', 'scaled');
    colormap(gca,jet(256));
    caxis(td_common_limits);
    colorbar;
    title('New Time Difference Matrix');

    axes(ax6);
    imagesc(Cost_matrix_corrected, 'CDataMapping', 'scaled');
    colormap(gca,jet(256));
    caxis(cost_common_limits);
    colorbar;
    title(['New Cost Matrix, Cost Value: ' num2str(Cost_value_corrected)]);

    if strcmp(mode, 'sim_ann')
        axes(ax7);
        plot(costs);
        xlabel('Iteration');
        ylabel('Cost');
        title('Simulated Annealing Cost over Iterations');
    end

    fprintf('trains_corrected:\n');
    for i = 1:length(trains_corrected)
        fprintf('Train %d: %s\n', i, mat2str(trains_corrected{i}));
    end


end