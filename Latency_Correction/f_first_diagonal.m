%% Function to get the shifts using the first diagonal method
% Author : Lucas RAVELOARINORO
% Date : July 2026

function shifts = f_first_diagonal(time_difference_matrix, row)
    % Align every train to the anchor row by accumulating upper-diagonal offsets.
    % Trains before the anchor get negative shifts, trains after the anchor get positive shifts.
    % This function computes the shifts for each spike train based on the time difference matrix and the specified anchor row.
    %
    % It aligns every train to the anchor row by accumulating upper-diagonal offsets. Trains before the anchor get negative shifts, trains after the anchor get positive shifts.
    %
    % Args:
    %     time_difference_matrix (matrix): A square matrix where the element at (i,j) represents the time difference between spike trains i and j. 
    %     row (int): The index of the anchor row to which all other trains will be aligned. 
    %
    % Returns:
    %     shifts (array): An array where each element represents the shift for the corresponding spike train.

    
    n_trains = size(time_difference_matrix, 1);
    shifts = zeros(1, n_trains);

    for train_idx = row-1:-1:1
        shifts(train_idx) = shifts(train_idx + 1) - time_difference_matrix(train_idx, train_idx + 1);
    end

    for train_idx = row+1:n_trains
        shifts(train_idx) = shifts(train_idx - 1) + time_difference_matrix(train_idx - 1, train_idx);
    end
end