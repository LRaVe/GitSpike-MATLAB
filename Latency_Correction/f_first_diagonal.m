

function shifts = f_first_diagonal(time_difference_matrix, row)
    % Align every train to the anchor row by accumulating upper-diagonal offsets.
    % Trains before the anchor get negative shifts, trains after the anchor get positive shifts.

    n_trains = size(time_difference_matrix, 1);
    shifts = zeros(1, n_trains);

    for train_idx = row-1:-1:1
        shifts(train_idx) = shifts(train_idx + 1) - time_difference_matrix(train_idx, train_idx + 1);
    end

    for train_idx = row+1:n_trains
        shifts(train_idx) = shifts(train_idx - 1) + time_difference_matrix(train_idx - 1, train_idx);
    end
end