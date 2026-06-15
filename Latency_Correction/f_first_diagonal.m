

function shifts = f_first_diagonal(time_difference_matrix, row)
    % Get the first diagonal of the time difference matrix for a given row
    % time_difference_matrix: square matrix of time differences
    % row: reference row for which to calculate the shifts
    
    shifts = zeros(1, size(time_difference_matrix, 2)); % Initialize shifts vector
    multipliers = zeros(1, size(time_difference_matrix, 1));
    for i=1:size(time_difference_matrix, 1)
        if i < row
            multipliers(i) = -1; % For rows before the reference row, the shift is -1
        else
            multipliers(i) = 1; % For rows after the reference row, the shift is +1
        end
    end
    sum_row = 0;
    for rows=1:(size(time_difference_matrix, 1))
        if rows ~= row
            if rows < row
                shifts(rows) = sum_row + multipliers(rows) * time_difference_matrix(rows, rows+1); % Extract the time difference for the given row and all columns
                sum_row = sum_row + time_difference_matrix(rows, rows+1);
            else
                shifts(rows) = sum_row + multipliers(rows) * time_difference_matrix(rows-1, rows); % Extract the time difference for the given row and all columns
                sum_row = sum_row + time_difference_matrix(rows-1, rows);
            end
        else
            shifts(rows) = 0; % For the given row, the shift is the sum of the previous costs
            sum_row = 0; 
        end
    end
end