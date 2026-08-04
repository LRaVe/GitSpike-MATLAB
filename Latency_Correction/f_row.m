%% Function to calculate the shifts with the row method
% Author: Agathe JULIEN
% Date: June 2026

function shifts=f_row(td_matrix,row)
    % This function gets the shifts for each spike train based on the time difference matrix and the specified anchor row.
    %
    % Args:
    %     td_matrix (matrix): A square matrix where the element at (i,j) represents the time difference between spike trains i and j.
    %     row (int): The index of the anchor row to which all other trains will be aligned. 
    %
    % Returns:
    %     shifts (array): An array where each element represents the shift for the corresponding spike train. 


    shifts=td_matrix(row,:);
end 
