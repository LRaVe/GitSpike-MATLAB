%% Function to plot the shifted spikes with the row method 
% Author: Agathe JULIEN
% Date: June 2026

function plot_shifts_row(trains,sortedOrders,shifts,sortedTimes)
    % This function plots the shifted spike trains based on the shifts calculated using the row method.
    %
    % Args:
    %     trains (cell array): A cell array where each cell contains the spike times for a specific spike train. 
    %     sortedOrders (array): An array of all spike orders sorted according to the sorted spike times.
    %     shifts (array): An array where each element represents the shift for the corresponding spike train. 
    %     sortedTimes (array): An array of all spike times sorted in ascending order. 

    
    sortedShifts=zeros(1,length(sortedTimes));
    j=1;
    for i=1:length(sortedTimes)
        if j~=length(shifts)
            sortedShifts(i)=sortedTimes(i)-shifts(j);
            j=j+1;
        else
            sortedShifts(i)=sortedTimes(i)-shifts(j);
            j=1;
        end
    end
    plot_synfire_trains(trains,sortedOrders,sortedShifts,'Shifts with Row Method',true,shifts);

end
