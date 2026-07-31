%% Function to plot the shifted spikes with the row method 
% Author: Agathe JULIEN
% Date: June 2026


function plot_shifts_row(trains,sortedOrders,shifts,sortedTimes)
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
