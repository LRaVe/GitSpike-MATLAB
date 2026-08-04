%% SPIKE-order profile and matrix computation
% Author: Agathe JULIEN
% Date: May 2026

function [sortedOrders, sortedTimes, SO_matrix]=order_spikes(tmin,tmax,spikes)
    % This function computes the SPIKE-order profile and matrix for a given set of spike trains.
    %
    % Args:
    %     tmin (float): The minimum time for the analysis window. 
    %     tmax (float): The maximum time for the analysis window. 
    %     spikes (cell array): A cell array where each cell contains the spike times for a specific spike train. 
    %
    % Returns:
    %     sortedOrders (array): An array of all spike orders sorted according to the sorted spike times.
    %     sortedTimes (array): An array of all spike times sorted in ascending order.
    %     SO_matrix (array): A symmetric matrix where each element (i,j) represents the SPIKE-order value between spike train i and spike train j. 


    n=length(spikes);
    results=cell(n,1); 
    SO_matrix=zeros(n,n);

    if n==0
        return;
    end

    if n==1
        results{1}=zeros(1,length(spikes{1})); 
        return;
    end

    for i=1:n
        aggregated=zeros(1,length(spikes{i}));

        % Sum pairwise orderings with all other spike trains
        for j=1:n
            if i~=j
                pairwise=pairwise_order(tmin,tmax,spikes,i,j);
                aggregated=aggregated+pairwise;
            end
        end

        results{i}=aggregated/(n-1); 
    end

    time = horzcat(spikes{:}); 
    value = horzcat(results{:}); 

    [sortedTimes,orderInd]=sort(time);
    sortedOrders=value(orderInd);



    if nargout>1
        for i=1:n-1
            for j=i+1:n
                so_profile_i=pairwise_order(tmin,tmax,spikes,i,j);
                so_profile_j=pairwise_order(tmin,tmax,spikes,j,i);
                num_pair_spikes=length(spikes{i})+length(spikes{j});

                if num_pair_spikes>0
                    so_sum=sum(so_profile_i)+sum(so_profile_j);
                    so_value=so_sum/num_pair_spikes;
                    SO_matrix(i,j)=so_value;
                    SO_matrix(j,i)=so_value;
                end
            end
        end
    end
end
