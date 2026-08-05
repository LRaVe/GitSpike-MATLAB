%% Function to compute spike-train-order values for multiple spike trains
% Author: Agathe JULIEN
% Date: May 2026

function [results,order_matrix]=order_trains(tmin,tmax,spikes)
    % This function computes the SPIKE-train-order values for multiple spike trains by calculating pairwise orderings between all distinct train pairs and aggregating the results.
    %
    % Args: 
    %     tmin (float): The minimum time for the analysis window.
    %     tmax (float): The maximum time for the analysis window.
    %     spikes (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %
    % Returns:
    %     results (cell array): A cell array where each cell contains the SPIKE-train-order values for the corresponding spike train, normalized by the number of other trains.
    %     order_matrix (matrix): A matrix where the element at (i,j) represents the mean SPIKE-train-order value of train i with respect to train j, and the element at (j,i) is the negative of that value.


    n=length(spikes);
    results=cell(n,1);
    order_matrix=zeros(n,n);

    % Initialize result vectors for each train
    for i=1:n
        results{i}=zeros(1,length(spikes{i}));
        order_matrix(i,i)=0;  % Diagonal is always 0 (train compared to itself)
    end

    % Compute pairwise orderings between all distinct train pairs
    for i=1:n-1
        for j=i+1:n
            [res_i,res_j]=pairwise_train_order(tmin,tmax,spikes,i,j);
            results{i}=results{i}+res_i;
            results{j}=results{j}+res_j;
            order_matrix(i,j)=mean(res_i);
            order_matrix(j,i)=-mean(res_i);
        end
    end

    % Normalize results by number of other trains (n-1)
    if n>1
        for i=1:n
            results{i}=results{i}/(n-1);
        end
    end
end