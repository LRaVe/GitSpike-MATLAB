%% Function to compute pairwise spike-train ordering for two spike trains
% Author: Agathe JULIEN
% Date: May 2026

function [res1,res2]=pairwise_train_order(tmin,tmax,spikes,spike_ind1,spike_ind2)
    % This function computes the pairwise SPIKE-train-order values for two distinct spike trains by comparing their spike times and determining the relative ordering of spikes within a specified coincidence window.
    % 
    % For each spike in the first train, it checks for coincidences with spikes in the second train. If a coincidence is found, both trains receive the same ordering value based on their relative timing.
    %
    % Args:
    %     tmin (float): The minimum time for the analysis window.
    %     tmax (float): The maximum time for the analysis window.
    %     spikes (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %     spike_ind1 (int): The index of the first spike train in the spikes cell array.
    %     spike_ind2 (int): The index of the second spike train in the spikes cell array.
    %
    % Returns:
    %     res1 (array): An array containing the SPIKE-train-order values for the first spike train, where each element corresponds to a spike in that train.
    %     res2 (array): An array containing the SPIKE-train-order values for the second spike train, where each element corresponds to a spike in that train.

    
    tol = 1e-10;

    n=length(spikes);
    if spike_ind1>n || spike_ind2>n || spike_ind1<1 || spike_ind2<1
        error('Index out of bounds');
    end

    s1=spikes{spike_ind1};
    s2=spikes{spike_ind2};
    res1=zeros(1,length(s1));
    res2=zeros(1,length(s2));

    for i=1:length(s1)
        for j=1:length(s2)
            if abs(s1(i)-s2(j))<coincidence_window(tmin,tmax,spikes,spike_ind1,spike_ind2,i,j)
                % Both trains get the same ordering sign based on relative timing
                if abs(s1(i)-s2(j)) <= tol
                    signValue=0;
                elseif s1(i)<s2(j)
                    signValue=1;
                elseif s1(i)>s2(j)
                    signValue=-1;
                else
                    signValue=0;
                end
                res1(i)=signValue;
                res2(j)=signValue;
                break;
            end
        end
    end
end