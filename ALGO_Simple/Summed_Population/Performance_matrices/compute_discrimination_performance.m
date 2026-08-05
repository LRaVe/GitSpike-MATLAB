%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function P = compute_discrimination_performance(D,labels)
% COMPUTE_DISCRIMINATION_PERFORMANCE Calculates classification performance from distance matrices.
%
%   Computes the overall discrimination performance score :math:`P` by evaluating 
%   the difference between the mean inter-stimulus distance and the mean intra-stimulus 
%   distance across all paired trials:
%
%   .. math::
%
%      P = \langle D_{\text{inter}} \rangle - \langle D_{\text{intra}} \rangle
%
%   A higher performance value :math:`P` indicates distinct, well-separated 
%   population responses for different stimuli with high repeatability within the same stimulus class.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Calculate performance from a distance matrix and trial labels
%      P = compute_discrimination_performance(D, labels);
%
%   :param D: Pairwise distance matrix of dimensions `[T x T]` computed between all pooled trials.
%   :type D: matrix of doubles
%   :param labels: Class label vector of length `T` associating each trial with its stimulus index.
%   :type labels: vector of integers
%
%   :returns: **P** -- Discrimination performance score :math:`P`.
%   :type P: double
%
%   .. important::
%      This algorithms has a MEX version.
    
    T = length(labels);
    
    intra = zeros(1,2*T);
    inter = zeros(1,(T*(T-1)/2)-(2*T));
    idx_intra = 1;
    idx_inter = 1;

    for i = 1:T
    
        for j = i+1:T
    
            %% =====================================
            %% same stimulus
            %% =====================================
    
            if labels(i) == labels(j)

                intra(idx_intra) = D(i,j);
                idx_intra = idx_intra + 1;
    
            %% =====================================
            %% different stimuli
            %% =====================================
    
            else
    
                inter(idx_inter) = D(i,j);
                idx_inter = idx_inter + 1;
    
            end
        end
    end

    %% =====================================
    %% discrimination performance
    %% =====================================
    
    P = mean(inter) - mean(intra);

    
end


