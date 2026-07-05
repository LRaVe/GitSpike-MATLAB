%% ALGO Computation
% Author: Maxime BELTOISE
% Date: July 2026


function [Ppair,Mpair] = compute_LL_pair_performance(D,s1,s2,R)
    %==========================================================================
    % COMPUTE_LL_PAIR_PERFORMANCE
    %
    % Computes the discrimination performance of ONE neuron for ONE pair of
    % stimuli according to Delis et al.
    %
    % Inputs
    %   D  : trial-by-trial distance matrix for one neuron
    %   s1 : first stimulus
    %   s2 : second stimulus
    %   R  : repetitions per stimulus
    %
    % Outputs
    %   Ppair : discrimination performance (Eq.11)
    %   Mpair : logical discrimination result (Eq.10)
    %==========================================================================
    
    alpha = 0.001;
    
    %% ------------------------------------------------------------
    %% Trial indices
    %% ------------------------------------------------------------
    
    idx1 = (s1-1)*R + (1:R);
    idx2 = (s2-1)*R + (1:R);
    
    %% ------------------------------------------------------------
    %% Distance distributions
    %% ------------------------------------------------------------
    
    % Within stimulus 1
    D11 = D(idx1,idx1);
    D11 = D11(triu(true(R),1));
    
    % Within stimulus 2
    D22 = D(idx2,idx2);
    D22 = D22(triu(true(R),1));
    
    % Between stimuli
    D12 = D(idx1,idx2);
    D12 = D12(:);
    
    %% ------------------------------------------------------------
    %% Wilcoxon rank-sum tests (Eq.10)
    %% ------------------------------------------------------------
    
    try
        p1 = ranksum(D11,D12);
        p2 = ranksum(D22,D12);
        p3 = ranksum(D11,D22);
    catch
        % very small sample safety
        p1 = 1;
        p2 = 1;
        p3 = 1;
    end
    
    Mpair = (p1<alpha) || (p2<alpha) || (p3<alpha);
    
    %% ------------------------------------------------------------
    %% Discrimination performance (Eq.11)
    %% ------------------------------------------------------------
    
    if ~Mpair
    
        Ppair = 0;
        return
    
    end
    
    meanWithin = mean([D11(:); D22(:)]);
    meanBetween = mean(D12);
    
    den = meanBetween + meanWithin;
    
    if den<=eps
    
        Ppair = 0;
    
    else
    
        Ppair = (meanBetween-meanWithin)/den;
    
    end

end