%% ALGO Computation
% Author: Maxime BELTOISE
% Date: July 2026

function [spikes,responseMatrix] = generate_LL_dataset(params)

    N = params.N;
    S = params.S;
    R = params.R;
    
    %% ------------------------------------------------------------------------
    %% Response matrix
    %% ------------------------------------------------------------------------
    
    switch lower(params.mode)
    
        case 'structured'
    
            responseMatrix = params.responseMatrix;
    
            if size(responseMatrix,1)~=N || size(responseMatrix,2)~=S
                error('responseMatrix must be of size N x S.');
            end
    
        case 'random'
    
            if isfield(params,'connectionProbability')
                p = params.connectionProbability;
            else
                p = 0.4;
            end
    
            responseMatrix = rand(N,S) < p;
    
            % ensure each neuron responds to at least one stimulus
            for n = 1:N
                if ~any(responseMatrix(n,:))
                    responseMatrix(n,randi(S)) = 1;
                end
            end
    
        otherwise
            error('Unknown mode.')
    end
    
    %% ------------------------------------------------------------------------
    %% Dataset generation
    %% ------------------------------------------------------------------------
    
    spikes = cell(N,S,R);
    
    meanRate = params.meanRate;
    Tmax     = params.Tmax;
    jitter   = params.jitter;
    jitterIntensity = params.jitterIntensity;
    
    if ~isfield(params,'sameResponse')
        sameResponse = true(N,1);     
    else
        sameResponse = logical(params.sameResponse(:));
    
        if length(sameResponse) ~= N
            error('sameResponse must have length N.');
        end
    end
    
    for n = 1:N
    
        %==============================================================
        % Case 1: A single template shared among all stimuli
        %==============================================================
        if sameResponse(n)
    
            nSpikes_template = poissrnd(meanRate*Tmax);
            sharedTemplate = sort(rand(nSpikes_template,1)*Tmax);
    
        end
    
        for s = 1:S
    
            if responseMatrix(n,s)
    
                %======================================================
                % Case 2: A template specific to each stimulus
                %======================================================
                if ~sameResponse(n)
    
                    nSpikes_template = poissrnd(meanRate*Tmax);
                    sharedTemplate = sort(rand(nSpikes_template,1)*Tmax);
    
                end
    
                for r = 1:R
    
                    stim_shift = (s-1)*0.002*Tmax;
    
                    noise = (2*rand(size(sharedTemplate))-1) * jitter * jitterIntensity(n);
    
                    trial = sharedTemplate + stim_shift + noise;
    
                    trial = trial(trial>=0 & trial<=Tmax);
    
                    spikes{n,s,r} = sort(trial);
    
                end
    
            else
    
                %======================================================
                % Pas de réponse : bruit de fond
                %======================================================
                for r = 1:R
    
                    nSpikes = poissrnd(meanRate*Tmax);
    
                    spikes{n,s,r} = sort(rand(nSpikes,1)*Tmax);
    
                end
    
            end
    
        end
    
    end
    
    %% ------------------------------------------------------------------------
    %% Display
    %% ------------------------------------------------------------------------
    
    fprintf('\nLL response matrix\n');
    disp(responseMatrix)

end

