%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function plot_SP_figure(spikes,params,plotParams)
% PLOT_SP_FIGURE Plots raster graphs for individual and pooled spike trains.
%
%   Generates a multi-panel figure displaying raster plots of individual neuron 
%   spike trains alongside pooled spike trains (Coding, Non-Coding, and All neurons). 
%   Spikes are color-coded based on neuron type (Coding in red, Individual/Noise in magenta, 
%   and Non-Coding in blue). Dashed separation lines visually divide the different functional groups.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Plot raster with pooled spike trains enabled
%      plotParams.stimuli = [1, 2];
%      plotParams.repetitions = [1];
%      plotParams.showPooling = true;
%      plot_SP_figure(spikes, params, plotParams);
%
%      % Basic plot showing only individual raster lines
%      plotParams.showPooling = false;
%      plot_SP_figure(spikes, params, plotParams);
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param params: Structure containing simulation parameters with fields:
%
%                  * **c** (*integer*): Number of coding neurons.
%                  * **nIndi** (*integer*): Number of individual/noise-contributing neurons.
%                  * **Tmax** (*double*): Upper temporal boundary of the window.
%   :type params: struct
%   :param plotParams: Structure containing figure configuration parameters with fields:
%
%                      * **stimuli** (*vector*): Indices of stimuli to plot (e.g., `[1, 2]`).
%                      * **repetitions** (*vector*): Indices of repetitions to plot (e.g., `[1]`).
%                      * **showPooling** (*logical*): Flag to toggle displaying pooled spike trains (`C`, `NC`, `All`).
%   :type plotParams: struct
%
%   .. note::
%      When `showPooling` is set to `true`, the figure adds three extra rows at the bottom 
%      representing pooled spike trains: **C** (Coding neurons), **NC** (Non-Coding neurons), 
%      and **All** (Entire population).

    [N,~,~] = size(spikes);
    
    codingNeurons = 1:params.c;
    indiNeurons = params.c+1:params.c+params.nIndi;
    nonCoding = params.c+1:N;
    
    stimuli = plotParams.stimuli;
    repetitions = plotParams.repetitions;
    
    showPooling = plotParams.showPooling;
    
    nPanels = length(stimuli)*length(repetitions);
    
    figure;
    set(gcf,'Color','w');
    
    plotIndex = 1;
    
    for s = stimuli
    
        for r = repetitions
    
            subplot(nPanels,1,plotIndex);
            hold on;
    
            %% =====================================
            %% AXES
            %% =====================================
    
            xlim([0 params.Tmax]);
    
            %% =====================================
            %% POSITIONS Y
            %% =====================================
    
            % individual neurons
            neuronY = N:-1:1;
    
            % pooled trains
            yC   = 0;
            yNC  = -1;
            yAll = -2;
    
            %% =====================================
            %% INDIVIDUAL SPIKE TRAINS
            %% =====================================
    
            for n = 1:N
    
                t = spikes{n,s,r};
    
                %% color
    
                if ismember(n,codingNeurons)
                    col = [1 0 0];
                elseif ismember(n,indiNeurons)
                    col = [1 0 1];
                else
                    col = [0 0 1];
                end
    
                %% position y
    
                y = neuronY(n);
    
                %% spikes
    
                for k = 1:length(t)
    
                    line([t(k) t(k)], [y-0.4 y+0.4], 'Color',col, 'LineWidth',1.5);
    
                end
            end
    
            %% =====================================
            %% DASHED LINE
            %% separation coding / indi / non-coding
            %% =====================================
    
            separationY = N - params.c + 0.5;
    
            yline(separationY,'--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

            separationY2 = N - (params.c + params.nIndi) + 0.5;
    
            yline(separationY2,'--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);

    
            %% =====================================
            %% POOLING
            %% =====================================
    
            if showPooling
    
                %% -----------------
                %% C
                %% -----------------
    
                pooledC = pool_neurons(spikes,codingNeurons,s,r);
    
                for k = 1:length(pooledC)
    
                    line([pooledC(k) pooledC(k)], [yC-0.4 yC+0.4], 'Color',[1 0 0], 'LineWidth',1.5);
    
                end
    
                %% separation line
    
                yline(yC+0.5,'k-','LineWidth',1);
    
                %% -----------------
                %% NC
                %% -----------------
    
                pooledNC = pool_neurons(spikes,nonCoding,s,r);
    
                for k = 1:length(pooledNC)
    
                    line([pooledNC(k) pooledNC(k)], [yNC-0.4 yNC+0.4], 'Color',[0 0 1], 'LineWidth',1.5);
    
                end
    
                %% separation line
    
                yline(yNC+0.5,'--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    
                %% -----------------
                %% ALL
                %% -----------------
    
                pooledAll = pool_neurons(spikes,1:N,s,r);
    
                for k = 1:length(pooledAll)
    
                    line([pooledAll(k) pooledAll(k)], [yAll-0.4 yAll+0.4], 'Color', 'k', 'LineWidth',1.5);
    
                end

                %% separation line
    
                yline(yAll+0.5,'--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    
                %% =====================================
                %% TICKS
                %% =====================================
    
                yticks([yAll yNC yC neuronY(end) neuronY(params.c+1) neuronY(1)]);
                
                yticklabels({'All', 'NC', 'C', sprintf('N%d',N),sprintf('N%d',params.c+1), 'N1'});
                    
            else
    
                yticks([neuronY(1) neuronY(end)]);
    
                yticklabels({'1', num2str(N)});
    
            end

            %% =====================================
            %% remove tick marks
            %% =====================================
        
            set(gca,'TickLength',[0 0]);
    
            %% =====================================
            %% LIMITS
            %% =====================================
    
            ylim([yAll-1 N+1]);
    
            %% =====================================
            %% STYLE
            %% =====================================
    
            xlabel('Time');
    
            ylabel('Spike trains');
    
            title(sprintf('S%d-R%d',s,r), 'FontWeight', 'bold');
    
            set(gca,'FontSize',11);
    
            box on;
    
            plotIndex = plotIndex + 1;
    
        end
    end
end

