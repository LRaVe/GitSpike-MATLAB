%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function plot_SP_figure(spikes,params,plotParams)
    
    [N,~,~] = size(spikes);
    
    codingNeurons = 1:params.c;
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
            %% separation coding / non-coding
            %% =====================================
    
            separationY = N - params.c + 0.5;
    
            yline(separationY,'--', 'Color',[0.4 0.4 0.4], 'LineWidth',1.2);
    
            %% =====================================
            %% POOLING
            %% =====================================
    
            if showPooling
    
                %% -----------------
                %% C
                %% -----------------
    
                pooledC = ...
                    pool_neurons(spikes,codingNeurons,s,r);
    
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
    
                yticks([yAll yNC yC neuronY(end) neuronY(1)]);
                
                yticklabels({'All', 'NC', 'C', sprintf('N%d',N), 'N1'});
                    
            else
    
                yticks([neuronY(1) neuronY(end)]);
    
                yticklabels({'1', num2str(N)});
    
            end
    
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

