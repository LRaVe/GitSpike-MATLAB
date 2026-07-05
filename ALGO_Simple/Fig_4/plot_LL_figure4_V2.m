function plot_LL_figure4_V2(spikes,result,params,Distances,threshold)
    
    %======================================================================
    % FIGURE 4 - PAPER EXACT REPRODUCTION (FIXED)
    %======================================================================
    
    [N,S,R] = size(spikes);
    
    num_rows = 5;
    num_cols = 7;
    fs = 10;
    linewid = 1.5;
    
    figure('Color','w','Renderer','painters');
    set(gcf,'Position',[50 50 1314 1000]);
    
    ColormapLL = [
        0 0 0;
        0 0 1;
        1 0 0;
        0 0.5 0;
        1 1 0;
    ];
    
    %% ============================================================
    %% PANEL B (SPIKE TRAINS)
    %% ============================================================
    
    for nc = 1:N
        for sc = 1:S
    
            subplot(num_rows,num_cols,nc*num_cols + sc);
            hold on;
    
            for rc = 1:R
                st = spikes{nc,sc,rc};
                if isempty(st), continue; end
                plot(st, rc*ones(size(st)), 'k|', 'MarkerSize', 4);
            end
    
            set(gca,'xtick',[],'ytick',[]);
            xlim([0 params.Tmax]);
            ylim([0 R+1]);
    
            if nc == 1
                title(sprintf('S%d',sc),'FontWeight','bold','FontSize',fs);
            end
    
            if sc == 1
                ylabel(sprintf('N%d',nc));
            end
        end
    end
    
    %% ============================================================
    %% PANEL C / D / E (exact paper structure)
    %% ============================================================
    
    maxD = 0;
    maxP = 0;
    
    for n = 1:N
        maxD = max(maxD, max(result.DistanceMatrix{n}(:)));
        maxP = max(maxP, max(result.Mn{n}(:)));
    end
    
    
    
    %% ============================================================
    %% PANEL C
    %% ============================================================
    
    for nc = 1:N
    
        subplot(num_rows,num_cols,nc*num_cols+5)
    
        imagesc(result.DistanceMatrix{nc})
    
        axis square
        colormap(gca,jet)
        caxis([0 maxD])
    
        set(gca,'xtick',[],'ytick',[])
    
        title(sprintf('D_%d',nc))
        
    end
    
    
    %% ============================================================
    %% PANEL D
    %% ============================================================
    
    Colormap = [
        0 0 0
        0 0 1
        1 0 0
        0 0.5 0
        1 1 0];
    
    for nc = 1:N
    
        subplot(num_rows,num_cols,nc*num_cols+6)

        disp(['Neuron ',num2str(nc)])
        disp(result.Discrimination{nc})
        disp(unique(result.Discrimination{nc}*nc))
    
        imagesc(result.Discrimination{nc}*nc)
    
        axis square
        colormap(gca,Colormap)
        caxis([0 4])
    
        set(gca,'xtick',[],'ytick',[])
    
        title(sprintf('M_%d',nc))
    
        DrawSeparators(S,1,linewid)
    
    end
    
    
    %% ============================================================
    %% PANEL E
    %% ============================================================
    
    for nc = 1:N
    
        subplot(num_rows,num_cols,nc*num_cols+7)
    
        imagesc(result.Mn{nc})
    
        axis square
        colormap(gca,jet)
        caxis([0 maxP])
    
        set(gca,'xtick',[],'ytick',[])
    
        title(sprintf('P_%d',nc))
    
        DrawSeparators(S,1,linewid)
    
    end
        
    %% ============================================================
    %% PANEL F (population performance)
    %% ============================================================
    
    subplot(num_rows,num_cols,7);
    imagesc(result.populationPerformance);
    axis square;
    colormap(jet);
    set(gca,'xtick',[],'ytick',[]);
    title('F');
    
    DrawSeparators(S,1,linewid);
    
    %% ============================================================
    %% PANEL G (best neuron matrix)
    %% ============================================================
    
    Colormap = [...
        0 0 0
        0 0 1
        1 0 0
        0 0.5 0
        1 1 0];

    subplot(num_rows,num_cols,6);
    
    imagesc(result.bestNeuronMatrix);
    
    axis square;
    
    colormap(gca,Colormap);
    
    caxis([0 4]);
    set(gca,'xtick',[],'ytick',[]);
    title('G');
    
    DrawSeparators(S,1,linewid);
    
    %% ============================================================
    %% PRINT SUMMARY
    %% ============================================================
    
    fprintf('\n===== LL RESULTS =====\n');
    fprintf('Best population: ');
    fprintf('%d ', result.bestPopulation);
    fprintf('\nPLL = %.4f\n', result.bestP);
    fprintf('======================\n');

end



