function plotPerfBehavIndiv (mouseObj)
    %loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
    %modelComparison = 'SBARpRn';
    %nPrev = 1;
    lineWidthArg = {'LineWidth',2};
    trialLim = [1 size(mouseObj.behav,1)]; %trialLim = [1 2000];
    absFlag = false; 

    figure; margins = [0.01 0.05]; 
    subplot_tight(4,1,1,margins); hold on; plot(mouseObj.behav.acc,'Color', fn_wheelColorsPT('Reinf'),lineWidthArg{:});
    hold on; plot(mouseObj.probe.probeIdxBin,mouseObj.probe.probeAcc,'--o','Color',fn_wheelColorsPT('Probe'),lineWidthArg{:})
    plot([1 length(mouseObj.behav.acc)],[0.5 0.5],'Color',[0.8 0.8 0.8])
    xlim([1 length(mouseObj.behav.acc)]);xticks([]); ylabel('Accuracy')
    yticks([0.5 1]); ylim([0.4 1]); xlim(trialLim)

    subplot_tight(4,1,2,margins); hold on; 
    if ~absFlag
        plot(mouseObj.behav.bias,'Color',fn_wheelColorsPT('Reinf'),lineWidthArg{:});
        plot(mouseObj.probe.probeIdxBin,mouseObj.probe.probeBias,'--o','Color',fn_wheelColorsPT('Reinf'),lineWidthArg{:});
    else
        plot(abs(mouseObj.behav.bias),'Color',fn_wheelColorsPT('Reinf'),lineWidthArg{:});
        plot(mouseObj.probe.probeIdxBin,abs(mouseObj.probe.probeBias),'--o','Color',fn_wheelColorsPT('Probe'),lineWidthArg{:});
    end
    plot([1 length(mouseObj.behav.acc)],[0.4 0.4],'Color',[0.8 0.8 0.8])
    plot([1 length(mouseObj.behav.acc)],[-0.4 -0.4],'Color',[0.8 0.8 0.8])
    plot([1 length(mouseObj.behav.acc)],[-0 -0],'Color',[0.8 0.8 0.8])
    
    xlim([1 length(mouseObj.behav.acc)]);ylabel('Bias');
    xticks([]); xlim(trialLim);xticks(0:500:length(mouseObj.behav.acc))
    if ~absFlag; yticks([-1 0 1]); ylim([-1 1]);
    else; yticks([0 0.5 1]); ylim([0 1]);
    end
    
    subplot_tight(4,1,3,margins); hold on;
    plot(abs(mouseObj.behav.bias),'Color',fn_wheelColorsPT('Reinf'),lineWidthArg{:});
    %{
    filename=  [loadPath filesep mouseObj.ops.mouse 'psytrack_' modelComparison '_nPrev' int2str(nPrev) '.mat'];
    if exist(filename)
        load(filename);
        bW = wMode(4,:); tempbW = -((1+exp(-bW)).^(-1)*2-1);
        subplot_tight(5,1,3,margins); hold on;
        if ~absFlag
            plot(mouseObj.behav.bias,lineWidthArg{:}); 
            plot(tempbW,lineWidthArg{:});
        else 
            plot(abs(mouseObj.behav.bias),lineWidthArg{:}); 
            plot(abs(tempbW),lineWidthArg{:}); 
        end 
        plot([1 length(mouseObj.behav.acc)],[0.4 0.4],'Color',[0.8 0.8 0.8])
        plot([1 length(mouseObj.behav.acc)],[-0.4 -0.4],'Color',[0.8 0.8 0.8])
        plot([1 length(mouseObj.behav.acc)],[-0 -0],'Color',[0.8 0.8 0.8])
        xlim([1 length(mouseObj.behav.acc)]);
        if ~absFlag;yticks([-1 0 1]); ylim([-1 1]);
        else; yticks([0 0.5 1]); ylim([0 1]); 
        end
            
        xticks(0:500:length(mouseObj.behav.acc))
        ylabel('Model wight bias'); legend('Behavioral Bias','Model Bias')
        xlim(trialLim)
        
        
        subplot_tight(5,1,4,margins); hold on;
        if ~absFlag
            plot(tempbW,lineWidthArg{:});
        else  
            plot(abs(tempbW),lineWidthArg{:}); 
        end 
        plot([1 length(mouseObj.behav.acc)],[0.4 0.4],'Color',[0.8 0.8 0.8])
        plot([1 length(mouseObj.behav.acc)],[-0.4 -0.4],'Color',[0.8 0.8 0.8])
        plot([1 length(mouseObj.behav.acc)],[-0 -0],'Color',[0.8 0.8 0.8])
        xlim([1 length(mouseObj.behav.acc)]);
        if ~absFlag;yticks([-1 0 1]); ylim([-1 1]);
        else; yticks([0 0.5 1]); ylim([0 1]); 
        end
            
        xticks(0:500:length(mouseObj.behav.acc))
        ylabel('Model wight bias'); legend('Behavioral Bias','Model Bias')
        xlim(trialLim)
        
        
        subplot_tight(5,1,5,margins); hold on; 
        plot(mouseObj.behav.acc,'Color', fn_wheelColorsPT('Reinf'),lineWidthArg{:});
        hold on; plot(mouseObj.probe.probeIdxBin,mouseObj.probe.probeAcc,'--o','Color',fn_wheelColorsPT('Probe'),lineWidthArg{:})
        plot(mouseObj.behav.modelAcc,'Color', fn_wheelColorsPT('Reinf'),lineWidthArg{:});
        plot([1 length(mouseObj.behav.acc)],[0.5 0.5],'Color',[0.8 0.8 0.8])
        xlim([1 length(mouseObj.behav.acc)]);xticks([]); ylabel('Accuracy')
        yticks([0.5 1]); ylim([0.4 1]); xlim(trialLim)
    end
    %}
    

end
