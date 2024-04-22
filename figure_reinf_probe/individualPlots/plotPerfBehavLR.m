function plotPerfBehavLR (mouseObj)
    %loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
    %modelComparison = 'SBARpRn';
    %nPrev = 1;
    lineWidthArg = {'LineWidth',2};
    trialLim = [1 size(mouseObj.behav,1)]; %trialLim = [1 2000];
    absFlag = false; 

    figure; 
    subplot(2,1,1); hold on; 
    plot(mouseObj.behav.acc1,'Color', matlabColors(1),lineWidthArg{:});
    plot(mouseObj.behav.acc2,'Color', matlabColors(2),lineWidthArg{:});
    hold on; 
    plot([1 length(mouseObj.behav.acc)],[0.5 0.5],'Color',[0.8 0.8 0.8])
    plot([1 length(mouseObj.behav.acc)],[0.7 0.7],'Color',[0.8 0.8 0.8] ,'LineStyle','--')
    xlim([1 length(mouseObj.behav.acc)]);xticks([]); ylabel('Accuracy')
    yticks([0.5 1]); ylim([0 1]); xlim(trialLim)

    subplot(2,1,2); hold on; 
    if ~absFlag
        plot(mouseObj.behav.bias,'Color',fn_wheelColorsPT('Reinf'),lineWidthArg{:});
    else
        plot(abs(mouseObj.behav.bias),'Color',fn_wheelColorsPT('Reinf'),lineWidthArg{:});
    end
    plot([1 length(mouseObj.behav.acc)],[0 0],'Color',[0.8 0.8 0.8])
    
    xlim([1 length(mouseObj.behav.acc)]);ylabel('Bias');
    xticks([]); xlim(trialLim);xticks(0:500:length(mouseObj.behav.acc))
    if ~absFlag; yticks([-1 0 1]); ylim([-1 1]);
    else; yticks([0 0.5 1]); ylim([0 1]);
    end
    
    
    

end
