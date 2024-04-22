function fn_plotChoiceByTrial(plotMouse,tempTrialsIdx)
    
    stim = plotMouse.behav.stimulus(tempTrialsIdx);% flagR = plotMouse.behav.stimulus(tempTrialsIdx) == 2;
    action = plotMouse.behav.action(tempTrialsIdx);
    trialNum = 1:length(tempTrialsIdx);
    colors = zeros(length(trialNum),3);
    tempResp = plotMouse.behav.responseType(tempTrialsIdx);
    colors(tempResp == 1,:) = repmat(fn_wheelColorsPT('correct'),[sum(tempResp == 1) 1]);
    colors(tempResp == 2,:) = repmat(fn_wheelColorsPT('incorrect'),[sum(tempResp == 2) 1]);
    hold on;
    scatter(trialNum,action,25,colors,'filled')

    scatter(trialNum(tempResp == 2),stim(tempResp == 2),25,fn_wheelColorsPT('correct'))
    %scatter(trialNum(flagR),action,10,matlabColors(2),'filled')
    %legend('Stim 1','Stim 2','AutoUpdate','off')
    yticks([1 2]);yticklabels({'L','R'});
    xlabel('Trials')

end