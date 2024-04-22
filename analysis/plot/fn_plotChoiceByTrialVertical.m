function fn_plotChoiceByTrialVertical(plotMouse,tempTrialsIdx,varargin)

    p = inputParser;
    p.KeepUnmatched = true;
    % ------ Varargin for superclass
    p.addParameter('probeIdx', []);
    p.addParameter('csize', 25);
    p.parse(varargin{:});

    
    stim = plotMouse.behav.stimulus(tempTrialsIdx);% flagR = plotMouse.behav.stimulus(tempTrialsIdx) == 2;
    action = plotMouse.behav.action(tempTrialsIdx);
    trialNum = 1:length(tempTrialsIdx); trialNum = trialNum - 0.5;
    colors = zeros(length(trialNum),3);
    tempResp = plotMouse.behav.responseType(tempTrialsIdx);
    colors(tempResp == 1,:) = repmat(fn_wheelColorsPT('correct'),[sum(tempResp == 1) 1]);
    colors(tempResp == 2,:) = repmat(fn_wheelColorsPT('incorrect'),[sum(tempResp == 2) 1]);

    hold on;
    if ~isempty(p.Results.probeIdx)
        probeIdx = [p.Results.probeIdx(1)-1 p.Results.probeIdx];
        probeArea(1,:) = [0.5 2.5 2.5 0.5]; probeArea(2,:) = ...
            [p.Results.probeIdx(1)-1 p.Results.probeIdx(1)-1 p.Results.probeIdx(end) p.Results.probeIdx(end)] ;
        fill(probeArea(1,:),probeArea(2,:),fn_wheelColorsPT('Probe',0.3),'LineStyle','None');
    end

    scatter(action,trialNum,p.Results.csize,colors,'filled')

    scatter(stim(tempResp == 2),trialNum(tempResp == 2),p.Results.csize,fn_wheelColorsPT('correct'))
    %scatter(trialNum(flagR),action,10,matlabColors(2),'filled')
    %legend('Stim 1','Stim 2','AutoUpdate','off')
    xticks([1 2]);xticklabels({'L','R'});
    ylabel('Trials')

end