function [wheelDownSample_all, RT_all,selFlag] = loadWheelClusterData(mouseMega, idx)
    if isempty(idx); idx = 1:mouseMega.nMouse; end
    nTrialsum = [0 cumsum(mouseMega.nTrials)];
    [wheelDownSample,badFlag] = getDownSample(mouseMega,idx);
    RT_all = fn_cell2mat(mouseMega.getProp('behav','field','reactionTime','matFlag',false,'idx',idx),1);
    wheelDownSample_all = fn_cell2mat(wheelDownSample,1);

    selFlag = [];
    for i = 1:length(idx); selFlag = [selFlag (nTrialsum(idx(i))+1 : nTrialsum(idx(i)+1))];end
end