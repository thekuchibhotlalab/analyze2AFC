function [tempWheel, badFlag] = fn_removeBadTrials(tempWheel)
    totalTrial = size(tempWheel,1);
    badFlag = tempWheel(:,end)==0;
    tempWheel(badFlag,:) = nan;
    disp(['Bad Trials = ' int2str(sum(badFlag)) ' out of ' int2str(totalTrial)])
end