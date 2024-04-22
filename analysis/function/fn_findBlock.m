function [startIdx,endIdx,len] = fn_findBlock(stateFlag,element,blockIntervalThre)
    idx = find(stateFlag == element);incre = diff(idx);
    startFlag = find(incre>blockIntervalThre); 
    if ~isempty(idx)
        startIdx = idx([1; startFlag+1]) ; endIdx = idx([startFlag; length(idx)]);
        len = endIdx - startIdx + 1;
    else 
        startIdx =[]; endIdx = []; len = []; 
    end
end