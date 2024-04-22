function strCategory = fn_multistrcmpCategory(allStr, targetStrCell)
strCategory = nan(length(targetStrCell),length(allStr));
for i = 1:length(targetStrCell)
    targetStr = targetStrCell{i};
    [idxAll, ~, ~]= fn_multistrcmp(allStr, targetStr);
    strCategory(i,idxAll) = i;
    
end

end