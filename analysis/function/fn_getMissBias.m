function bias = fn_getMissBias(stimulus,responeType)
biasBin = 50;
miss = (responeType == 0);
bias = nan(size(miss));
for j = biasBin:length(miss)
    tempMiss = miss((j-biasBin+1):j);
    tempStim = stimulus((j-biasBin+1):j);
    acc_s1 = sum(tempMiss & tempStim==1) / sum(tempStim==1);
    acc_s2 = sum(tempMiss & tempStim==2) / sum(tempStim==2);
    if sum(tempMiss) < biasBin*0.2
        bias(j - biasBin/2) = nan;
    else
        bias(j - biasBin/2) = (acc_s1 - acc_s2) / (acc_s1 + acc_s2);
    end
end

end




