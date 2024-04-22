function [bias,acc,dprime,crit,acc_L,acc_R,stimL,arL,arR] = fn_getAccBiasSmooth(stimulus,responseType,biasBin)
    if nargin<3; biasBin = 30; end
    correct = (responseType == 1);
    miss = (responseType == 0) | isnan(responseType);
    bias = nan(size(correct)); acc_L = nan(size(correct)); acc_R = nan(size(correct)); 
    acc = nan(size(correct)); dprime = nan(size(correct)); crit = nan(size(correct));
    stimL = nan(size(correct)); arL = nan(size(correct)); arR = nan(size(correct)); 
    nTrialThreshold = ceil(biasBin)/4;
    
    for j = biasBin:length(correct)
        tempCorrect = correct((j-biasBin+1):j);
        tempMiss = miss((j-biasBin+1):j);
        tempStim = stimulus((j-biasBin+1):j);
        
        tempStimL = tempStim==1 | tempStim==3;
        tempStimR = tempStim==-1 | tempStim==2 | tempStim==4;
        
        stimL(j - biasBin/2) = sum(tempStimL)/(sum(tempStimL)+sum(tempStimR));
        tempStim(tempStimL) = 1; tempStim(tempStimR) = 2;

        if sum(tempStimL) + sum(tempStimR) >= nTrialThreshold  
            [tempBias,tempAcc,tempDP,tempCrit,tempL,tempR,tempArL,tempArR] = fn_getAccBias(tempStim,tempCorrect,tempMiss);     
            bias(j - biasBin/2) = tempBias; % / (acc_s1 + acc_s2);
            acc_L(j - biasBin/2) = tempL; acc_R(j - biasBin/2) = tempR;
            acc(j - biasBin/2) = tempAcc; dprime(j - biasBin/2) = tempDP; crit(j - biasBin/2) = tempCrit;
            arL(j - biasBin/2) = tempArL; arR(j - biasBin/2) = tempArR;
        end
    end
end