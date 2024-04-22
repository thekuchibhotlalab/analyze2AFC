function [bias,acc_L,acc_R] = fn_getBiasModel(stimulus,correct,biasBin)
    if nargin<3; biasBin = 30; end
    bias = nan(size(correct)); acc_L = nan(size(correct)); acc_R = nan(size(correct));
    for j = biasBin:length(correct)
        tempCorrect = correct((j-biasBin+1):j);
        tempStim = stimulus((j-biasBin+1):j);
        acc_s1 = nanmean(tempCorrect (tempStim==-1) );
        acc_s2 = nanmean(tempCorrect (tempStim==1) );
        bias(j - biasBin/2) = (acc_s1 - acc_s2); % / (acc_s1 + acc_s2);
        acc_L(j - biasBin/2) = acc_s1; acc_R(j - biasBin/2) = acc_s2;
    end
end