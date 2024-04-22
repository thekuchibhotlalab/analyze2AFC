function [bias,acc] = fn_generateChoiceFromAgent(trialLength,stimW, biasW)
    nRep = 200; 
    randStim = rand(nRep,trialLength); 
    randStim(randStim<0.5) = -1; randStim(randStim>=0.5) = 1;
    randStim1 = randStim == -1; randStim2 = randStim == 1;

    if nargin == 1
        stimW = zeros(size(randStim));
        biasW = zeros(size(randStim));
    else
        if isempty(stimW); stimW = zeros(1,trialLength); end
        if isempty(biasW); biasW = zeros(1,trialLength); end
        if length(stimW) == 1; stimW = repmat(stimW,[1 trialLength]); end
        if length(biasW) == 1; biasW = repmat(biasW,[1 trialLength]); end
        if size(stimW,2) == 1 && length(stimW)~=1; stimW = stimW'; end
        if size(biasW,2) == 1 && length(biasW)~=1; biasW = biasW'; end
        stimW = repmat(stimW,[nRep 1]); 
        biasW = repmat(biasW,[nRep 1]); 
    end

    choiceProb = 1./(1+exp(-stimW.*randStim+biasW));
    randChoice = rand(size(randStim));c2 = (randChoice<=choiceProb); c1 = (randChoice>choiceProb);
    randChoice(c1) = -1; randChoice(c2) = 1;
    acc1 = (randChoice==randStim).*randStim1; acc1(randStim2) = nan;
    acc2 =  (randChoice==randStim).*randStim2; acc2(randStim1) = nan;
    acc1 = smoothdata(acc1,2,'movmean',30); acc2 = smoothdata(acc2,2,'movmean',30);
    acc = acc1/2+acc2/2; bias = smoothdata(acc2-acc1,2,'gaussian',30);

end