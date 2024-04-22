function [bias] = fn_generateRandomChoice(acc,binSize)
    w = -log(1/acc-1);
    nSample = 5000*binSize;
    nCounterBalance = 4;

    randStim = zeros(nCounterBalance,nSample/nCounterBalance); 
    for i = 1:nSample/nCounterBalance
        randStim(:,i) = randperm(nCounterBalance)/nCounterBalance;
    end
    randStim = randStim(:)';
    randStim(randStim<=0.5) = -1; randStim(randStim>0.5) = 1;
    logW = ones(1,nSample)*w;
    choiceProb = 1./(1+exp(-logW.*randStim));
    randStim1 = randStim == -1; randStim2 = randStim == 1;
    
    randChoice = rand(1,nSample);c2 = (randChoice<=choiceProb); c1 = (randChoice>choiceProb);
    randChoice(c1) = -1; randChoice(c2) = 1;
    acc1 = (randChoice==randStim).*randStim1; acc1(randStim2) = nan;
    acc2 =  (randChoice==randStim).*randStim2; acc2(randStim1) = nan;
    acc1 = reshape(acc1,binSize,[]);acc1 = nanmean(acc1,1);
    acc2 = reshape(acc2,binSize,[]);acc2 = nanmean(acc2,1);
    bias = acc1 - acc2; 
    



end