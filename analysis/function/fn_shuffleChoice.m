function fn_shuffleChoice(mouseObj)
    actionBiasThre = 0.3; 
    nShuffle  = 1000; biasBin = 30;
    %shuffIdx = cellfun(@(x)randperm(x),num2cell(ones(nShuffle,1)*length(mouseObj.behav.action)),'UniformOutput',false);
    permChoice = cellfun(@(x)(round(rand(x,1)+1)),num2cell(ones(nShuffle,1)*length(mouseObj.behav.action)),'UniformOutput',false);
    
    
    permResponseType = cellfun(@(x)double((x==mouseObj.behav.stimulus)),permChoice,'UniformOutput',false);
    for i = 1:length(permResponseType); permResponseType{i}(permResponseType{i}==0) = 2; end 
    
    permAcc = cellfun(@(x)smoothdata(x==1,1,'movmean',biasBin),permResponseType,'UniformOutput',false);
    
    [permBias,~,~] = cellfun(@(x)fn_getBias(mouseObj.behav.stimulus,x,biasBin),...
        permResponseType,'UniformOutput',false);

    nBlockL = []; nBlockR = []; lenBlockL = []; lenBlockR = [];
    for i =1:length(permBias)
        [~, blockL, blockR, ~]=fn_detectBlock(permBias{i},0, actionBiasThre,'blockLenThre',10,'blockIntervalThre',10);
        nBlockL(i) = length(blockL.len);nBlockR(i) = length(blockR.len);
        lenBlockL = cat(1,lenBlockL,blockL.len); lenBlockR = cat(1,lenBlockR,blockR.len);
    end
    meanNBlock = mean([(nBlockL) (nBlockR)]);
    prcNBlock = prctile([nBlockL nBlockR],95);
    
    meanBlockLen = mean([(lenBlockL);(lenBlockR)]);
    prcBlockLen = prctile([lenBlockL;lenBlockR],95);
    
    tempBias = fn_cell2mat(permBias,2);
    meanBias = mean(tempBias,2);
    stdBias = std(tempBias,0,2);
    
    figure; 
    subplot(2,2,1); plot(permBias{1}); hold on; title('Example of bias + threhold ')
    plot(meanBias+actionBiasThre,'Color',[0.8 0.8 0.8]); plot(meanBias-actionBiasThre,'Color',[0.8 0.8 0.8])
    ylim([-0.6 0.6])
    subplot(2,2,2); plot(stdBias); title('Std of Bias');ylim([0 0.35])
    subplot(2,2,3); histogram([nBlockL nBlockR],0:2:30);
    xlabel('nBlock'); ylabel('frequency of 1000 rep')
    title(['mean= ' int2str(meanNBlock) ', 95 prctile=' int2str(prcNBlock)])
    subplot(2,2,4); histogram([lenBlockL;lenBlockR],0:2:40); 
    xlabel('block length'); ylabel('frequency of blocks, 1000 rep')
    title(['mean= ' int2str(meanBlockLen) ', 95 prctile=' int2str(prcBlockLen)])
end