function plotBlockTransitionFreq(mouseObj)
    if ~isobject(mouseObj); return; end
    if ~strcmp(class(mouseObj),'wheel2AFCmega'); return; end 
    
    %------------------MAKE PLPOTS BY ALIGNING BY TRIAL-------------
    global nAni; 
    nAni = mouseObj.nMouse;
    bins = 0:500:3000;
    binAxis = (bins(1:end-1) + bins(2:end))/2;
    
    %bias = abs(getProp(mouseObj,'behav','field','bias','matFlag',true));
    bias = abs(getProp(mouseObj,'behav','field','bias','matFlag',true));
    stateFlag = fn_cell2matFillNan(cellfun(@(x)(x.biasBlock.stateFlag'),mouseObj.mouseCell,'UniformOutput',false));
    
    highBiasCount = getBiasFreq(bias,stateFlag~=0,[0.4 1],bins);
    nanCount = sum(~isnan(highBiasCount),2);
    semCount = nanstd(highBiasCount,0,2)./sqrt(nanCount);
    figure; subplot(1,4,1);%scatter(1:size(highBiasCount,1),nanmean(highBiasCount,2))
    errorbar(binAxis,nanmean(highBiasCount,2),semCount,'LineWidth',2);
    title('Trials in high-biased block');ylim([0 1]);yticks([0 0.5 1]);xticks([0 1500 3000])
    
    lowBiasCount = getBiasFreq(bias,stateFlag~=0,[0.2 0.4],bins);
    nanCount = sum(~isnan(lowBiasCount),2);
    semCount = nanstd(lowBiasCount,0,2)./sqrt(nanCount);
    subplot(1,4,2); %scatter(1:size(lowBiasCount,1),nanmean(lowBiasCount,2))
    errorbar(binAxis,nanmean(lowBiasCount,2),semCount,'LineWidth',2);
    title('Trials in low-biased block');ylim([0 1]);yticks([0 0.5 1]);xticks([0 1500 3000])
    
    unBiasCount = getBiasFreq(bias,[],[0 0.2],bins);
    nanCount = sum(~isnan(unBiasCount),2);
    semCount = nanstd(unBiasCount,0,2)./sqrt(nanCount);
    subplot(1,4,3); %scatter(1:size(unBiasCount,1),nanmean(unBiasCount,2))
    errorbar(binAxis,nanmean(unBiasCount,2),semCount,'LineWidth',2);
    title('Trials in unbiased block');ylim([0 1]); yticks([0 0.5 1]);xticks([0 1500 3000])
    
    trans = (getProp(mouseObj,'biasBlock','field','trans'));
    transCount = fn_cell2mat(cellfun(@(x)getTransFreq(x,bins),trans,'UniformOutput',false),1)';
    transCount(isnan(highBiasCount)) = nan;
    nanCount = sum(~isnan(highBiasCount),2);
    semCount = nanstd(transCount,0,2)./sqrt(nanCount);
    subplot(1,4,4); %scatter(1:size(unBiasCount,1),nanmean(unBiasCount,2))
    errorbar(binAxis,nanmean(transCount,2),semCount,'LineWidth',2);
    title('nTransitions');ylim([0 8]);yticks([0 4 8]);xticks([0 1500 3000])
    
    %biasStart = nanmean(highBiasCount(1:2,:),1);
    %biasEnd = nanmean(highBiasCount((end-2):(end-1),:),1);
    figure; 
    subplot(1,4,1);plotStartEndbar(highBiasCount); ylim([0 1]); yticks([0 0.5 1])
    subplot(1,4,2);plotStartEndbar(lowBiasCount);ylim([0 1]); yticks([0 0.5 1])
    subplot(1,4,3);plotStartEndbar(unBiasCount);ylim([0 1]); yticks([0 0.5 1])
    subplot(1,4,4);plotStartEndbar(transCount);ylim([0 8]); yticks([0 4 8])
    
    blockL = getProp(mouseObj,'biasBlock','field','blockL');
    lenL = cellfun(@(x)(x.len),blockL,'UniformOutput',false);
    trialL = cellfun(@(x)(x.start/2+x.end/2),blockL,'UniformOutput',false);
    biasL = cellfun(@(x)(x.bias),blockL,'UniformOutput',false);
    
    blockR = getProp(mouseObj,'biasBlock','field','blockR');
    lenR = cellfun(@(x)(x.len),blockR,'UniformOutput',false);
    trialR = cellfun(@(x)(x.start/2+x.end/2),blockR,'UniformOutput',false);
    biasR = cellfun(@(x)(x.bias),blockR,'UniformOutput',false);

    blockU = getProp(mouseObj,'biasBlock','field','blockU');
    lenU = cellfun(@(x)(x.len),blockU,'UniformOutput',false);   
    trialU = cellfun(@(x)(x.start/2+x.end/2),blockU,'UniformOutput',false);

    lenByTrial = []; biasByTrial = [];
    for i = 1:nAni
        for j = 1:length(bins)-1
            temp = [lenL{i}(trialL{i} <= bins(j+1) & trialL{i} >= bins(j)+1) lenR{i}(trialR{i} <= bins(j+1) & trialR{i} >= bins(j)+1 )];
            lenByTrial(j,i) = nanmean(temp); 
            temp = [biasL{i}(trialL{i} <= bins(j+1) & trialL{i} >= bins(j)+1) biasR{i}(trialR{i} <= bins(j+1) & trialR{i} >= bins(j)+1 )];
            biasByTrial(j,i) = nanmean(temp); 
        end
    end
    nanCount = sum(~isnan(lenByTrial),2);
    semLen = nanstd(lenByTrial,0,2)./sqrt(nanCount);

    figure; subplot(1,2,1);
    errorbar(binAxis,nanmean(lenByTrial,2),semLen,'LineWidth',2);
    title('Average Block Length');xticks([0 1500 3000])
    subplot(1,2,2);
    plotStartEndbar(lenByTrial);

    nanCount = sum(~isnan(biasByTrial),2);
    semBias = nanstd(biasByTrial,0,2)./sqrt(nanCount);

    figure; subplot(1,2,1);
    errorbar(binAxis,nanmean(biasByTrial,2),semBias,'LineWidth',2);
    title('Average Block Bias Level');xticks([0 1500 3000])
    subplot(1,2,2);
    plotStartEndbar(biasByTrial);

    %------------------MAKE PLPOTS BY ALIGNING BY PERFORMANCE-------------
    % PLOT THE TIME SPENT IN BIASED BLOCKS WITH ALIGNED LEARNING CURVE
    probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
    outCell = mouseObj.objFun('binProbeByTrialFromLearningOnset',{[reinfThre probeThre],probeTrialBin});
    % reinf data
    reinfData = mouseObj.loadReinf; reinfData.stateFlag = stateFlag;
    reinfAlignLoc = cell2mat(outCell{3});
    [reinfDataAlign,reinfAlignPoint] = attachNan(reinfData, reinfAlignLoc);
    % plot axis
    reinfAxis = (1:size(reinfDataAlign.acc,1)) - reinfAlignPoint;
    reinfStart = reinfAlignPoint - 1000; reinfEnd = reinfAlignPoint + 2000; 

    highBiasCount = getBiasFreq(reinfDataAlign.bias,reinfDataAlign.stateFlag~=0,[0.4 1],bins+reinfStart);
    nanCount = sum(~isnan(highBiasCount),2);
    semCount = nanstd(highBiasCount,0,2)./sqrt(nanCount);
    figure; subplot(1,3,1); hold on;%scatter(1:size(highBiasCount,1),nanmean(highBiasCount,2))
    errorbar(binAxis-1000,nanmean(highBiasCount,2),semCount,'LineWidth',2);
    unBiasCount = getBiasFreq(reinfDataAlign.bias,[],[0 0.2],bins+reinfStart);
    nanCount = sum(~isnan(unBiasCount),2);
    semCount = nanstd(unBiasCount,0,2)./sqrt(nanCount);
    errorbar(binAxis-1000,nanmean(unBiasCount,2),semCount,'LineWidth',2);
    title('Trials in high/unbiased block');ylim([0 1]); yticks([0 0.5 1]);xticks(-1000:1000:2000)

    lowBiasCount = getBiasFreq(reinfDataAlign.bias,reinfDataAlign.stateFlag~=0,[0.2 0.4],bins+reinfStart);
    nanCount = sum(~isnan(lowBiasCount),2);
    semCount = nanstd(lowBiasCount,0,2)./sqrt(nanCount);
    subplot(1,3,2); %scatter(1:size(lowBiasCount,1),nanmean(lowBiasCount,2))
    errorbar(binAxis-1000,nanmean(lowBiasCount,2),semCount,'LineWidth',2);
    title('Trials in low-biased block');ylim([0 1]);yticks([0 0.5 1]);xticks(-1000:1000:2000)

    

    transAlign = cellfun(@(x,y)(x-y),trans,num2cell(reinfAlignLoc),'UniformOutput',false);
    transCount = fn_cell2mat(cellfun(@(x)getTransFreq(x,bins-1000),transAlign,'UniformOutput',false),1)';
    transCount(isnan(highBiasCount)) = nan;
    semCount = nanstd(transCount,0,2)./sqrt(nanCount);
    subplot(1,3,3); %scatter(1:size(unBiasCount,1),nanmean(unBiasCount,2))
    errorbar(binAxis-1000,nanmean(transCount,2),semCount,'LineWidth',2);
    title('nTransitions');ylim([0 8]);yticks([0 4 8]);xticks(-1000:1000:2000)

    figure; 
    subplot(1,4,1);plotStartEndbar(highBiasCount); ylim([0 1]); yticks([0 0.5 1])
    subplot(1,4,2);plotStartEndbar(lowBiasCount);ylim([0 1]); yticks([0 0.5 1])
    subplot(1,4,3);plotStartEndbar(unBiasCount);ylim([0 1]); yticks([0 0.5 1])
    subplot(1,4,4);plotStartEndbar(transCount);ylim([0 8]); yticks([0 4 8])

    lenByTrial = []; biasByTrial = []; newBins = bins - 1000; 
    for i = 1:nAni
        newTL = trialL{i}-reinfAlignLoc(i); newTR = trialR{i}-reinfAlignLoc(i);
        for j = 1:length(newBins)-1
            temp = [lenL{i}(newTL <= newBins(j+1) & newTL >= newBins(j)+1) lenR{i}(newTR <= newBins(j+1) & newTR >= newBins(j)+1 )];
            lenByTrial(j,i) = nanmean(temp); 
            temp = [biasL{i}(newTL <= newBins(j+1) & newTL >= newBins(j)+1) biasR{i}(newTR <= newBins(j+1) & newTR >= newBins(j)+1 )];
            biasByTrial(j,i) = nanmean(temp); 
        end
    end
    nanCount = sum(~isnan(lenByTrial),2);
    semLen = nanstd(lenByTrial,0,2)./sqrt(nanCount);

    figure; subplot(1,2,1);
    errorbar(binAxis-1000,nanmean(lenByTrial,2),semLen,'LineWidth',2);
    title('Average Block Length');xticks([0 1500 3000])
    subplot(1,2,2);
    plotStartEndbar(lenByTrial);

    nanCount = sum(~isnan(biasByTrial),2);
    semBias = nanstd(biasByTrial,0,2)./sqrt(nanCount);

    figure; subplot(1,2,1);
    errorbar(binAxis-1000,nanmean(biasByTrial,2),semBias,'LineWidth',2);
    title('Average Block Bias Level');xticks([0 1500 3000]); ylim([0.2 0.5]); yticks(0.2:0.1:0.5);
    subplot(1,2,2);
    plotStartEndbar(biasByTrial);
    
   
    
    %{
    reinfThre = [0.6 0.75 0.85]; accBin = 300;
    
    highBiasCount = []; lowBiasCount = []; unBiasCount = [];
    for i = 1:nAni
        [~,acc] = fn_getAccBiasSmooth(mouseObj.mouseCell{i}.behav.stimulus,mouseObj.mouseCell{i}.behav.responseType,accBin);
        reinfThreTrial = cellfun(@(x)find(acc>x,1), num2cell(reinfThre),'UniformOutput',false);
        emptyFlag = cellfun(@isempty,reinfThreTrial);
        reinfThreTrial(emptyFlag) = {length(acc)}; reinfThreTrial = fn_cell2mat(reinfThreTrial,1);
        reinfThreTrial = [0;reinfThreTrial;length(acc)];

        highBiasCount(:,i) = getBiasFreq(bias(:,i),stateFlag(:,i)~=0,[0.4 1],reinfThreTrial);
        lowBiasCount(:,i) = getBiasFreq(bias(:,i),stateFlag(:,i)~=0,[0.2 0.4],reinfThreTrial);
        unBiasCount(:,i) = getBiasFreq(bias(:,i),[],[0 0.2],reinfThreTrial);
    end 
    figure; subplot(1,3,1);plot(nanmean(highBiasCount,2));title('Trials in high-biased block');
    xticks(1:4); xticklabels({'acc<0.6', '0.6-0.75','0.7-0.85','>0.85'})
    subplot(1,3,2);plot(nanmean(lowBiasCount,2));title('Trials in low-biased block')
    xticks(1:4); xticklabels({'acc<0.6', '0.6-0.75','0.7-0.85','>0.85'})
    subplot(1,3,3);plot(nanmean(unBiasCount,2));title('Trials in unbiased block')
    xticks(1:4); xticklabels({'acc<0.6', '0.6-0.75','0.7-0.85','>0.85'})
    %}
end

function tempBiasCount = getBiasFreq(bias,stateFlag,biasThre,bins)
    if ~isempty(stateFlag)
        biasFlag = bias > biasThre(1) & bias < biasThre(2) & stateFlag;
    else
        biasFlag = bias > biasThre(1) & bias < biasThre(2);
    end
    trialFlag = ~isnan(bias);
    tempBiasCount =[];
    for i = 1:length(bins)-1
        tempFlag = (bins(i)+1):bins(i+1);% if there is no trials in this bin 
        if ~isempty(tempFlag)
            tempB = biasFlag(tempFlag,:); tempTotal = trialFlag(tempFlag,:);
            tempTotal = sum(tempTotal==1,1); tempB = sum(tempB==1,1);
            tempB(tempTotal<200) = nan;
            tempBiasCount(i,:) = tempB./tempTotal;
        else; tempBiasCount(i,:) = nan;
        end
    end 
    
end

function plotStartEndbar(countStartEnd)
    global nAni;
    biasStart = zeros(1,nAni); biasEnd = zeros(1,nAni); nBlock = 2; 
    for i = 1:nAni
        % keepFlag = find(~isnan(countStartEnd(:,i)));
        % startIdx = keepFlag(1:nBlock);  endIdx = keepFlag(end-nBlock+1:end);
        startIdx = 1:nBlock; 
        endIdx = size(countStartEnd,1)-nBlock+1:size(countStartEnd,1);
        biasStart(i) = nanmean(countStartEnd(startIdx,i));
        biasEnd(i) = nanmean(countStartEnd(endIdx,i));
    end
    nanFlag = isnan(biasEnd); 
    biasStart(nanFlag) = []; biasEnd(nanFlag) = [];
    [h,p] = ttest(biasStart,biasEnd);
    bar([nanmean(biasStart) nanmean(biasEnd)] ,'EdgeColor',[0 0 0],'FaceColor','None'); hold on;
    for i = 1:length(biasStart)
        f = plot([1 2],[biasStart(i) biasEnd(i) ],'Color',[0.6 0.6 0.6],'Marker','.','MarkerSize',15,...
            'MarkerFaceColor',[0.6 0.6 0.6],'LineWidth',0.5);
    end
    legend(f,['pBef = ' num2str(p,'%.2e') ],'Location','Best')
    xlim([0 3]);
    xticks([1 2 3]); xticklabels({'start','end'}); 
    
    
end

function transCount = getTransFreq(trans,bins)
    transIdx = trans(:,2); transCount = histcounts(transIdx,bins);
    %transCount = transCount ./ sum(transCount(1:2));
end

function [probeData,maxAlignPoint] = attachNan(probeData, alignPoint)
    attachDim = 1; 
    maxAlignPoint = max(alignPoint);
    tempFieldNames = fieldnames(probeData);
    for i = 1:length(tempFieldNames)
        tempField = probeData.(tempFieldNames{i});
        if isnumeric(tempField)
            tempMatSize = size(tempField);
            tempMatSize(attachDim) = tempMatSize(attachDim) + maxAlignPoint;
            tempMat = nan(tempMatSize);
            for j = 1:length(alignPoint)
                startPoint = maxAlignPoint-alignPoint(j)+1;
                tempMat(startPoint:startPoint+size(tempField,attachDim)-1,j) = tempField(:,j);

            end
            probeData.(tempFieldNames{i}) = tempMat;
        end        
    end
end