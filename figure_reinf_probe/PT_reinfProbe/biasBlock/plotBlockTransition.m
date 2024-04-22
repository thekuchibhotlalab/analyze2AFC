function plotBlockTransition(mouse,varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('transPlotXlim', 50);
    p.parse(varargin{:});

    if ~isobject(mouse); return; end
    if strcmp(class(mouse),'wheel2AFC')
        trialBin = mouse.ops.learningCurveBin;
        mat.action = smoothdata(mouse.behav.action==1,'movmean',trialBin); 
        mat.stimulus = smoothdata(mouse.behav.stimulus==1,'movmean',trialBin); 
        mat.bias = mouse.behav.bias/2+0.5;
        mat.modelBias = mouse.behav.modelBias/2+0.5;
        mat.trans = mouse.biasBlock.trans;
        mat.transID = mouse.biasBlock.transID;
    
        makeSinglePlot({'U2L'},mat,[mouse.ops.mouse ' ']);
        makeSinglePlot({'U2R'},mat,[mouse.ops.mouse ' ']);
        makeSinglePlot({'L2U'},mat,[mouse.ops.mouse ' ']);
        makeSinglePlot({'R2U'},mat,[mouse.ops.mouse ' ']);
    elseif strcmp(class(mouse),'wheel2AFCmega')
        trialBin = mouse.mouseCell{1}.ops.learningCurveBin;
        tempTrans= getProp(mouse,'biasBlock','field','trans'); 
        tempTransID= getProp(mouse,'biasBlock','field','transID');
        tempTrialNum = cumsum([0 cellfun(@length,getProp(mouse,'behav','field','day'))]);
        
        action = getProp(mouse,'behav','field','action');
        stimulus = getProp(mouse,'behav','field','stimulus');
        bias = getProp(mouse,'behav','field','bias');
        responseType = getProp(mouse,'behav','field','responseType');
        modelBias = getProp(mouse,'behav','field','modelBias');
        
        mat.trans = [];mat.transID = {};
        mat.action = []; mat.stimulus = []; mat.bias = []; mat.modelBias = []; mat.responseType = [];
        for i = 1:length(tempTrans)
            mat.trans = cat(1,mat.trans,tempTrans{i} + tempTrialNum(i));
            mat.transID = cat(1,mat.transID,tempTransID{i});
            mat.action = cat(1,mat.action,action{i});
            mat.stimulus = cat(1,mat.stimulus,stimulus{i});
            mat.bias = cat(1,mat.bias,bias{i});
            mat.responseType = cat(1,mat.responseType,responseType{i});
            mat.modelBias = cat(1,mat.modelBias,modelBias{i});
        end
        
        mat.rewardRateL = smoothdata(mat.responseType==1 & mat.action==1,'movmean',trialBin)./smoothdata(mat.action==1,'movmean',trialBin); 
        mat.rewardRateR = smoothdata(mat.responseType==1 & mat.action==2,'movmean',trialBin)./smoothdata(mat.action==2,'movmean',trialBin); 

        mat.action = smoothdata(mat.action==1,'movmean',trialBin);
        mat.stimulus = smoothdata(mat.stimulus==1,'movmean',trialBin); 
        mat.rewardRate = smoothdata(mat.responseType==1,'movmean',trialBin); 
        mat.bias = mat.bias/2+0.5;
        mat.modelBias = mat.modelBias/2+0.5;
        
        blockL = getProp(mouse,'biasBlock','field','blockL');
        lenL = cellfun(@(x)(x.len),blockL,'UniformOutput',false);
        startL = cellfun(@(x)(x.start),blockL,'UniformOutput',false);
        endL = cellfun(@(x)(x.end),blockL,'UniformOutput',false);
        nL_first = cellfun(@(x,y)(sum(x<y/2)),endL,num2cell(mouse.nTrials));
        nL_last = cellfun(@(x,y)(sum(x>y/2)),endL,num2cell(mouse.nTrials));
        nL = cellfun(@(x)(length(x.len)),blockL);
        
        blockR = getProp(mouse,'biasBlock','field','blockR');
        lenR = cellfun(@(x)(x.len),blockR,'UniformOutput',false);
        startR = cellfun(@(x)(x.start),blockR,'UniformOutput',false);
        endR = cellfun(@(x)(x.end),blockR,'UniformOutput',false);
        nR_first = cellfun(@(x,y)(sum(x<y/2)),endR,num2cell(mouse.nTrials));
        nR_last = cellfun(@(x,y)(sum(x>y/2)),endR,num2cell(mouse.nTrials));
        nR = cellfun(@(x)(length(x.len)),blockR);
        
        blockU = getProp(mouse,'biasBlock','field','blockU');
        lenU = cellfun(@(x)(x.len),blockU,'UniformOutput',false);
        startU = cellfun(@(x)(x.start),blockU,'UniformOutput',false);
        endU = cellfun(@(x)(x.end),blockU,'UniformOutput',false);
        nU_first = cellfun(@(x,y)(sum(x<y/2)),endU,num2cell(mouse.nTrials));
        nU_last = cellfun(@(x,y)(sum(x>y/2)),endU,num2cell(mouse.nTrials));
        nU = cellfun(@(x)(length(x.len)),blockU);
        
        for i = 1:length(mouse.nTrials)
            if mouse.nTrials(i)<2000
                nL_first(i) = length(endL{i});
                nR_first(i) = length(endR{i});
                nU_first(i) = length(endU{i});
                nL_last(i) = 0; nR_last(i) = 0; nU_last(i) = 0;

            else
                tempThre = mouse.nTrials(i)/2;
                nL_first(i) = sum((endL{i}+startL{i})/2<tempThre);
                nR_first(i) = sum((endR{i}+startR{i})/2<tempThre);
                nU_first(i) = sum((endU{i}+startU{i})/2<tempThre);
                nL_last(i) = sum((endL{i}+startL{i})/2>=tempThre);
                nR_last(i) = sum((endR{i}+startR{i})/2>=tempThre);
                nU_last(i) = sum((endU{i}+startU{i})/2>=tempThre);
                %disp(mouse.mouseCell{i}.behav.acc(round(mouse.nTrials(i)/2)))
            end
        end


        y1 = [nR; nU; nL;]; y1 = y1./repmat(sum(y1,1),size(y1,1),1); [~,sortIdx1] = sort(y1(1,:),'descend');
        lenLsum = cellfun(@sum,lenL); lenRsum = cellfun(@sum,lenR); lenUsum = cellfun(@sum,lenU);
        y2 = [lenRsum; lenUsum; lenLsum;]; y2 = y2./repmat(sum(y2,1),size(y2,1),1); [~,sortIdx2] = sort(y2(1,:),'descend');

        figure; subplot(1,2,1); bar(y1(:,sortIdx1)','stacked')
        subplot(1,2,2); bar(y2(:,sortIdx2)','stacked')
        
        y1 = [nR_first; nU_first; nL_first;]; y1 = y1./repmat(sum(y1,1),size(y1,1),1); [~,sortIdx1] = sort(y1(1,:),'descend');
        y2 = [nR_last; nU_last; nL_last;]; y2 = y2./repmat(sum(y2,1),size(y2,1),1); [~,sortIdx2] = sort(y2(1,:),'descend');

        figure; subplot(1,2,1); bar(y1(:,sortIdx1)','stacked'); title('proportion of blocks in first half of trials')
        subplot(1,2,2); bar(y2(:,sortIdx2)','stacked');  title('proportion of blocks in last half of trials')
        
        lenPrefer = {}; lenUnPrefer = {}; nPrefer = []; nUnPrefer = []; 
        for i = 1:length(nL)
            lenPrefer{i} = []; lenUnPrefer{i} = []; 
            if sum(lenL{i}) > sum(lenR{i})
                lenPrefer{i} = [lenPrefer{i} lenL{i}];
                lenUnPrefer{i} = [lenUnPrefer{i} lenR{i}];
                nPrefer = [nPrefer nL(i)];
                nUnPrefer = [nUnPrefer nR(i)];
            else
                lenPrefer{i} = [lenPrefer{i} lenR{i}];
                lenUnPrefer{i} = [lenUnPrefer{i} lenL{i}];
                nPrefer = [nPrefer nR(i)];
                nUnPrefer = [nUnPrefer nL(i)];          
            end
        end
        lenPreferFlat = fn_cell2mat(lenPrefer,2); lenUnPreferFlat = fn_cell2mat(lenUnPrefer,2); 

        figure; xAxis = 25 * 2.^(-1.25:0.5:6.25); 
        h = fn_plotHistLine(cell2mat(lenU),'histCountArgIn',{xAxis,'Normalization','count'});
        set(h,'Color',[matlabColors(3) 0.9]); set(h,'LineWidth',2)
        set(gca, 'XScale', 'log'); xlim([xAxis(1) xAxis(end)]); xticks(25 * 2.^(-1:1:6))
        title('Distribution of unbiased block length')
        xlabel('Length of unbiased blocks'); ylabel('Number of blocks')
        

        % PLOT LENGTH OF PREFERRED AND UNPREFERRED BLOCKS
        mouseRand = rand(1,1000000); mouseRand(mouseRand<0.5) = -1; mouseRand(mouseRand>0.5) = 1;
        biasRand = smoothdata(mouseRand,'movmean',mouse.mouseCell{1}.ops.learningCurveBin);
        [~, ~, ~,bareaL] = fn_getBlockOnOff(biasRand>0,biasRand);
        [~, ~, ~,bareaR] = fn_getBlockOnOff(biasRand<0,biasRand);
        barea = [bareaL bareaR]; areaThre = prctile(barea,95);

        motorlenL = {}; motorlenR = {};
        for i = 1:mouse.nMouse
            a = nanmean(abs(bias{i}));
            biasW = -log(1/(0.5+a/2)-1);
            %biasW = log(2*a/(1-a)+sqrt((4*a*a+4+4*a)/(1-a)^2)/2);
            [tempBias,~] = fn_generateChoiceFromAgent(mouse.nTrials(i),[], biasW);
            motorlenL{i} = []; motorlenR{i} = [];
            
            for j = 1:size(tempBias,1)
                [~, ~, idxL,areaL] = fn_getBlockOnOff(tempBias(j,:)>0,tempBias(j,:));
                [~, ~, idxR,areaR] = fn_getBlockOnOff(tempBias(j,:)<0,tempBias(j,:));
                selFlagL = (areaL > areaThre); tempLen = cellfun(@length,idxL); 
                motorlenL{i} = [motorlenL{i} tempLen(selFlagL)];
                selFlagR = (areaR > areaThre); tempLen = cellfun(@length,idxR); 
                motorlenR{i} = [motorlenR{i} tempLen(selFlagR)];
            end
        end
        motorlenLflat = fn_cell2mat(motorlenL,2); motorlenRflat = fn_cell2mat(motorlenR,2);

        xAxis = 25 * 2.^(-0.75:0.5:7.25); 
        xMid = 25 * 2.^(-0.5:0.5:7); 
        
        [nML,~,~] = histcounts(motorlenLflat,xAxis, 'Normalization','count'); nML = nML/200; 
        [nMR,~,~] = histcounts(motorlenRflat,xAxis, 'Normalization','count'); nMR = nMR/200; 
        [nUnP,~,~] = histcounts(lenUnPreferFlat,xAxis, 'Normalization','count');
        [nP,edges,~] = histcounts(lenPreferFlat,xAxis, 'Normalization','count');

        nML = cat(1,nML,nML); nML = nML(:); nML = [0; nML; 0];
        nMR = cat(1,nMR,nMR); nMR = nMR(:); nMR = [0; nMR; 0];
        nUnP = cat(1,nUnP,nUnP); nUnP = nUnP(:); nUnP = [0; nUnP; 0];
        nP = cat(1,nP,nP); nP = nP(:); nP = [0; nP; 0];
        edges = cat(1,edges,edges); edges = edges(:);

        figure; hold on;
        plot([0 0],[0 25],'Color',[0.6 0.6 0.6],'LineWidth',0.7)
        plot([flipud(-log(edges)+log(edges(1))); log(edges)-log(edges(1))] ,[flipud(nML); nMR],...
            'Color',[0.6 0.6 0.6],'LineWidth',2);
        plot([flipud(-log(edges)+log(edges(1))); log(edges)-log(edges(1))] ,[flipud(nUnP); nP],...
            'Color',fn_wheelColorsPT('bias'),'LineWidth',2);
        xticks([-log(xMid(end:-2:2))+log(edges(1)) log(xMid(2:2:end))-log(edges(1))]); 
        xticklabels([xMid(end:-2:2) xMid(2:2:end)]); 
        xlim([-log(edges(end))+log(edges(1)), log(edges(end))-log(edges(1))]); ylim([0 25])

        figure; subplot(2,3,1);hold on;  
        h2 = fn_plotHistLine(motorlenLflat,'histCountArgIn',{xAxis,'Normalization','count'},'yaxisNormalize',200);
        set(h2,'Color',[0.6 0.6 0.6]); set(h2,'LineWidth',2);
        h1 = fn_plotHistLine(lenUnPreferFlat,'histCountArgIn',{xAxis,'Normalization','count'});
        set(h1,'Color',fn_wheelColorsPT('bias')); set(h1,'LineWidth',2)
        legend('Motor biased agent','Mouse behavior');  ylim([0 15]);
        set(gca,'XScale','log'); xlim([xAxis(1) xAxis(end)]); xticks(25 * 2.^(0:1:6))

        subplot(2,3,2);hold on;   
        x1 = cellfun(@(x,y)(sum(x)/y),lenUnPrefer,num2cell(mouse.nTrials)); x1(isnan(x1)) = 0;
        x2 = cellfun(@(x,y)(sum(x)/y),motorlenL,num2cell(mouse.nTrials))/200; x2(isnan(x2)) = 0;
        fn_plotComparison({x2,x1},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn', {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},...
            'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); ylabel('Proportion of biased trials')
        xlim([0.6 2.4]); xticks(1:2);xticklabels({'Motor biased agent','Mouse behavior'}); ylim([0 0.4])
        

        subplot(2,3,3);hold on; 
        x1 = cellfun(@(x,y)(length(x)/y*1000),lenUnPrefer,num2cell(mouse.nTrials)); x1(isnan(x1)) = 0;
        x2 = cellfun(@(x,y)(length(x)/y*1000),motorlenL,num2cell(mouse.nTrials))/200; x2(isnan(x2)) = 0;
        fn_plotComparison({x2,x1},'paired',true,'compType','errorbarWithDot'); ylabel('# blocks per 1000 trials')
        xlim([0.6 2.4]); xticks(1:2);xticklabels({'Motor biased agent','Mouse behavior'})
       
        subplot(2,3,4);hold on;     
        h2 = fn_plotHistLine(motorlenRflat,'histCountArgIn',{xAxis,'Normalization','count'},'yaxisNormalize',200);
        set(h2,'Color',[0.6 0.6 0.6]); set(h2,'LineWidth',2);
        h1 = fn_plotHistLine(lenPreferFlat,'histCountArgIn',{xAxis,'Normalization','count'});
        set(h1,'Color',fn_wheelColorsPT('bias')); set(h1,'LineWidth',2)
        legend('Motor biased agent','Mouse behavior')
        set(gca,'XScale','log'); xlim([xAxis(1) xAxis(end)]); xticks(25 * 2.^(0:1:6))
        ylabel('# of biased blocks, unpreffered direction'); ylim([0 25])

        subplot(2,3,5);hold on;   
        x1 = cellfun(@(x,y)(sum(x)/y),lenPrefer,num2cell(mouse.nTrials)); x1(isnan(x1)) = 0;
        x2 = cellfun(@(x,y)(sum(x)/y),motorlenR,num2cell(mouse.nTrials))/200; x2(isnan(x2)) = 0;
        fn_plotComparison({x2,x1},'paired',true,'dotType','side','compType','errorbarWithDot','scatterArgIn', {5,'MarkerEdgeColor','none','MarkerFaceColor','none'},...
            'errorbarArgIn', {'Color',[0.2 0.2 0.2],'LineWidth',1.5,'LineStyle','none'}); ylabel('Proportion of biased trials')
        xlim([0.6 2.4]); xticks(1:2);xticklabels({'Motor biased agent','Mouse behavior'})
        
        subplot(2,3,6);hold on; 
        x1 = cellfun(@(x,y)(length(x)/y*1000),lenPrefer,num2cell(mouse.nTrials)); x1(isnan(x1)) = 0;
        x2 = cellfun(@(x,y)(length(x)/y*1000),motorlenR,num2cell(mouse.nTrials))/200; x2(isnan(x2)) = 0;
        fn_plotComparison({x2,x1},'paired',true,'compType','errorbarWithDot'); ylabel('# blocks per 1000 trials')
        xlim([0.6 2.4]); xticks(1:2);xticklabels({'Motor biased agent','Mouse behavior'})

        figure;hold on;     
        h1 = fn_plotHistLine(lenPreferFlat,'histCountArgIn',{xAxis,'Normalization','count'});
        set(h1,'Color',[matlabColors(1) 0.9]); set(h1,'LineWidth',2)
        h2 = fn_plotHistLine(lenUnPreferFlat,'histCountArgIn',{xAxis,'Normalization','count'});
        set(h2,'Color',[matlabColors(2) 0.9]); set(h2,'LineWidth',2);
        legend('Preferred','UnPreferred'); ylabel('# of biased blocks, unpreffered direction')
        set(gca,'XScale','log'); xlim([xAxis(1) xAxis(end)]); xticks(25 * 2.^(0:1:6)); 

        figure; 
        subplot(1,2,1);hold on ; cdfplot(lenPreferFlat); cdfplot(lenUnPreferFlat)
        subplot(1,2,2); hold on ; cdfplot(nPrefer); cdfplot(nUnPrefer)
        
        makeSinglePlot({'U2L'},mat,'allMouse ','plotXLim',p.Results.transPlotXlim);
        makeSinglePlot({'U2R'},mat,'allMouse ','plotXLim',p.Results.transPlotXlim);
        makeSinglePlot({'L2U'},mat,'allMouse ','plotXLim',p.Results.transPlotXlim);
        makeSinglePlot({'R2U'},mat,'allMouse ','plotXLim',p.Results.transPlotXlim);
        makeSinglePlot({'L2R'},mat,'allMouse ','plotXLim',p.Results.transPlotXlim);
        makeSinglePlot({'R2L'},mat,'allMouse ','plotXLim',p.Results.transPlotXlim);
    end
        
end

function makeSinglePlot(transType,mat,extraTitle,varargin)
    p = inputParser;
    p.KeepUnmatched = true;
    p.addParameter('plotXLim', 30);
    p.parse(varargin{:});

    transFlag = fn_multistrcmp(mat.transID,transType);
    transIdx = mat.trans(transFlag,:);
    if isempty(transIdx); return; end
    
    preTransMax = max(transIdx(:,2) - transIdx(:,1));
    postTransMax = max(transIdx(:,4) - transIdx(:,2));
    
    actionTransPlot = fn_align2idx(mat.action, transIdx,preTransMax,postTransMax);
    stimulusTransPlot = fn_align2idx(mat.stimulus, transIdx,preTransMax,postTransMax);
    rewardTransPlot = fn_align2idx(mat.rewardRateL, transIdx,preTransMax,postTransMax);
    rewardTransPlot2 = fn_align2idx(mat.rewardRateR, transIdx,preTransMax,postTransMax);
    biasTransPlot = fn_align2idx(mat.bias, transIdx,preTransMax,postTransMax);
    modelTranPlot = fn_align2idx(mat.modelBias, transIdx,preTransMax,postTransMax);
    
    tempStimWindow = 40;
    try
        tempStim = nanmean(stimulusTransPlot(:,preTransMax-tempStimWindow:preTransMax),2);
    catch
        tempStimWindow = 20;
        tempStim = nanmean(stimulusTransPlot(:,preTransMax-tempStimWindow:preTransMax),2);
    end
    tempThre = 0.2;
    removeFlag = tempStim > 0.5+tempThre | tempStim <0.5-tempThre;
    actionTransPlot(removeFlag,:) = [];
    stimulusTransPlot(removeFlag,:) = [];
    biasTransPlot(removeFlag,:) = [];
    modelTranPlot(removeFlag,:) = [];
    %figure; histogram(tempStim,0.3:0.02:0.7);
    plotBin = [-p.Results.plotXLim p.Results.plotXLim];

    % stimulus prob normalize to show -1 to 1 on the bias axis
    normStimulusTransPlot = (stimulusTransPlot -0.5) * 0.6 + 0.5;
    % reward rate normally fluctuate between 0 - 1, normalize to -0.5 - 0.5
    rewardTransPlot = rewardTransPlot + 0.5 - repmat(nanmean(rewardTransPlot(:,preTransMax-30:preTransMax+30),2),[1 size(rewardTransPlot,2)]);
    rewardTransPlot2 = rewardTransPlot2 + 0.5 - repmat(nanmean(rewardTransPlot2(:,preTransMax-30:preTransMax+30),2),[1 size(rewardTransPlot2,2)]);

    disp([transType ' nTransitions ' int2str(size(rewardTransPlot,1))])
    figure; 
    % subplot 1
    plotRedBlue(modelTranPlot,1,[extraTitle strjoin(transType) ' model bias'],plotBin);
    
    % subplot 2
    plotRedBlue(biasTransPlot,2,[extraTitle strjoin(transType) ' behavioral bias'],plotBin);
    yticklabels(-0.6:0.3:0.6);
    
    % subplot 3
    plotRedBlue(actionTransPlot,3,[extraTitle strjoin(transType) ' choice prob'],plotBin);
    % subplot 4
    plotRedBlue(stimulusTransPlot,4,[extraTitle strjoin(transType) ' stim prob'],plotBin);
    %fn_plotMeanErrorbar(-preTransMax:postTransMax,normStimulusTransPlot,matlabColors(5),matlabColors(5,0.2), {'LineWidth',2},{});
    fn_plotMeanErrorbar(-preTransMax:postTransMax,rewardTransPlot,matlabColors(6),matlabColors(6,0.2), {'LineWidth',2},{});
    fn_plotMeanErrorbar(-preTransMax:postTransMax,rewardTransPlot2,matlabColors(1),matlabColors(6,0.2), {'LineWidth',2},{});
    figure; 
    plotRedBlue(rewardTransPlot,1,[extraTitle strjoin(transType) ' L reward rate'],plotBin);
    plotRedBlue(rewardTransPlot2,2,[extraTitle strjoin(transType) ' R reward rate'],plotBin);


    function plotRedBlue(plotData,plotPos,dataTitle,plotBin)
        ylimm = [0.2 0.8]; 
        subplot(2,4,plotPos);
        plotDataTemp = plotData; plotDataTemp(isnan(plotData)) = 0.5;
        imagesc(plotDataTemp); colorbar;
        xlim([preTransMax+plotBin(1) preTransMax+plotBin(2)]); caxis(ylimm)
        colormap(redblue); %colorbar
        xticks(preTransMax+[plotBin(1) 0 plotBin(2)])
        xticklabels(strsplit(int2str([plotBin(1) 0 plotBin(2)])))
        title(dataTitle); xlabel('Trials to Transition'); ylabel('nTransitions')
        axis off; 

        subplot(2,4,4+plotPos); hold on
        plot([0 0],ylimm,'--','LineWidth',2,'Color',[0.8 0.8 0.8]);
        plot([-preTransMax postTransMax],[0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8]);
        fn_plotMeanErrorbar(-preTransMax:postTransMax,plotData,fn_wheelColorsPT('bias'),fn_wheelColorsPT('bias',0.2), {'LineWidth',2},{});
        %plot(-preTransMax:postTransMax,nanmean(plotData,1),'LineWidth',2,'Color',matlabColors(1));
        xlim(plotBin); ylim(ylimm); xticks(linspace(plotBin(1),plotBin(2),5)); 
        yticks(linspace(ylimm(1),ylimm(2),5))


    end
end

function dataTrans = fn_align2idx(data, transIdx, preTransMax, postTransMax)
    nTrans = size(transIdx,1);

    dataTrans = nan(nTrans,preTransMax+postTransMax+1);
    
    for i = 1:nTrans
        tempStartIdx = preTransMax- (transIdx(i,2) - transIdx(i,1));
        tempDataTrans = data(transIdx(i,1):transIdx(i,4));
        dataTrans(i,tempStartIdx+(1:length(tempDataTrans)) ) = tempDataTrans;
    end
    for i = 1:size(dataTrans,1); tempIdx(i) = find(~isnan(dataTrans(i,:)),1); end
    [sortResults,sortIdx] = sort(tempIdx,'ascend');
    dataTrans = dataTrans(sortIdx,:);

    % TAKE OUT THE UNBIASED STATES THAT HAVE LESS THAN 5 REAL VALUES 
    dataTrans(sortResults>=preTransMax,:) = []; 
end

