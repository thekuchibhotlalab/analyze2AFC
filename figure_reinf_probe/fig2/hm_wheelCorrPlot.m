%% LOAD DATA
clear; 
mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
opsParam.biasBlockType = 'AUC';
allMouse = fn_getObjPT_bin40(mice,opsParam); 
mouseMega = wheel2AFCmega(allMouse);
probeThre = nan; reinfThre = 0.70; probeTrialBin = 400;
[label,attributesOrig] = fn_loadWheelCluster();
fastLabel = label'; fastLabel(fastLabel<=2) = 0; fastLabel(fastLabel>2) = 1;
labelCell = {};attributesCell= {};tempTrialSum = [0 cumsum(mouseMega.nTrials)];
for i = 1:mouseMega.nMouse
    labelCell{i} = fastLabel(tempTrialSum(i)+1:tempTrialSum(i+1)); 
    attributesCell{i} = attributesOrig(tempTrialSum(i)+1:tempTrialSum(i+1),:);
    allMouse{i}.behav.label = labelCell{i};

    allMouse{i}.behav.onsetTime = attributesCell{i}(:,1); allMouse{i}.behav.onsetSpeed = attributesCell{i}(:,2); allMouse{i}.behav.onsetDist = attributesCell{i}(:,3); 
    allMouse{i}.behav.totalTime = attributesCell{i}(:,4); allMouse{i}.behav.totalSpeed = attributesCell{i}(:,5); allMouse{i}.behav.totalDist = attributesCell{i}(:,6); 
end
mouseMega = wheel2AFCmega(allMouse);

%% Figure 2I -- reaction time
figure; subplot(1,2,1)
comp = {'L','R'}; a = makeCorrPlot(allMouse,'reactionTime',comp,'balanced');
subplot(1,2,2); comp = {'L','R'}; a = makeCorrPlot(allMouse,'reactionTime',comp,'balanced');
ylim([-0.15 0.15]); xlim([-0.5 0.5])
%% Figure 2J -- onset time
figure; subplot(1,2,1)
comp = {'L','R'}; a = makeCorrPlot(allMouse,'onsetTime',comp,'balanced');
subplot(1,2,2); comp = {'L','R'}; a = makeCorrPlot(allMouse,'onsetTime',comp,'balanced');
ylim([-150 150]); xlim([-0.5 0.5])

%% Figure 2M -- average speed
figure; subplot(1,2,1)
comp =  {'L','R'}; a = makeCorrPlot(allMouse,'totalSpeed',comp,'balanced');
subplot(1,2,2); comp = {'L','R'}; a = makeCorrPlot(allMouse,'totalSpeed',comp,'balanced');
ylim([-2 2]); xlim([-0.5 0.5])


%% Figure 2L -- ITI movement plot
comp = {'U'};
[biasDiffL,biasDiffDirL, itemDiffL] = corrPlot(allMouse,'itiMoveL',comp);
[biasDiffR,biasDiffDirR, itemDiffR] = corrPlot(allMouse,'itiMoveR',comp);
[biasDiffLR,biasDiffDirLR, itemDiffLR] = corrPlot(allMouse,'itiMove',comp);
biasDiffDirL = fn_cell2mat(biasDiffDirL,2); itemDiffL = fn_cell2mat(itemDiffL,2);
biasDiffDirR = fn_cell2mat(biasDiffDirR,2); itemDiffR = fn_cell2mat(itemDiffR,2);
biasDiffDirLR = fn_cell2mat(biasDiffDirLR,2); itemDiffLR = fn_cell2mat(itemDiffLR,2);
%figure; subplot(1,2,1);
%[lm] = fn_plotScatterCorr({biasDiffDirL},{itemDiffL});
%xlabel('Increase in leftward bias'); ylabel('Increase of ITI wheel movement in L direction')
%subplot(1,2,2); [lm] = fn_plotScatterCorr({biasDiffDirR},{itemDiffR});
%xlabel('Increase in leftward bias'); ylabel('Increase of ITI wheel movement in R direction')

figure; subplot(1,2,1);
[lm] = fn_plotScatterCorr({biasDiffDirLR},{itemDiffLR},'scatterColors', {[0.6 0.6 0.6]},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]});
xlabel('Increase in left bias'); ylabel('Increase of left ITI movement'); xlim([-1 1]); ylim([-30 30]); yticks(-30:15:30)
tempR = corrcoef([biasDiffDirLR',itemDiffLR'],'rows','complete'); disp([ 'r=' num2str(tempR(2,1)) ' slope=' num2str(lm{1}.Coefficients{2,1})])

subplot(1,2,2);
[lm] = fn_plotScatterCorr({biasDiffDirLR},{itemDiffLR},'scatterColors', {[0.6 0.6 0.6]},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]});
xlabel('Increase in left bias'); ylabel('Increase of left ITI movement'); xlim([-0.5 0.5]); ylim([-10 10])
%% Supplemental figure -- Plot everything in 70 trial-bin (the max bin number to keep half of blocks)

%% Supplemental figure -- change in label
figure; comp = {'L','R'}; b = makeCorrPlot(allMouse,'label',comp,'balanced');

%% Supplemental Figure -- use unbiased as a control
figure; comp = {'U'};
subplot(1,4,1); a = makeCorrPlot(allMouse,'reactionTime',comp,'balanced');
subplot(1,4,2); a = makeCorrPlot(allMouse,'onsetTime',comp,'balanced'); 
subplot(1,4,3); a = makeCorrPlot(allMouse,'totalSpeed',comp,'balanced');
subplot(1,4,4); a = makeCorrPlot(allMouse,'label',comp,'balanced');


%% unincluded Figure -- split into fast and slow trials
comp = {'L','R'};
makeCorrPlotFastSlow(allMouse,'reactionTime',comp)
makeCorrPlotFastSlow(allMouse,'onsetTime',comp); 
makeCorrPlotFastSlow(allMouse,'totalSpeed',comp);
makeCorrPlotFastSlow(allMouse,'label',comp);

%% unincluded Figure -- evaluate directional effect of blocks
comp = {'L'};
a = makeCorrPlot(allMouse,'reactionTime',comp,'');
a = makeCorrPlot(allMouse,'onsetTime',comp,''); 
a = makeCorrPlot(allMouse,'totalSpeed',comp,'');
a = makeCorrPlot(allMouse,'label',comp,'');
comp = {'R'};
a = makeCorrPlot(allMouse,'reactionTime',comp,'');
a = makeCorrPlot(allMouse,'onsetTime',comp,''); 
a = makeCorrPlot(allMouse,'totalSpeed',comp,'');
a = makeCorrPlot(allMouse,'label',comp,'');


%% FUNCTIONS

function itemDiffBalanced = makeCorrPlot(allMouse,itemName,comp,oneFigFlag, colorIn)
    if ~exist('oneFigFlag','var'); oneFigFlag = false; end; if ~exist('colorIn','var'); colorIn = [0.6 0.6 0.6]; end 
    [biasDiff, ~,itemDiff, itemDiffBalanced,itemDiffL,itemDiffR] = corrPlot(allMouse,itemName,comp);

    biasDiff = fn_cell2mat(biasDiff,2); itemDiff = fn_cell2mat(itemDiff,2); itemDiffBalanced = fn_cell2mat(itemDiffBalanced,2);
    itemDiffL = fn_cell2mat(itemDiffL,2); itemDiffR = fn_cell2mat(itemDiffR,2);
    if strcmp(oneFigFlag,'balanced')
        [lm] = fn_plotScatterCorr({biasDiff},{itemDiffBalanced},'scatterColors', {colorIn}, 'scatterArgIn',{'filled','MarkerFaceAlpha',0.2},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]}); 
        title([itemName ', balanced']); xlabel('Change in bias level'); ylabel(['Change in ' itemName]);
        tempR = corrcoef([biasDiff',itemDiffBalanced'],'rows','complete'); disp([ 'r=' num2str(tempR(2,1)) ' slope=' num2str(lm{1}.Coefficients{2,1})])
    elseif strcmp(oneFigFlag,'all')
        [lm] = fn_plotScatterCorr({biasDiff},{itemDiff},'scatterColors', {colorIn}, 'scatterArgIn',{'filled','MarkerFaceAlpha',0.2},'corrPlotArgIn',{'LineWidth',1.5,'Color',[0 0 0]}); 
        title([itemName ', all trial']); xlabel('Change in bias level'); ylabel(['Change in ' itemName]);
    else

        figure; subplot(2,2,1); [lm] = fn_plotScatterCorr({biasDiff},{itemDiff}); title([itemName ', all']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
        subplot(2,2,2); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffBalanced}); title([itemName ', balanced']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
        subplot(2,2,3); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffL}); title([itemName ', L']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
        subplot(2,2,4); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffR}); title([itemName ', R']);
        xlabel('Change in bias level'); ylabel(['Change in ' itemName])
    end 
end

function  makeCorrPlotFastSlow(allMouse,itemName,comp)
    [biasDiff, ~,~, ~,~,~,itemDiffBalancedFast,itemDiffBalancedSlow] = corrPlot(allMouse,itemName,comp);

    biasDiff = fn_cell2mat(biasDiff,2); 
    itemDiffBalancedFast = fn_cell2mat(itemDiffBalancedFast,2);
    itemDiffBalancedSlow = fn_cell2mat(itemDiffBalancedSlow,2);
   
    figure; subplot(1,2,1); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffBalancedFast}); title([itemName ', fast']);
    xlabel('Change in bias level'); ylabel(['Change in ' itemName])
    subplot(1,2,2); [lm] = fn_plotScatterCorr({biasDiff},{itemDiffBalancedSlow}); title([itemName ', slow']);
    xlabel('Change in bias level'); ylabel(['Change in ' itemName])

end



function  [biasDiff,biasDiffDir, itemDiff,itemDiffBalanced,itemDiffL,itemDiffR,itemDiffBalancedSlow,itemDiffBalancedFast] = corrPlot(allMouse,item,comp)
nMouse = length(allMouse); 
allTbl_ctxt = cell(1,nMouse); countL = cell(1,nMouse); countR = cell(1,nMouse); countBlock = cell(1,nMouse);
totalBlockRemoved = 0; totalBlock = 0;
for i = 1:nMouse
    allTbl_ctxt{i} = allMouse{i}.behav; 
    nTrials = size(allMouse{i}.behav,1); 

    binSize = 70; moveBin = 10;  
    %blockStartFlag = contains(allMouse{i}.biasBlock.transID,comp);
    blockStart = []; blockEnd = []; 
    if contains('L',comp)
        blockStart = [blockStart allMouse{i}.biasBlock.blockL.start];
        blockEnd = [blockEnd allMouse{i}.biasBlock.blockL.end];
    end
    if contains('R',comp)
        blockStart = [blockStart allMouse{i}.biasBlock.blockR.start];
        blockEnd = [blockEnd allMouse{i}.biasBlock.blockR.end];
    end
    if contains('U',comp)
        blockStart = [blockStart allMouse{i}.biasBlock.blockU.start];
        blockEnd = [blockEnd allMouse{i}.biasBlock.blockU.end];
    end
    %blockStart = (allMouse{i}.biasBlock.trans(blockStartFlag,2)); 
    %blockEnd = (allMouse{i}.biasBlock.trans(blockStartFlag,4));

    biasDiffDir{i} = []; biasDiff{i} = []; itemDiff{i} = []; itemDiffBalanced{i} = []; itemDiffL{i} = [];itemDiffR{i} = [];
    itemDiffL{i} = [];itemDiffR{i} = []; itemDiffBalancedSlow{i} = []; itemDiffBalancedFast{i} = [];
    
    for j =1:length(blockStart)
        totalBlock = totalBlock+1;
        tempLen = blockEnd(j) - blockStart(j)+1;
        if tempLen >= binSize*2
            tempCrop = ceil(mod(tempLen,moveBin)/2); 
            tempStart = blockStart(j)+tempCrop; tempEnd = tempStart+floor(tempLen/moveBin)*moveBin-1;
            tempBin = tempStart:moveBin:tempEnd; tempMaxIdx = find(tempBin == (max(tempBin) - binSize*2));
            for k = 1:tempMaxIdx
                idx1 = tempBin(k):tempBin(k)+binSize-1; idx2 = tempBin(k)+binSize:tempBin(k)+binSize*2-1;
                [tempBias1,~,~,~,~,~] = fn_getAccBias(allMouse{i}.behav.stimulus(idx1),...
                    allMouse{i}.behav.responseType(idx1)==1);
                [tempBias2,~,~,~,~,~] = fn_getAccBias(allMouse{i}.behav.stimulus(idx2),...
                    allMouse{i}.behav.responseType(idx2)==1);
                tempBiasDiff = abs(tempBias2) - abs(tempBias1);
                tempBiasDiffDir = (tempBias2) - (tempBias1);
                tempItemDiff = nanmean(allMouse{i}.behav.(item)(idx2))-nanmean(allMouse{i}.behav.(item)(idx1));
                
                tempTable1 = allMouse{i}.behav(idx1,:); tempTable2 = allMouse{i}.behav(idx2,:);
                tempItem1 = nanmean(tempTable1.(item)(tempTable1.action==1))/2 + nanmean(tempTable1.(item)(tempTable1.action==2))/2;
                tempItem2 = nanmean(tempTable2.(item)(tempTable2.action==1))/2 + nanmean(tempTable2.(item)(tempTable2.action==2))/2;
                tempItemDiffBalanced = tempItem2 - tempItem1;
                tempItemDiffL = nanmean(tempTable2.(item)(tempTable2.action==1)) - nanmean(tempTable1.(item)(tempTable1.action==1));
                tempItemDiffR = nanmean(tempTable2.(item)(tempTable2.action==2)) - nanmean(tempTable1.(item)(tempTable1.action==2));


                tempItem1Slow = nanmean(tempTable1.(item)(tempTable1.action==1 & tempTable1.label==0))/2 +...
                    nanmean(tempTable1.(item)(tempTable1.action==2 & tempTable1.label==0))/2;
                tempItem1Fast = nanmean(tempTable1.(item)(tempTable1.action==1 & tempTable1.label==1))/2 +...
                    nanmean(tempTable1.(item)(tempTable1.action==2 & tempTable1.label==1))/2;
                tempItem2Slow = nanmean(tempTable2.(item)(tempTable2.action==1 & tempTable2.label==0))/2 +...
                    nanmean(tempTable2.(item)(tempTable2.action==2 & tempTable2.label==0))/2;
                tempItem2Fast = nanmean(tempTable2.(item)(tempTable2.action==1 & tempTable2.label==1))/2 +...
                    nanmean(tempTable2.(item)(tempTable2.action==2 & tempTable2.label==1))/2;

                tempItemDiffBalancedFast = tempItem2Fast - tempItem1Fast;
                tempItemDiffBalancedSlow = tempItem2Slow - tempItem1Slow;

                biasDiff{i} = [biasDiff{i} tempBiasDiff];  biasDiffDir{i} = [biasDiffDir{i} tempBiasDiffDir]; itemDiff{i} = [itemDiff{i} tempItemDiff]; 
                itemDiffL{i} = [itemDiffL{i} tempItemDiffL];itemDiffR{i} = [itemDiffR{i} tempItemDiffR];
                itemDiffBalanced{i} = [itemDiffBalanced{i} tempItemDiffBalanced];
                itemDiffBalancedFast{i} = [itemDiffBalancedFast{i} tempItemDiffBalancedFast];
                itemDiffBalancedSlow{i} = [itemDiffBalancedSlow{i} tempItemDiffBalancedSlow];
            end

        else
            totalBlockRemoved = totalBlockRemoved+1;
        end 
    end

end
disp([totalBlock totalBlockRemoved])

end

