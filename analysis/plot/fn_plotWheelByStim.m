function fn_plotWheelByStim(mouseMega,exampleAnimal,trialnum,colorSmallDots)
global downSample csizeS csizeL toneOn
downSample = 50; csizeS = 20; csizeL = 20; toneOn = 3000;

tempWheel = mouseMega.mouseCell{exampleAnimal}.behavVar.wheelPos_aligned(trialnum,:);
tempStim = mouseMega.mouseCell{exampleAnimal}.behav.stimulus(trialnum);
tempRT = mouseMega.mouseCell{exampleAnimal}.behav.reactionTime(trialnum);

tempSpeedRaw = diff(tempWheel,1,2); tempSpeed = fn_binSum(tempSpeedRaw, downSample,2); 
tempSpeedAbs = fn_binSum(abs(tempSpeedRaw), downSample,2);
[tempIdx1,tempIdx2] = find(tempSpeedAbs>0 & tempSpeed==0);
for i = 1:length(tempIdx1)
    tempSpeedFirst = tempSpeedRaw(tempIdx1(i),(tempIdx2(i)-1)*downSample+1:(tempIdx2(i))*downSample);
    tempFirstIdx = find(tempSpeedFirst,1,'first'); tempFirstFlag = tempSpeedFirst(tempFirstIdx)>0;
    if tempFirstFlag; tempSpeed(tempIdx1(i),tempIdx2(i)) = 1; else; tempSpeed(tempIdx1(i),tempIdx2(i)) = -1; end
end

% deal with the problem tath 
nSample = round(tempSpeed/downSample);

xStart = 1500; xEnd = 4500; 
xaxis = xStart/downSample:500/downSample:xEnd/downSample;
%{
figure;
subplot(1,2,1);hold on; 
plotWheelAlign2Stim (tempRT, tempSpeed, tempStim==1, 'L')
xlim(xaxis([1 end])); xticks(xaxis); xticklabels(-1.5:0.5:1.5); ylim([0 sum(tempStim==1)])

subplot(1,2,2);hold on; 
plotWheelAlign2Stim (tempRT, tempSpeed, tempStim==2, 'R');
xlim(xaxis([1 end])); xticks(xaxis); xticklabels(-1.5:0.5:1.5); ylim([0 sum(tempStim==2)])

%}
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8);
subplot(1,2,1);hold on; 
tempFlagL = find(tempStim==1); tempFlagL = tempFlagL(1:50);
plotSpeedShiftL = plotWheelAlign2Choice (tempRT, tempSpeed, tempFlagL, 'L',colorSmallDots);
xlim(xaxis([1 end])); xticks(xaxis); xticklabels(-1.5:0.5:1.5); ylim([0 length(tempFlagL)])

subplot(1,2,2);hold on; 
tempFlagR = find(tempStim==2); tempFlagR = tempFlagR(1:50);
plotSpeedShiftR = plotWheelAlign2Choice (tempRT, tempSpeed, tempFlagR, 'R',colorSmallDots);
xlim(xaxis([1 end])); xticks(xaxis); xticklabels(-1.5:0.5:1.5); ylim([0 length(tempFlagR)])


fn_figureSmartDim('hSize',0.25,'widthHeightRatio',0.75); hold on; 
plot([toneOn/downSample toneOn/downSample],[-1 1],'--','Color',[0.8 0.8 0.8],'LineWidth',2);

fn_plotMeanErrorbar(1:size(plotSpeedShiftL,2),plotSpeedShiftL/downSample,...
    fn_wheelColorsPT('L'), fn_wheelColorsPT('L',0.2),{'LineWidth',2}, {});
fn_plotMeanErrorbar(1:size(plotSpeedShiftR,2),plotSpeedShiftR/downSample,...
    fn_wheelColorsPT('R'), fn_wheelColorsPT('R',0.2),{'LineWidth',2}, {});

plot(nanmean(plotSpeedShiftL,1)/downSample,'LineWidth',2,'Color',fn_wheelColorsPT('L')); 
plot(nanmean(plotSpeedShiftR,1)/downSample,'LineWidth',2,'Color',fn_wheelColorsPT('R'));
xlim(xaxis([1 end])); xticks(xaxis); xticklabels(-1.5:0.5:1.5);ylim([-0.6 0.6]); ylabel('wheel speed (deg/ms)')

end


function plotWheelAlign2Stim (RT, speed, selFlag, stimDirection)
global downSample csize toneOn
    plotSpeed =  speed(selFlag,:); [sortRT,sortIdx] = sort(RT(selFlag));plotSpeed = plotSpeed(sortIdx,:);
    plotSpeed(plotSpeed>downSample) = downSample; plotSpeed(plotSpeed<-downSample) = -downSample;
    tempColorAxis = linspace(0,1,downSample);

    [idxL1,idxL2]= find(plotSpeed>=1); [idxR1,idxR2]= find(plotSpeed<=-1);
    sortRT = round(sortRT*1000/downSample); tempArea = [];
    for i = 1:length(sortRT); tempArea = cat(2, tempArea, [ cat(2,sortRT(i)+2000/downSample,sortRT(i)+2000/downSample);cat(2,i-1,i)]); end
    tempArea = cat(2,tempArea,cat(1,2000/downSample*ones(1,length(sortRT)+1),length(sortRT):-1:0));
    fill(tempArea(1,:),tempArea(2,:),[0.8 0.8 0.8],'LineStyle','None');

    % SCATTER OF THE LEFT ACTIONS
    if strcmp(stimDirection,'L'); colorStr = 'correct'; else; colorStr = 'incorrect'; end
    plotColor = zeros(length(idxL1),3);
    for i = 1:length(idxL1)
        tempColorNum = abs(plotSpeed(idxL1(i),idxL2(i)));
        plotColor(i,:) = fn_wheelColorsPT(colorStr,tempColorAxis(tempColorNum));
    end
    scatter(idxL2,idxL1-0.5,csize,plotColor,'square','filled')
    
    % SCATTER OF THE RIGHT ACTIONS 
    if strcmp(stimDirection,'R'); colorStr = 'correct'; else; colorStr = 'incorrect'; end
    plotColor = zeros(length(idxR1),3);
    for i = 1:length(idxR1)
        tempColorNum = abs(plotSpeed(idxR1(i),idxR2(i)));
        plotColor(i,:) = fn_wheelColorsPT(colorStr,tempColorAxis(tempColorNum));
    end
    scatter(idxR2,idxR1-0.5,csize,plotColor,'square','filled')


end


function plotSpeedShift = plotWheelAlign2Choice (RT, speed, selFlag, stimDirection,colorSmallDots)
global downSample csizeS csizeL toneOn

    %preTone = 3000; 
    toneFrame = round(toneOn/downSample); RTshift = 0; 
    plotSpeed = speed(selFlag,:); [sortRT,sortIdx] = sort(RT(selFlag));
    plotSpeed(plotSpeed>downSample) = downSample; plotSpeed(plotSpeed<-downSample) = -downSample;
    sortRT = round(sortRT*1000/downSample,2) + RTshift; plotSpeed = plotSpeed(sortIdx,:);
    tempColorAxis = linspace(0.4,1,downSample);

    [idxL1,idxL2]= find(plotSpeed>=1); [idxR1,idxR2]= find(plotSpeed<=-1);  
    largeCircleFlag = findLargeCircle(idxL1,idxL2,sortRT);

    % Shift the speed with respect to decision time
    plotSpeedShift = plotSpeed; 
    for i = 1:size(plotSpeedShift,1)
        tempShift = -round(sortRT(i));
        plotSpeedShift(i,:) = circshift(plotSpeedShift(i,:),tempShift);
    end

    % Plot the gray area
    tempArea = [];
    for i = 1:length(sortRT); tempArea = cat(2, tempArea, [ cat(2,-sortRT(i)+round(toneOn/downSample),-sortRT(i)+round(toneOn/downSample));cat(2,i-1,i)]); end
    tempArea = cat(2,tempArea,cat(1,round(toneOn/downSample)*ones(1,length(sortRT)+1),length(sortRT):-1:0));
    fill(tempArea(1,:),tempArea(2,:),[0.8 0.8 0.8],'LineStyle','None');

    % SCATTER OF THE LEFT ACTIONS
    if strcmp(stimDirection,'L'); colorStr = 'correct'; else; colorStr = 'incorrect'; end
    plotColor = zeros(length(idxL1),3); colorStr2 = 'ITI'; plotColorITI = zeros(length(idxR1),3);
    for i = 1:length(idxL1)
        tempColorNum = abs(plotSpeed(idxL1(i),idxL2(i)));
        plotColor(i,:) = fn_wheelColorsPT(colorStr,tempColorAxis(tempColorNum));
        plotColorITI(i,:) = fn_wheelColorsPT(colorStr2,tempColorAxis(tempColorNum)); 

    end   
    tempX = idxL2-sortRT(idxL1); tempY = idxL1-0.5; 
    scatter(tempX(largeCircleFlag),tempY(largeCircleFlag),csizeL,plotColor(largeCircleFlag,:),'square','filled') % scatter large dots
    if strcmp('color',colorSmallDots)
        scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColor(~largeCircleFlag,:),'square','filled') % scatter small dots
    else
        %scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColor(~largeCircleFlag,:) ...
        %* [0.3 0.3 0.3; 0.59 0.59 0.59; 0.11 0.11 0.11],'square','filled') % scatter small dots
        scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColorITI(~largeCircleFlag,:)...
        ,'square','filled') % scatter small dots
    end 
    % SCATTER OF THE RIGHT ACTIONS 
    if strcmp(stimDirection,'R'); colorStr = 'correct'; else; colorStr = 'incorrect'; end
    plotColor = zeros(length(idxR1),3); colorStr2 = 'ITI'; plotColorITI = zeros(length(idxR1),3);
    for i = 1:length(idxR1)
        tempColorNum = abs(plotSpeed(idxR1(i),idxR2(i)));
        plotColor(i,:) = fn_wheelColorsPT(colorStr,tempColorAxis(tempColorNum));    
        plotColorITI(i,:) = fn_wheelColorsPT(colorStr2,tempColorAxis(tempColorNum)); 
    end
    

    tempX = idxR2-sortRT(idxR1); tempY = idxR1-0.5; 
    largeCircleFlag = findLargeCircle(idxR1,idxR2,sortRT);
    scatter(tempX(largeCircleFlag),tempY(largeCircleFlag),csizeL,plotColor(largeCircleFlag,:),'filled','s') % scatter large dots
    if strcmp('color',colorSmallDots)
        scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColor(~largeCircleFlag,:),'square','filled') % scatter small dots
    else
        %scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColor(~largeCircleFlag,:) ...
        %* [0.3 0.3 0.3; 0.59 0.59 0.59; 0.11 0.11 0.11],'square','filled') % scatter small dots
        scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColorITI(~largeCircleFlag,:)...
        ,'square','filled') % scatter small dots
    end

    function largeCircleFlag = findLargeCircle(idx1,idx2,sortRT)
        largeCircleFlag = false(size(idx1)); 
        tempTrial = unique(idx1);
        for j = 1:length(tempTrial) % loop through all the trials
            if j==149
                disp('check')
            end
            tempIdx = find(idx1==tempTrial(j)); tempIdx2 = sort(idx2(tempIdx),'ascend');
            tempIdxStart = find(tempIdx2 >= toneFrame,1,'first');  
            tempIdxEnd = find(tempIdx2 <= toneFrame+sortRT(j)+1,1,'last'); endFlag = false; 
            while ~endFlag
                if ~isempty(tempIdxStart) && ~isempty(tempIdxEnd) &&  length(tempIdx2)>tempIdxEnd && tempIdx2(tempIdxEnd)+1 == tempIdx2(tempIdxEnd+1); tempIdxEnd = tempIdxEnd+1;
                else; endFlag = true;end
            end
            try
                tempLargeIdx = ( idx2>=tempIdx2(tempIdxStart) &  idx2<= tempIdx2(tempIdxEnd) ) & idx1 == tempTrial(j);
            catch
                tempLargeIdx = []; %disp('un oh');
            end
            largeCircleFlag(tempLargeIdx) = true;
        end
    end

end