function fn_plotWheelByTrial(mouseMega,exampleAnimal,trialnum,colorSmallDots,varargin)
global downSample csizeS csizeL toneOn
downSample = 50; csizeS = 2; csizeL = 6; toneOn = 3000;


tempWheel = mouseMega.mouseCell{exampleAnimal}.behavVar.wheelPos_aligned(trialnum,:);
tempResp = mouseMega.mouseCell{exampleAnimal}.behav.responseType(trialnum);
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

xStart = 500; xEnd = 5500; 
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
fn_figureSmartDim('hSize',0.25,'widthHeightRatio',1.8);hold on; 
plotWheelAlign2Choice (tempRT, tempSpeed, tempStim,colorSmallDots,varargin{:});
xlim(xaxis([1 end])); xticks(xaxis); xticklabels(-2.5:0.5:2.5); ylim([0 length(tempStim)]); yticks([0 length(tempStim)])

end


function plotSpeedShift = plotWheelAlign2Choice (RT, speed, tempStim,colorSmallDots,varargin)

    p = inputParser;
    p.KeepUnmatched = true;
    % ------ Varargin for superclass
    p.addParameter('probeIdx', []);
    p.parse(varargin{:});
    
    global downSample csizeS csizeL toneOn

    %preTone = 3000; 
    toneFrame = round(toneOn/downSample); RTshift = 1; 
    plotSpeed =  speed; sortRT = RT;
    plotSpeed(plotSpeed>downSample) = downSample; plotSpeed(plotSpeed<-downSample) = -downSample;
    sortRT = round(sortRT*1000/downSample) + RTshift;
    tempColorAxis = linspace(0.4,1.6,downSample);

    [idx1,idx2]= find(plotSpeed~=0); 
    largeCircleFlag = findLargeCircle(idx1,idx2,sortRT);
    plotColor = zeros(length(idx1),3);

    % Shift the speed with respect to decision time
    plotSpeedShift = plotSpeed; 
    for i = 1:size(plotSpeedShift,1)
        tempShift = -sortRT(i);
        plotSpeedShift(i,:) = circshift(plotSpeedShift(i,:),tempShift);
    end

    % COLORING SCHEME FOR TRIAL-BY-TRIAL COLOR
    correctFlag = false(size(idx1)); incorrectFlag = false(size(idx1)); 
    tempTrial = unique(idx1);
    for i = 1:length(tempTrial) % loop through all the trials
        tempIdx = find(idx1==tempTrial(i)); tempIdx2 = idx2(tempIdx);
        selStim = tempStim(i); selStim = -selStim*2+3; % convert sti to -1 and 1
        % START HERE TOM
        for k = 1:length(tempIdx)
            selSpeed = plotSpeed(i,tempIdx2(k)); selDir = sign(selSpeed);
            tempColorNum = abs(selSpeed);
            if selStim*selDir==1; 
                correctFlag(tempIdx(k)) = true; 
                plotColor(tempIdx(k),:) = fn_wheelColorsPT('correct',tempColorAxis(tempColorNum));
            else; 
                plotColor(tempIdx(k),:) = fn_wheelColorsPT('incorrect',tempColorAxis(tempColorNum));
                incorrectFlag(tempIdx(k)) = true; 
            end
        end
    end

    % Plot the gray area
    tempArea = [];
    for i = 1:length(sortRT); tempArea = cat(2, tempArea, [ cat(2,-sortRT(i)+round(toneOn/downSample),-sortRT(i)+round(toneOn/downSample));cat(2,i-1,i)]); end
    tempArea = cat(2,tempArea,cat(1,round(toneOn/downSample)*ones(1,length(sortRT)+1),length(sortRT):-1:0));
    fill(tempArea(1,:),tempArea(2,:),[0.8 0.8 0.8],'LineStyle','None');

    if ~isempty(p.Results.probeIdx)
        probeIdx = [p.Results.probeIdx(1)-1 p.Results.probeIdx];
        probeArea = tempArea; probeFlag = probeArea(2,:);
        probeFlag = fn_findAny(probeFlag,probeIdx);
        probeArea = probeArea(:,probeFlag);
        fill(probeArea(1,:),probeArea(2,:),fn_wheelColorsPT('Probe',0.3),'LineStyle','None');
    end


    tempX = idx2-sortRT(idx1); tempY = idx1-0.5; 
    scatter(tempX(largeCircleFlag),tempY(largeCircleFlag),csizeL,plotColor(largeCircleFlag,:),'filled') % scatter large dots
    if strcmp('color',colorSmallDots)
        scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColor(~largeCircleFlag,:),'filled') % scatter small dots
    else
        scatter(tempX(~largeCircleFlag),tempY(~largeCircleFlag),csizeS,plotColor(~largeCircleFlag,:) ...
        * [0.3 0.3 0.3; 0.59 0.59 0.59; 0.11 0.11 0.11],'filled') % scatter small dots
    end 
    

    function largeCircleFlag = findLargeCircle(idx1,idx2,sortRT)
        largeCircleFlag = false(size(idx1)); 
        tempTrial = unique(idx1);
        for j = 1:length(tempTrial) % loop through all the trials
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