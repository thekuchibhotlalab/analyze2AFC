function outStruct = fn_preprocessWheel(wheelTrace, wheelTraceCrop,timeBef,timeAft,behav,ops)
    
    nTrials = size(wheelTraceCrop,1);
    if contains(ops.mouse,{'zz113','zz115'}); tempMax = 100;
    else; tempMax = 35; end
    tempDiv = 35/tempMax; % convert to deg/ms
    % ------------------------- PRETRIAL SAMPLING ---------------------------------------------   
    tempSpeed = diff(wheelTraceCrop,1,2); preMoveTemp = tempSpeed(:,1:timeBef);
    outStruct.preMoveAbs = nanmean(abs(preMoveTemp),2)*tempDiv*1000; outStruct.preMove = nanmean(preMoveTemp,2)*tempDiv*1000; 
    % note: use nanmean to account for difference in premove time, *1000 for sec

    % separate for left and right
    preMoveL = nan(nTrials,1);preMoveR = nan(nTrials,1);
    for i = 1:nTrials
        preMoveTemp = tempSpeed(i,1:timeBef);
        flagL = preMoveTemp>0; flagR = preMoveTemp<0;
        preMoveL(i) = nanmean(preMoveTemp(flagL));
        preMoveR(i) = nanmean(preMoveTemp(flagR));
    end
    outStruct.preMoveL = preMoveL*tempDiv*1000; outStruct.preMoveR = preMoveR*tempDiv*1000; 

    % ------------------------- ITI SAMPLING ---------------------------------------------   
    % Take all the time in a trial that is 2sec after response is reached
    tempRT = round(behav.reactionTime*1000); sampleAfterRT = 2000; 
    itiMoveL = nan(nTrials,1); itiMoveR = nan(nTrials,1); itiMove = nan(nTrials,1); itiMoveAbs = nan(nTrials,1);

    for i = 1:nTrials
        itiMoveTemp = diff(wheelTrace{i},1);
        itiMoveTemp = itiMoveTemp(timeBef+tempRT(i)+sampleAfterRT:end);
        itiMoveAbs(i) = nanmean(abs(itiMoveTemp)); itiMove(i) = nanmean(itiMoveTemp);
        flagL = itiMoveTemp>0; itiMoveL(i) = nanmean(itiMoveTemp(flagL));
        flagR = itiMoveTemp<0; itiMoveR(i) = nanmean(itiMoveTemp(flagR));
    end
    outStruct.itiMoveL = [nan;itiMoveL(1:end-1)*tempDiv*1000]; outStruct.itiMoveR = [nan;itiMoveR(1:end-1)*tempDiv*1000]; 
    outStruct.itiMove = [nan;itiMove(1:end-1)*tempDiv*1000]; outStruct.itiMoveAbs = [nan;itiMoveAbs(1:end-1)*tempDiv*1000]; 
    dayChange = diff([0; behav.day],1)==1;
    outStruct.itiMoveL(dayChange)= nan; outStruct.itiMoveR(dayChange)= nan;
    outStruct.itiMove(dayChange)= nan; outStruct.itiMoveAbs(dayChange)= nan;
    choiceFlag = behav.action; choiceFlag(choiceFlag == 2) = -1; 
    outStruct.itiChoiceDir = outStruct.itiMove .* choiceFlag;


    % ------------------------- TOTAL CHOICE PERIOD MOVEMENT ---------------------------------------------   
    timeAroundChoice = 500;  choiceMove = nan(nTrials,1); choiceMoveAbs = nan(nTrials,1); choiceMoveCorr = nan(nTrials,1);
    stimFlag = behav.stimulus; stimFlag(stimFlag == 2) = -1; 
    for i = 1:nTrials
        if ~isempty(wheelTrace{i})
            choiceTemp = diff(wheelTrace{i},1);
            choiceTemp = choiceTemp(timeBef+tempRT(i)-timeAroundChoice:timeBef+tempRT(i)+timeAroundChoice);
            choiceMove(i) = nansum((choiceTemp))*tempDiv; choiceMoveAbs(i) = nansum(abs(choiceTemp))*tempDiv; 
            choiceMoveCorr(i) = nansum((choiceTemp))*stimFlag(i)*tempDiv; 
        else 
            choiceMove(i) = nan; choiceMoveAbs(i) = nan; choiceMoveCorr(i) = nan; 
        end
    end
    outStruct.choiceMove = choiceMove; outStruct.choiceMoveAbs = choiceMoveAbs; outStruct.choiceMoveCorr = choiceMoveCorr;

    % ------------------------- PREPARE FOR DOWNSAMPLE FOR ATTRIBUTES ---------------------------------------------   
    downSampleRate = 20; tempTimeDownSampleIdx = timeBef+1:downSampleRate:timeBef+2521;
    rt = floor(behav.reactionTime*1000/downSampleRate); rt = min(rt,floor(2500/downSampleRate));
    
    wheelPos = wheelTraceCrop(:,tempTimeDownSampleIdx); % crop the time to 0 to 2.5 sec
    % ------------------------- DOWNSAMPLE FOR ATTRIBUTES ---------------------------------------------   
    % normalize the max position of the wheel to 1; downsample every 50ms 
    wheelPosThre = wheelPos * tempDiv;
    %try for i = 1:nTrials; wheelPosThre(i,rt(i)+2:end) = wheelPosThre(i,rt(i)+1); end
    %catch; error('ERROR -- RT>2.5s found'); end
    %[wheelPos, badFlag] = fn_removeBadTrials(wheelPos);
    %wheelPosThre(wheelPosThre>1) = 1; wheelPosThre(wheelPosThre<-1) = -1;
    % fix the wheel movement after the last bout
    tempWheelChange = [zeros(size(wheelPosThre,1),1) diff(wheelPosThre,1,2)];
    for i = 1:nTrials
        % find one-block or multiple block movements
        [onsetIdx, offIdx,blockIdx] = fn_getBlockOnOff(tempWheelChange(i,:)~=0);
        % find last bout and fix wheel position after last bout
        lastBoutIdx = find(onsetIdx<rt(i)+1,1,'last');
        if ~isempty(lastBoutIdx); wheelPosThre(i,offIdx(lastBoutIdx)+1:end) = wheelPosThre(i,offIdx(lastBoutIdx)); end
    end
    % ------------------------- GET SIX ATTRIBUTES ---------------------------------------------  
    tempWheelChange = [zeros(size(wheelPosThre,1),1) diff(wheelPosThre,1,2)];
    wheelOneBlockIdx = false(nTrials,1); wheelMultiBlockIdx = false(nTrials,1);
    onsetDist = nan(nTrials,1); onsetSpeed = nan(nTrials,1);totalDist = nan(nTrials,1);  overshoot = nan(nTrials,1);
    onsetTime = nan(nTrials,1); totalTime = nan(nTrials,1); totalSpeed = nan(nTrials,1);  onsetFrame = nan(nTrials,1);
    for i = 1:nTrials
        % find one-block or multiple block movements
        [onsetIdx, offIdx,blockIdx] = fn_getBlockOnOff(tempWheelChange(i,:)~=0);
        % correct for short blocks
        tempSelFlag = [];
        for j = 1:length(blockIdx); tempSelFlag(j) = sum(abs(tempWheelChange(i,blockIdx{j}))); end
        moveBlockThre = 0.2; onsetIdx(tempSelFlag< moveBlockThre) = [];
        offIdx(tempSelFlag< moveBlockThre) = []; blockIdx(tempSelFlag< moveBlockThre) = [];
        if ~isempty(blockIdx)
            onsetDist(i) = sum(abs(tempWheelChange(i,blockIdx{1}))) ; 
            onsetSpeed(i) = onsetDist(i) / length(blockIdx{1}) * 1000/downSampleRate ; % convert to ms for speed
            totalDist(i) = sum(abs(tempWheelChange(i,:))) ; 
            
            totalTime(i) = (offIdx(end) - onsetIdx(1)) * downSampleRate; % convert to ms
            totalSpeed(i) = totalDist(i)/totalTime(i) * 1000/downSampleRate;
            onsetTime(i) = onsetIdx(1) * downSampleRate; onsetFrame(i) = onsetIdx(1);
            overshoot(i) = max(abs(wheelPosThre(i,:)));
        end 
            
        
        if length(onsetIdx)>1; wheelMultiBlockIdx(i) = 1;
        else; wheelOneBlockIdx(i) = 1; end      

    end
    % Construct matrix for clustering
    outStruct.onsetDistN = onsetDist; outStruct.onsetSpeedN = onsetSpeed; outStruct.onsetTimeN = onsetTime; 
    outStruct.totalDistN = totalDist; outStruct.totalSpeedN = totalSpeed; outStruct.totalTimeN = totalTime; 
    outStruct.overshoot = overshoot;

end