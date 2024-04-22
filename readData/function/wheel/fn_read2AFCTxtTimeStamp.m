function [wheelPos_aligned,trialData] = fn_read2AFCTxtTimeStamp(txtFilename)

soundOnFlag = false; hit = [];
wheelSoundOn = {}; wheelPreSound = {};  
tempSound = {};  tempPreSound = {}; 

if ~isfile(txtFilename)
    error(['Error - txtFilename ' txtFilename ' does not exist!']); 
end

txtData = splitlines(fileread(txtFilename));
trialStartIdx = find(contains(txtData,'Trial')); trialStartIdx = [trialStartIdx;length(txtData)];
nTrials = sum(contains(txtData,'Trial'));
trialData = struct(); % trialdata records the onset of each trial, sounds and wheel timestamp within that trial

for i = 1:length(trialStartIdx)-1
    tempTrialData = txtData(trialStartIdx(i):(trialStartIdx(i+1)-1));
    
    % find sound time points
    tempSoundIdx = find(contains(tempTrialData,'Sound'));
    [trialData.timeStart(i), trialData.wheelStart(i) ] = findTimeWheel(tempTrialData{tempSoundIdx(1)});
    for j = 1:length(tempSoundIdx); [trialData.sound{i}(j), ~] = findTimeWheel(tempTrialData{tempSoundIdx(j)}); end
    trialData.sound{i} = trialData.sound{i} -  trialData.timeStart(i);
    % find choice time points 
    tempChoiceIdx = find(contains(tempTrialData,'choice'));
    if ~isempty(tempChoiceIdx)
        [trialData.choiceTime(i), ~] = findTimeWheel(tempTrialData{tempChoiceIdx}); 
        trialData.choiceTime(i) = trialData.choiceTime(i) -  trialData.timeStart(i);
    else 
        trialData.choiceTime(i) = nan;
    end
    % find wheel time points
    tempWheelIdx = contains(tempTrialData,'Time'); tempWheelIdx(tempSoundIdx) = 0; tempWheelIdx = find(tempWheelIdx);
    trialData.wheelTime{i} = zeros(1,length(tempWheelIdx)); trialData.wheelPos{i} = zeros(1,length(tempWheelIdx));
    for j = 1:length(tempWheelIdx); [trialData.wheelTime{i}(j), trialData.wheelPos{i}(j)] = findTimeWheel(tempTrialData{tempWheelIdx(j)}); end
    trialData.wheelPos{i} = trialData.wheelPos{i} - trialData.wheelStart(i);
    trialData.wheelTime{i} = trialData.wheelTime{i} - trialData.timeStart(i);
    % define wheel off as either sounds off or choice made
    if isnan(trialData.choiceTime(i)); trialData.timeEnd(i) = trialData.sound{i}(end);
    else; trialData.timeEnd(i) = trialData.choiceTime(i);
    end 
end
% extact preTrial period wheel movement
preTrialThre = 3000; % 2 seconds 

% extract the pretrial for the first trial
tempTrialData = txtData(1:(trialStartIdx(1)-1));
tempWheelIdx = find(contains(tempTrialData,'Time'));
tempWheelTime = zeros(1,length(tempWheelIdx)); tempWheelPos = zeros(1,length(tempWheelIdx));
for j = 1:length(tempWheelIdx); [tempWheelTime(j), tempWheelPos(j)] = findTimeWheel(tempTrialData{tempWheelIdx(j)}); end
preTrialTimeThre = (trialData.timeStart(1) - preTrialThre); preTrialIdx = tempWheelTime >= preTrialTimeThre;
trialData.wheelTime_preTrial{1} = tempWheelTime(preTrialIdx) - trialData.timeStart(1);
trialData.wheelPos_preTrial{1} = tempWheelPos(preTrialIdx) - trialData.wheelStart(1);
% extract the pretrial for the later trials
for i = 1:(nTrials-1)
    preTrialTimeThre = (trialData.timeStart(i+1) - trialData.timeStart(i) - preTrialThre);
    preTrialIdx = trialData.wheelTime{i} >= preTrialTimeThre;
    trialData.wheelTime_preTrial{i+1} = trialData.wheelTime{i}(preTrialIdx) - trialData.timeStart(i+1) + trialData.timeStart(i);
    trialData.wheelPos_preTrial{i+1} = trialData.wheelPos{i}(preTrialIdx) - trialData.wheelStart(i+1) + trialData.wheelStart(i);
end 
% extract sound period movement activity
%soundOnThre = 5000;
%bins = -preTrialThre:soundOnThre;
wheelPos_aligned = cell(1,nTrials);
%wheelPos_aligned = nan(nTrials,length(bins));
for i = 1:(nTrials)
    soundOnThre = max(5000,max(trialData.wheelTime{i}));
    %tempIdx = find(trialData.wheelTime{i} <= trialData.choiceTime(i));
    %tempIdx = find(trialData.wheelTime{i} <= soundOnThre);
    tempIdx = find(trialData.wheelTime{i});
    %tempTime = [-preTrialThre trialData.wheelTime_preTrial{i} trialData.wheelTime{i}(tempIdx) soundOnThre+1];
    tempTime = [-preTrialThre trialData.wheelTime_preTrial{i} trialData.wheelTime{i}(tempIdx) soundOnThre+1];
    if ~isempty(trialData.wheelPos_preTrial{i}); tempFirst = trialData.wheelPos_preTrial{i}(1);
    else; tempFirst = 0; end
    if ~isempty(tempIdx); tempLast = trialData.wheelPos{i}(tempIdx(end));       
    else; tempLast = 0;  end
    
    tempPos = [ tempFirst trialData.wheelPos_preTrial{i} trialData.wheelPos{i}(tempIdx) tempLast];
    
    bins = -preTrialThre:soundOnThre;
    wheelPos_aligned{i} = nan(1,length(bins));
    for j = 1:length(tempTime)-1
        tempBinIdx1 = find(tempTime(j) == bins); tempBinIdx2 = find((tempTime(j+1)-1) == bins);
        %wheelPos_aligned(i,tempBinIdx1:tempBinIdx2) = tempPos(j);
        wheelPos_aligned{i}(tempBinIdx1:tempBinIdx2) = tempPos(j);
    end
end
end

function [timePoint, wheel] = findTimeWheel(txtLine)
    splitLine = strsplit(txtLine,';'); 
    
    timeFlag = contains(splitLine,'Time');
    if any(timeFlag); tempTxt = splitLine{timeFlag}; tempTxt = strsplit(tempTxt,' = '); 
        timePoint = str2double(tempTxt{2}); end
    wheelFlag = contains(splitLine,'Wheel');
    if any(wheelFlag); tempTxt = splitLine{wheelFlag}; tempTxt = strsplit(tempTxt,' = '); 
        wheel = str2double(tempTxt{2}); end
end