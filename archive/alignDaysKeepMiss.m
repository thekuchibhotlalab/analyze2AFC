function outmat = alignDaysKeepMiss(trialData)

learningCurveBin = 50;

% convert stimulus type from 1 or 2 to 1 or -1 for plotting
stimulus = trialData.stimulus;
outmat.stimulus = fn_cell2mat(cellfun(@(x)(smoothdata(-x*2+3,'movmean',...
    learningCurveBin)),stimulus,'UniformOutput',false),1);

responseType = trialData.responseType; 

actionRate = cellfun(@(x)(smoothdata(double(x~=0),'movmean',...
    learningCurveBin)),responseType,'UniformOutput',false);

accuracy = cellfun(@(x)(smoothdata(double(x==1),'movmean',...
    learningCurveBin)),responseType,'UniformOutput',false);

accuracyCellByDay = cellfun(@(x,y)(x./y),accuracy,actionRate,'UniformOutput',false);
outmat.accuracyTop = fn_cell2mat(cellfun(@(x)(prctile(x,90)),accuracyCellByDay,'UniformOutput',false),1);
outmat.accuracyMean = fn_cell2mat(cellfun(@nanmean,accuracyCellByDay,'UniformOutput',false),1);

[actionRate, ~] = fn_cell2mat(actionRate,1);
[accuracy, cellSize] = fn_cell2mat(accuracy,1);
outmat.accuracy = accuracy./actionRate;
outmat.actionRate = actionRate;

[bias, acc_L, acc_R] = cellfun(@fn_getBias,trialData.stimulus,trialData.responseType,'UniformOutput',false);

bias_absMean = fn_removeEmptyEntryCell(bias);
bias_absMean = cellfun(@(x)(nanmean(abs(x))), bias_absMean,'UniformOutput',false);
outmat.bias_absMean = fn_cell2mat(bias_absMean,1);

accMean_L = fn_removeEmptyEntryCell(acc_L); accMean_R = fn_removeEmptyEntryCell(acc_R);
outmat.accMean_L = fn_cell2mat(cellfun(@nanmean, accMean_L,'UniformOutput',false),1);
outmat.accMean_R = fn_cell2mat(cellfun(@nanmean, accMean_R,'UniformOutput',false),1);

outmat.bias = fn_cell2mat(bias,1); 
outmat.acc_L = fn_cell2mat(acc_L,1); outmat.acc_R = fn_cell2mat(acc_R,1);

missBias = cellfun(@fn_getMissBias,trialData.stimulus,trialData.responseType,'UniformOutput',false);
outmat.missBias = fn_cell2mat(missBias,1);
%--------PROBE ORGANIZE BY DAY---------
[probeData, trialNum,probeTrialNum,trialNumNoMiss,probeTrialNumNoMiss,reinfDataBef, reinfDataAft]...
    = cellfun(@fn_getProbe,trialData.stimulus,trialData.responseType,trialData.context,'UniformOutput',false);

outmat.probeData = fn_cell2mat(probeData,1);
outmat.reinfDataBef = fn_cell2mat(reinfDataBef,1);
outmat.reinfDataAft = fn_cell2mat(reinfDataAft,1);

trialNum = fn_cell2mat(trialNum,1);
probeTrialNum = fn_cell2mat(probeTrialNum,1);

cumsumTrials = [0;cumsum(trialNum)];
probeTrialNum = probeTrialNum + cumsumTrials(1:end-1);
probeTrialNum(isnan(probeTrialNum)) = [];

dayLen = cellSize(1,:); dayLen(dayLen==0) = [];

outmat.trialNum = trialNum;
outmat.probeTrialNum = probeTrialNum;
outmat.probeDayNum = find(~isnan(probeTrialNum));
outmat.dayLen = dayLen;

end