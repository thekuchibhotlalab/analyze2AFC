function [accL_All,accR_All,missL_All,missR_All] = testProbeCSV(mouse)
sep = '\';
loadPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\probeCSV';
csvData = readtable([loadPath sep mouse '.csv']);

day = csvData.Day;

switch mouse
    case 'zz075'
        acqDayL = [10 12 13]; acqDayR = [2 3 4 5];
    case 'zz076'
        acqDayL = [3 4 5 6 7 8]; acqDayR = [9 10 11 12 13 14 15 16 17 18];
    case 'zz077'
        acqDayL = [3 4 5 6 7 8]; acqDayR = [2 3 4];
    case 'zz048'
        acqDayL = [1 2 3 4 5 6 7 8 9 10]; acqDayR = [8 9 10];
    case 'zz081'
        acqDayL = [1 3]; acqDayR = [1 2 3];
    case 'zz082'
        acqDayL = [ 2 3 4 5]; acqDayR = 1;
    case 'zz083'
        acqDayL = [10 11 12]; acqDayR = [ 3 9 10 11];
end


for i = 1:length(day)
    tempL = (csvData.PLCorr + csvData.PLInCorr + csvData.PLInMiss);
    tempR = (csvData.PRCorr + csvData.PRInCorr + csvData.PRInMiss);
    if  tempL> 10; csvData.PLCorr = csvData.PLCorr - tempL + 10; end
    if  tempR> 10; csvData.PRCorr = csvData.PRCorr - tempR + 10; end
    
end
accL = (csvData.PLCorr)./(csvData.PLCorr + csvData.PLInCorr);
accR = (csvData.PRCorr)./(csvData.PRCorr + csvData.PRInCorr);

missL = (csvData.PLCorr + csvData.PLInCorr)./(csvData.PLCorr + csvData.PLInCorr + csvData.PLInMiss);
missR = (csvData.PRCorr + csvData.PRInCorr)./(csvData.PRCorr + csvData.PRInCorr + csvData.PRInMiss);

%pre-acq L
if acqDayL(1) > csvData.Day(1)
    accLPre = accL(csvData.Day(1):acqDayL(1)-1); missLPre = missL(csvData.Day(1):acqDayL(1)-1); 
else; accLPre = nan; missLPre = nan; 
end 
%pre-acq R
if isempty(acqDayR); accRPre = accR;
elseif acqDayR(1) > csvData.Day(1)
    accRPre = accR(csvData.Day(1):acqDayR(1)-1); missRPre = missR(csvData.Day(1):acqDayR(1)-1); 
else; accRPre = nan; missRPre = nan; 
end 

% Acq L
accLAcq = accL(acqDayL); accRAcq = accR(acqDayR); missLAcq = missL(acqDayL); missRAcq = missR(acqDayR); 

%post-acq L
if acqDayL(end) < csvData.Day(end)
    accLPost = accL(acqDayL(end)+1:end); missLPost = missL(acqDayL(end)+1:end); 
else; accLPost = nan; missLPost = nan;
end 
%post-acq R
if ~isempty(acqDayR) && acqDayR(end) < csvData.Day(end)
    accRPost = accR(acqDayR(end)+1:end); missRPost = missR(acqDayR(end)+1:end); 
else; accRPost = nan; missRPost = nan; 
end 

maxLen = max([length(accLPre) length(accLAcq) length(accLPost)]);
accL_All = cat(2,fn_attachNan(accLPre,maxLen,1),fn_attachNan(accLAcq,maxLen,1),fn_attachNan(accLPost,maxLen,1));
maxLen = max([length(accRPre) length(accRAcq) length(accRPost)]);
accR_All = cat(2,fn_attachNan(accRPre,maxLen,1),fn_attachNan(accRAcq,maxLen,1),fn_attachNan(accRPost,maxLen,1));

maxLen = max([length(missLPre) length(missLAcq) length(missLPost)]);
missL_All = cat(2,fn_attachNan(missLPre,maxLen,1),fn_attachNan(missLAcq,maxLen,1),fn_attachNan(missLPost,maxLen,1));
maxLen = max([length(missRPre) length(missRAcq) length(missRPost)]);
missR_All = cat(2,fn_attachNan(missRPre,maxLen,1),fn_attachNan(missRAcq,maxLen,1),fn_attachNan(missRPost,maxLen,1));

end

