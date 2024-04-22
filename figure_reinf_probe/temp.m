behav = allMouse{1}.behav;


[probeOn, probeOff] = fn_getBlockOnOff(behav.goodProbe==1);
onIdx = find(probeOn); offIdx = find(probeOff);
% REDO the bef and aft selection since miss is removed
for i = 1:length(onIdx)
    tempProbeIdx = onIdx(i):offIdx(i); 
    obj.probe.probeData(i,:) = fn_getIdxAccBias(obj.behav.stimulus(tempProbeIdx), obj.behav.responseType(tempProbeIdx));
    tempBefIdx = tempProbeIdx - length(tempProbeIdx); tempAftIdx = tempProbeIdx + length(tempProbeIdx);
    obj.probe.reinfDataBef(i,:) = fn_getIdxAccBias(obj.behav.stimulus(tempBefIdx), obj.behav.responseType(tempBefIdx));
    obj.probe.reinfDataAft(i,:) = fn_getIdxAccBias(obj.behav.stimulus(tempAftIdx), obj.behav.responseType(tempAftIdx));
    obj.probe.trialIdx(i) = mean(tempProbeIdx);


    obj.probe.biasBef(i) = nan;
    obj.probe.biasAft(i) = nan;

end


%% dsfsdfsf
[days,nDays, dayFlag] = fn_unique(behav.day);
for i= 1:nDays
    tempBehav = behav.day(dayFlag{i},1);
    tempProbeFlag = tempBehav.context==3;
    
    
end