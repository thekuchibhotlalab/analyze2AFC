function [probeFlag, probeStartFlag, reinfBefFlag, reinfAftFlag] = fn_getProbeFlag(ctxt,nDay)
    
goodProbeFlag = true; probeIdx = find(ctxt==3); 
probeFlag = zeros(size(ctxt)); reinfBefFlag = zeros(size(ctxt)); reinfAftFlag = zeros(size(ctxt)); 
tempProbeStartIdx = zeros(size(ctxt)); probeStartFlag = zeros(size(ctxt));
if ~isempty(probeIdx)
    if mod(probeIdx,10) ~= 0; goodProbeFlag = false; disp(['Bugged Probe Trials=' int2str(probeIdx)]); end

    if (probeIdx(end)- probeIdx(1))==19; reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20; 
    else; reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
    end   
    reinfAftIdx(reinfAftIdx > length(ctxt)) = []; reinfBefIdx(reinfBefIdx <=0 ) = [];
    if goodProbeFlag
        probeFlag(probeIdx) = nDay; reinfBefFlag(reinfBefIdx) = nDay; reinfAftFlag(reinfAftIdx) = nDay; 
        tempProbeStartIdx(probeIdx-1) = nDay; tempIdx = (tempProbeStartIdx - probeFlag) == nDay;
        probeStartFlag(tempIdx==1) = nDay;
    end
end

%{
if ~isempty(probeIdx)
    if mod(nProbeRaw,10) == 2 && nProbeRaw==22
        probeIdx([1,12]) = [];
        reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        disp('nProbe = 22')
    elseif mod(nProbeRaw,10) == 1 && (nProbeRaw==11 || nProbeRaw==21)
        probeIdx(1) = [];
        if nProbeRaw == 11
            reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        elseif nProbeRaw == 21
            reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20;
        end
        disp('nProbe = 11 or 21')
        
    elseif nProbeRaw==20
        if (nProbeRaw(end)- nProbeRaw(1))~=19
            reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        else 
            reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20;
        end   
        disp('nProbe = 20')
    elseif nProbeRaw==10
        reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
    else
        disp(['WARNING -- Exception in Probe! nProbe = ' int2str(nProbeRaw)]);
        probeTrialNum = nan; probeTrialNumNoMiss = nan;
        goodProbeFlag = false;
    end
    if goodProbeFlag
        probeData = fn_getIdxAccBias(stim, response, probeIdx);
        reinfDataBef = fn_getIdxAccBias(stim, response, reinfBefIdx);
        reinfDataAft = fn_getIdxAccBias(stim, response, reinfAftIdx);

        probeTrialNum = round(mean(probeIdx));

        probeTrialNumNoMiss = round(mean(trialNumNoMissIdx(probeIdx)));
    end
else
    probeTrialNum = nan; probeTrialNumNoMiss = nan;
end


probeFlag(probeIdx) = nDay;






 nProbeRaw = length(probeIdx);

trialNumNoMissIdx = cumsum((response~=0));

%probeData = []; reinfDataBef = []; reinfDataAft = [];
probeData = nan(1,6); reinfDataBef = nan(1,6); reinfDataAft = nan(1,6);
goodProbeFlag = true;
% Two probe sessions 
if ~isempty(probeIdx)
    if mod(nProbeRaw,10) == 2 && nProbeRaw==22
        probeIdx([1,12]) = [];
        reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        disp('nProbe = 22')
    elseif mod(nProbeRaw,10) == 1 && (nProbeRaw==11 || nProbeRaw==21)
        probeIdx(1) = [];
        if nProbeRaw == 11
            reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        elseif nProbeRaw == 21
            reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20;
        end
        disp('nProbe = 11 or 21')
        
    elseif nProbeRaw==20
        if (nProbeRaw(end)- nProbeRaw(1))~=19
            reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
        else 
            reinfBefIdx = probeIdx-20; reinfAftIdx = probeIdx+20;
        end   
        disp('nProbe = 20')
    elseif nProbeRaw==10
        reinfBefIdx = probeIdx-10; reinfAftIdx = probeIdx+10;
    else
        disp(['WARNING -- Exception in Probe! nProbe = ' int2str(nProbeRaw)]);
        probeTrialNum = nan; probeTrialNumNoMiss = nan;
        goodProbeFlag = false;
    end
    if goodProbeFlag
        probeData = fn_getIdxAccBias(stim, response, probeIdx);
        reinfDataBef = fn_getIdxAccBias(stim, response, reinfBefIdx);
        reinfDataAft = fn_getIdxAccBias(stim, response, reinfAftIdx);

        probeTrialNum = round(mean(probeIdx));

        probeTrialNumNoMiss = round(mean(trialNumNoMissIdx(probeIdx)));
    end
else
    probeTrialNum = nan; probeTrialNumNoMiss = nan;
end
if isempty(ctxt); trialNum = nan; trialNumNoMiss = nan; %probeTrialNum = []; probeTrialNumNoMiss = []; 
else; trialNum = length(ctxt); trialNumNoMiss = max(trialNumNoMissIdx);
end 
%}
end