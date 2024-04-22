%% LOAD DATA
clear; 
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
learningCurveBin = 40;
mouseMega = wheel2AFCmega(allMouse);

%%
plotBlockTransition(mouseMega)
%%

for i = 1:mouseMega.nMouse
    figure; windowSize = 400;
    tempProbeIdx = mouseMega.mouseCell{i}.probe.onIdx;
    for j = 1:20
        subplot(4,5,j); hold on; 
        tempStart = (j-1)*windowSize+1:1:j*windowSize; 
        if tempStart(end) < size(mouseMega.mouseCell{i}.behav,1)
            probeSelIdx = tempProbeIdx > tempStart(1) & tempProbeIdx < tempStart(end);
    
            scatter(mouseMega.mouseCell{i}.behav.acc(tempStart),...
                mouseMega.mouseCell{i}.behav.bias(tempStart),10,matlabColors(1),'filled');
            plot(mouseMega.mouseCell{i}.behav.acc(tempStart),...
                mouseMega.mouseCell{i}.behav.bias(tempStart),'Color',matlabColors(1));
            if any(probeSelIdx)
                tempProbe = sum(mouseMega.mouseCell{i}.probe.probeData(probeSelIdx,:),1);
                tempReinf = sum(mouseMega.mouseCell{i}.probe.befData(probeSelIdx,:) + mouseMega.mouseCell{i}.probe.aftData(probeSelIdx,:),1);
                [pBias,pAcc] = fn_getAccBiasByCount(tempProbe);
                [rBias,rAcc] = fn_getAccBiasByCount(tempReinf);
                scatter([pAcc],[pBias],20,matlabColors(2),'filled');
                scatter([rAcc],[rBias],20,matlabColors(6),'filled');
            end
            xlim([0.2 1]); ylim([-1 1])
        end
    end 

end

%%

for i = 1:mouseMega.nMouse
    figure; windowSize = 200;
    tempProbeIdx = mouseMega.mouseCell{i}.probe.onIdx;
    for j = 1:12
        subplot(3,4,j); hold on; 
        tempStart = (j-1)*windowSize+1:2:j*windowSize; 
        
        probeSelIdx = tempProbeIdx > tempStart(1) & tempProbeIdx < tempStart(end);
        try
            scatter(mouseMega.mouseCell{i}.behav.acc1(tempStart),...
                mouseMega.mouseCell{i}.behav.acc2(tempStart),10,matlabColors(1),'filled');
            plot(mouseMega.mouseCell{i}.behav.acc1(tempStart),...
                mouseMega.mouseCell{i}.behav.acc2(tempStart),'Color',matlabColors(1));
            if ~isempty(probeSelIdx)
                tempProbe = sum(mouseMega.mouseCell{i}.probe.probeData(probeSelIdx,:));
                tempReinf = sum(mouseMega.mouseCell{i}.probe.befData(probeSelIdx,:) + mouseMega.mouseCell{i}.probe.aftData(probeSelIdx,:));
                [pBias,pAcc,~,pAcc1,pAcc2] = fn_getAccBiasByCount(tempProbe);
                [rBias,rAcc,~,rAcc1,rAcc2] = fn_getAccBiasByCount(tempReinf);
                scatter([pAcc1],[pAcc2],20,matlabColors(2),'filled');
                scatter([rAcc1],[rAcc2],20,matlabColors(6),'filled');
            end
        catch
        end
        xlim([0 1]); ylim([0 1])
    end 

end