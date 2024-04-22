function tempObj = fn_getObjMultiTask(mice)


for i = 1:length(mice)
    ops = opsSharedArg(struct());
    ops = opsIndivArg(ops,mice{i});
    
    load([ops.datapath mice{i} '.mat']);
    
    tempObj = wheel2AFC(trialData,'ops',ops,'mouse',mice{i} );
    
    tempObj = tempObj.correctMultiTaskStim(); 
    tempObj = tempObj.removeMiss(); 
    tempObj = tempObj.getAcc();
    try; tempObj = tempObj.calculateProbe(); catch; disp('Probe is corrupt'); end
    tempObj = tempObj.getBiasBlock(0.2);
end
    
    
end

function ops = opsSharedArg(ops)
    
% ------------------ required arguments ------------------
ops.datapath = 'C:\Users\zzhu34\Documents\tempdata\octoData\trialData\';
ops.behavType = 'wheel2AFC';

ops.preprocessedInput = false;
ops.psytrackPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
ops.correctWheelTrial = true;
ops.learningCurveBin = 200;
end

function ops = opsIndivArg(ops,mouse)
    
switch mouse
    case 'zz097'; ops.correctionMultiTaskName = {'FM_Dir_pip'}; ops.selectProtocol ={'PT_pip','FM_Dir_pip'};
    case 'zz098'; ops.correctionMultiTaskName = {'FM_Dir_pip'}; ops.selectProtocol ={'PT_pip','FM_Dir_pip'};
    case 'zz101'; ops.correctionMultiTaskName = {'PT_pip'}; ops.selectProtocol ={'PT_pip','FM_Dir_Dur_pip','Two_task'};
    case 'zz107'; ops.correctionMultiTaskName = {'FM_Dir_pip'}; ops.selectProtocol ={'puretone','PT_pip','FM_Dir_pip'};
    case 'zz111'; ops.correctionMultiTaskName = {'FM_Dir_pip'}; ops.selectProtocol ={'puretone','PT_pip','FM_Dir_pip'};
    case 'zz115'; ops.correctionMultiTaskName = {'FM_Dir_pip'}; ops.selectProtocol ={'puretone','PT_pip','FM_Dir_pip'};
    case 'zz121'; ops.correctionMultiTaskName = {'PT_pip'}; ops.selectProtocol ={'PT_pip','FM_Dir_pip','Two_task'};
    case 'zz122'; ops.correctionMultiTaskName = {'PT_pip'}; ops.selectProtocol ={'PT_pip','FM_Dir_pip','Two_task'};
    case 'msb01'; ops.selectProtocol ={'PT_pip','FM_Dir_Dur_pip'};
end

end