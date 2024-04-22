function allMouse = fn_getObjTask1(mice)

allMouse = {};
for i = 1:length(mice)
    ops = opsSharedArg(struct());
    ops = opsIndivArg(ops,mice{i});
    
    load([ops.datapath mice{i} '.mat']);
    
    tempObj = wheel2AFC(trialData,'ops',ops,'mouse',mice{i} );
    
    tempObj = tempObj.correctMultiTaskStim(); 
    tempObj = tempObj.removeMiss(); 
    %tempObj = tempObj.getAcc();
    %tempObj = tempObj.calculateProbe();
    %tempObj = tempObj.getBiasBlock(0.2);
    allMouse{i} = tempObj;
end
    
    
end

function ops = opsSharedArg(ops)
    
% ------------------ required arguments ------------------
ops.datapath = 'C:\Users\zzhu34\Documents\tempdata\octoData\trialData\';
ops.behavType = 'wheel2AFC';

ops.preprocessedInput = false;
ops.psytrackPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\psyTrackData\psyTrackFit\';
ops.correctWheelTrial = true;
ops.learningCurveBin = 500;
    
end

function ops = opsIndivArg(ops,mouse)
    
switch mouse
    case 'zz054'; ops.selectProtocol ={'puretone'};
    case 'zz062'; ops.selectProtocol ={'puretone'};
    case 'zz063'; ops.selectProtocol ={'puretone'};
    case 'zz066'; ops.selectProtocol ={'puretone'};
    case 'zz067'; ops.selectProtocol ={'puretone'};
    case 'zz068'; ops.selectProtocol ={'puretone'};
    case 'zz069'; ops.selectProtocol ={'puretone'}; 
    case 'zz097'; ops.selectProtocol ={'PT_pip'};
    case 'zz098'; ops.selectProtocol ={'PT_pip'};
    case 'zz101'; ops.selectProtocol ={'FM_Dir_Dur_pip'};
    case 'zz102'; ops.selectProtocol ={'FM_Dir_Dur_pip'};
    case 'zz105'; ops.selectProtocol ={'FM_Dir_Dur_pip'};
    case 'zz107'; ops.selectProtocol ={'puretone'};
    case 'zz109'; ops.selectProtocol ={'puretone'};
    case 'zz111'; ops.selectProtocol ={'puretone'};
    case 'zz112'; ops.selectProtocol ={'puretone'};
    case 'zz113'; ops.selectProtocol ={'puretone'};
    case 'zz115'; ops.selectProtocol ={'puretone'};
    case 'zz121'; ops.selectProtocol ={'FM_Dir_pip'};
    case 'zz122'; ops.selectProtocol ={'FM_Dir_pip'};
    case 'zz123'; ops.selectProtocol ={'FM_Dir_pip'};

end

end