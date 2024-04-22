function [label4,attributesOrig] = fn_loadWheelCluster()

pythonPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\wheelData\';

load([pythonPath filesep 'wheelAttributes_allMouse.mat' ]);
load([pythonPath filesep 'wheel_cluster_allMouse.mat' ]); 
load([pythonPath filesep 'wheel_cluster_allMouse_tsne.mat' ]);

% mice ={'zz107','zz109','zz111','zz112','zz113','zz115'};
% task = 'puretone'; nMouse = length(mice);
% allMouse = fn_getObjPT_bin30(mice);
% mouseMega = wheel2AFCmega(allMouse);

label4 = double(label4 + 1); label5 = double(label5 + 1); 
attributesOrig = tempAttributes;
%label4(tempNan) = nan; label5(tempNan) = nan;
%attributesOrig(tempNan,:) = nan;
end 