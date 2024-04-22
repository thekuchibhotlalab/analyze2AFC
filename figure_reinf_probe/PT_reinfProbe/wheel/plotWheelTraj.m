function plotWheelTraj(mouse,selectProtocol)
%plotWheelTraj(mouse,selectProtocol)
%InputArg:
%   mouse: 'zz067'
%   selectProtocol: {'puretone'} (default)

if nargin == 1; selectProtocol = {'puretone'}; end

global rootPath sep;
sep = '\';
rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath sep 'trialData\']; 
figPath = [rootPath sep 'figure\' sep 'wheelTraj' sep]; mkdir(figPath);

%selectProtocol = {'FM_One_Oct','FM_One_Oct_Prob','FM_Half_Oct','FM_Half_Oct_Prob'};
trialData = fn_selectProtocol(loadPath, mouse, selectProtocol);

wheelTraj = cellfun(@getWheelTraj, trialData.wheelSoundOn, trialData.wheelSoundOnCheckFlag,trialData.action,'UniformOutput',false);
wheelTrajL = cellfun(@(x,y)(nanmean(x(y==1,:),1)), wheelTraj, trialData.action,'UniformOutput',false);
wheelTrajL = fn_cell2matFillNan(wheelTrajL);
wheelTrajR = cellfun(@(x,y)(nanmean(x(y==2,:),1)), wheelTraj, trialData.action,'UniformOutput',false);
wheelTrajR = fn_cell2matFillNan(wheelTrajR);
wheelTrajM = cellfun(@(x,y)(nanmean(x(y==0,:),1)), wheelTraj, trialData.action,'UniformOutput',false);
wheelTrajM = fn_cell2matFillNan(wheelTrajM);


wheelTrajLP = cellfun(@(x,y,z)(nanmean(x(y==1 & z==3,:),1)), wheelTraj, trialData.action, trialData.context,'UniformOutput',false);
wheelTrajLP = fn_cell2matFillNan(wheelTrajLP);

wheelTrajRP = cellfun(@(x,y,z)(nanmean(x(y==2 & z==3,:),1)), wheelTraj, trialData.action, trialData.context,'UniformOutput',false);
wheelTrajRP = fn_cell2matFillNan(wheelTrajRP);

[out1, out2] = fn_sqrtInt(length(trialData.action));
f = fn_figureSmartDim('hSize',0.1*out1,'widthHeightRatio',out2/out1,'visible','off'); 

for i = 1:length(trialData.action)
    subplot_tight(out1,out2,i,0.04); hold on
    plot(wheelTrajL(i,1:100),'Color',matlabColors(1));
    plot(wheelTrajLP(i,1:100),'Color',matlabColors(1)*0.5 + [1 1 1] * 0.5);
    plot(wheelTrajR(i,1:100),'Color',matlabColors(2));
    plot(wheelTrajRP(i,1:100),'Color',matlabColors(2)*0.5 + [1 1 1] * 0.5);
    ylim([-100 100])
    if mod(i,out2) ~= 1; yticks([]); end
    if i <= (out1-1)*out2; xticks([]); end
end
saveas(f,[figPath sep mouse '.png']);
close all;
end


function wheelTraj = getWheelTraj(rawWheelTraj, checkFlag, action)
%getWheelTraj - Description
%
% Syntax: wheelTraj = getWheelTraj(rawWheelTraj)
%
% Long description
if ~isempty(rawWheelTraj) && checkFlag
    wheelTraj = rawWheelTraj; 
    wheelTraj(rawWheelTraj>0) = 1;
    wheelTraj(rawWheelTraj<0) = -1;
    wheelTraj = cumsum(wheelTraj,2);

else
    wheelTraj = nan(length(action),100);
end

end

