function [tbl,RT,label,attribute,selFlag] = makeRegressorTable(mouseMega,selMouse,labelInput,attributeInput)
bias = fn_cell2mat(mouseMega.getProp('behav','field','bias','matFlag',false,'idx',selMouse),1);
day = fn_cell2mat(mouseMega.getProp('behav','field','bias','matFlag',false,'idx',selMouse),1);
choice = fn_cell2mat(mouseMega.getProp('behav','field','action','matFlag',false,'idx',selMouse),1);
reward = fn_cell2mat(mouseMega.getProp('behav','field','reward','matFlag',false,'idx',selMouse),1);
context = fn_cell2mat(mouseMega.getProp('behav','field','context','matFlag',false,'idx',selMouse),1);
[wheelDownSample, RT,selFlag] = loadWheelClusterData(mouseMega, selMouse);
label = labelInput(selFlag); label(label<=2) = 0; label(label>2) = 1; label = double(label)';
attribute = attributeInput(selFlag,:);
trials = 1:length(day);
choice(choice==2) = -1; choiceXbias = choice.*bias;
context(context==2) = 1; 

nTrialHistory = 3;
rewardH = nan(size(reward,1),nTrialHistory); 
for i = 1:nTrialHistory
    rewardH(:,i) = circshift(reward,i);
    rewardH(1:i,i) = nan;
end

% generate models that predicts the reaction time
tbl = mouseMega.loadAnimalBehav({'day','action','actionRate'},'idx',selMouse); tbl = fn_cell2mat(tbl,1);
tbl = addvars(tbl,context,bias,choiceXbias,rewardH(:,1),rewardH(:,2),rewardH(:,3),trials','NewVariableNames',...
    {'context','bias','choiceXbias','rewardH1','rewardH2','rewardH3','trials'}); 

end