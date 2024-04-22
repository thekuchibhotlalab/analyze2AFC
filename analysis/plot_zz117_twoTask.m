%% GET OBJECT FOR PT pip COHORTS
clear; 
mouse ='zz117'; 

task = 'FM_Dir_sim'; mouseObj{1} = getObj(mouse,task);
task = 'PT_pip_sim'; mouseObj{2} = getObj(mouse,task);
task = 'Two_task_sim_training'; mouseObj{3} = getObj(mouse,task);

learningCurveBin = 100;
acc_simbyday = {};bias_simbyday = {};
for i = 1:2
    %tempDay = unique(mouseObj{i}.behav.day);
    acc_simbyday{i} = mouseObj{i}.behav.acc;
    bias_simbyday{i} = mouseObj{i}.behav.bias;
end 


tempDay = unique(mouseObj{3}.behav.day);
taskFlag = mouseObj{3}.behav.stimulus == 1 | mouseObj{3}.behav.stimulus == 2;
[tempBias1, tempAccL,  tempAccR] = fn_getBias(mouseObj{3}.behav.stimulus(taskFlag),mouseObj{3}.behav.responseType(taskFlag),learningCurveBin);
acc_simbysession{1} = (tempAccL + tempAccR)/2; bias_simbysession{1} = tempBias1;
[tempBias2, tempAccL,  tempAccR]  = fn_getBias(mouseObj{3}.behav.stimulus(~taskFlag)-2,mouseObj{3}.behav.responseType(~taskFlag),learningCurveBin);
acc_simbysession{2} = (tempAccL + tempAccR)/2; bias_simbysession{2} = tempBias2;

% acc
figure; subplot(1,2,1);hold on; xlimit = max([length(acc_simbyday{1}) length(acc_simbyday{2})]);
plot(acc_simbyday{1},'LineWidth',2);  plot(acc_simbyday{2},'LineWidth',2); 
plot([1 xlimit],[0.5 0.5],'Color',[0.8 0.8 0.8],'LineWidth',2)
ylim([0.3 1])


subplot(1,2,2); hold on; xlimit = max([length(acc_simbysession{1}) length(acc_simbysession{2})]);
plot(acc_simbysession{1},'LineWidth',2);  plot(acc_simbysession{2},'LineWidth',2); 
plot([1 xlimit],[0.5 0.5],'Color',[0.8 0.8 0.8],'LineWidth',2)
ylim([0.3 1])

% bias
figure; subplot(1,2,1);hold on; xlimit = max([length(bias_simbyday{1}) length(bias_simbyday{2})]);
plot(bias_simbyday{1},'LineWidth',2);  plot(bias_simbyday{2},'LineWidth',2); 
plot([1 xlimit],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
ylim([-1 1])


subplot(1,2,2); hold on; xlimit = max([length(bias_simbysession{1}) length(bias_simbysession{2})]);
plot(bias_simbysession{1},'LineWidth',2);  plot(bias_simbysession{2},'LineWidth',2); 
plot([1 xlimit],[0 0],'Color',[0.8 0.8 0.8],'LineWidth',2)
ylim([-1 1])



function mouseObj = getObj(mouse,task)
rootPath = 'C:\Users\zzhu34\Documents\tempdata\octoData\';
loadPath = [rootPath filesep 'trialData\'];
load([loadPath mouse '.mat']);
mouseObj = wheel2AFC(trialData,'opsPath',['C:\Users\zzhu34\Documents\tempdata\octoData\opsFiles\' task ],'mouse',mouse );
mouseObj = mouseObj.removeMiss(); 
mouseObj = mouseObj.getAcc();


end