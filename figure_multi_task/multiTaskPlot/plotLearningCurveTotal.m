%% SECTION 1
clear; 
%mice ={'msb01'};
%mice ={'zz101'}; 
%mice ={'msb08'}; 
mice ={'th02'}; 
%mice ={'sz06'};
mouseObj = fn_getObjMultiTask(mice);

%% SECTION 2

[accMulti,biasMulti,~,accL,accR] = mouseObj.getAccMultiTask; 

% plotting
figure; subplot(3,1,1); hold on; 
xlimm = [0 length(accMulti{1})];
plot(xlimm, [0.5 0.5],'Color',[0.8 0.8 0.8]);
plot(accMulti{1},'LineWidth',2,'Color',matlabColors(1)); plot(accMulti{2},'LineWidth',2,'Color',matlabColors(2));
xlim(xlimm); ylabel('accruacy')

subplot(3,1,2);  hold on; 
xlimm = [0 length(biasMulti{1})];
plot(xlimm, [0 0],'Color',[0.8 0.8 0.8])
plot(biasMulti{1},'LineWidth',2,'Color',matlabColors(1)); plot(biasMulti{2},'LineWidth',2,'Color',matlabColors(2))
xlim(xlimm); ylabel('bias')

subplot(3,1,3);  hold on; 
xlimm = [0 length(biasMulti{1})];
plot(xlimm, [0 0],'Color',[0.8 0.8 0.8])
plot(mouseObj.behav.actionRate,'LineWidth',2,'Color',matlabColors(1)); 
xlim(xlimm); ylabel('actionRate')

%% PLOT 
mouseList = {'sz06','sz07','sz08','sz10'};
stimStruct = struct(); stimStruct.sz06.choiceDir = [-1 -1 1 1]; stimStruct.sz06.dur = [20,80,50,50]; stimStruct.sz06.dir = [1,-1,0,0];
stimStruct.sz06.freqMean = [8,8,6,10]; stimStruct.sz06.freqStart = [4,16,6,10]; stimStruct.sz06.freqEnd = [16,4,6,10];

stimStruct.sz07.choiceDir = [-1 -1 1 1];stimStruct.sz07.dur = [80,20,50,50]; stimStruct.sz07.dir = [1,-1,0,0];
stimStruct.sz07.freqMean = [8,8,10,6]; stimStruct.sz07.freqStart = [4,16,10,6]; stimStruct.sz07.freqEnd = [16,4,10,6];

stimStruct.sz08.choiceDir = [-1 -1 1 1];stimStruct.sz08.dur = [20,80,50,50]; stimStruct.sz08.dir = [-1,1,0,0];
stimStruct.sz08.freqMean = [8,8,6,10]; stimStruct.sz08.freqStart = [16,4,6,10]; stimStruct.sz08.freqEnd = [4,16,6,10];

stimStruct.sz10.choiceDir = [-1 -1 1 1];stimStruct.sz10.dur = [20,80,50,50]; stimStruct.sz10.dir = [-1,1,0,0];
stimStruct.sz10.freqMean = [8,8,10,6]; stimStruct.sz10.freqStart = [16,4,10,6]; stimStruct.sz10.freqEnd = [4,16,10,6];

%selTrial = {1501:2000,851:1250,1001:2000,1801:2100}; % periods of good PT performance
selTrial = {4501:5000,3301:3600,2501:3300,1801:2100}; % periods of longest stable performance

figure; plotAccRearrange(accL,accR,mice{1},selTrial{strcmp(mouseList,mice{1})},stimStruct);


%% plot for kishore, msb08 and msb09
accMulti{2}(1:8900) = nan;
figure; xSel = 5801:11000; 
subplot(1,2,1); hold on; 
plot([0 (xSel(end)-xSel(1)+1)], [0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8]);
plot((xSel-xSel(1)+1),accMulti{1}(xSel),'LineWidth',2); 
plot((xSel-xSel(1)+1),accMulti{2}(xSel),'LineWidth',2); 
legend({['Chance' newline 'performance'],'Task 1','Task 2'}); title('Sequential training of task 1 and 2')
xticks((xSel(1)-xSel(1)):1000:(xSel(end)-xSel(1)+1)); xlim([0 xSel(end)-xSel(1)+1]); xlabel('Trials');ylabel('Accuracy')
ylim([0.25 1]); yticks(0.25:0.25:1)

subplot(1,2,2); hold on; xSel = 16001:17500; 
plot([0 1500], [0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8]);
plot(accMulti{1}(xSel),'LineWidth',2); plot(accMulti{2}(xSel),'LineWidth',2)
plot([0 (xSel(end)-xSel(1)+1)], [0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8]);
legend({['Chance' newline 'performance'],'Task 1','Task 2'}); title('Interleaved performance of task 1 and 2')
ylim([0.25 1]); yticks(0.25:0.25:1); xlabel('Trials');ylabel('Accuracy')



figure; xSel = 5801:19000; 
hold on; 
plot([0 (xSel(end)-xSel(1)+1)], [0.5 0.5],'LineWidth',2,'Color',[0.8 0.8 0.8]);
plot((xSel-xSel(1)+1),accMulti{1}(xSel),'LineWidth',2); 
plot((xSel-xSel(1)+1),accMulti{2}(xSel),'LineWidth',2); 
legend({['Chance' newline 'performance'],'Task 1','Task 2'}); title('Sequential training of task 1 and 2')
xticks((xSel(1)-xSel(1)):1000:(xSel(end)-xSel(1)+1)); xlim([0 xSel(end)-xSel(1)+1]); xlabel('Trials');ylabel('Accuracy')
ylim([0.25 1]); yticks(0.25:0.25:1)

%%
function plotAccRearrange(accL,accR,mouse,selTrial,stimStruct)
tempAcc = [1-mean(accL{1}(selTrial)) mean(accR{1}(selTrial)) 1-mean(accL{2}(selTrial)) mean(accR{2}(selTrial))];
allFields = fields(stimStruct.(mouse)); nFields = length(allFields);
for i = 1:nFields
    subplot(1,6,i);
    [tempSort,tempIdx] = sort(stimStruct.(mouse).(allFields{i}),'ascend'); bar(tempAcc(tempIdx))
    [ tempAxis, b, c] = unique(tempSort);


    xticks(1:length(tempSort)); xticklabels(strsplit(num2str(tempSort)))
    title(allFields{i}); 
    if i==1; ylabel('p(choose R)'); end
end


end
