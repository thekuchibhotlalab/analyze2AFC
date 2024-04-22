%% LOAD DATA
clear;
mice ={'zz054','zz062','zz063','zz066','zz067','zz068','zz069','zz107','zz109','zz111','zz112','zz113','zz115'};
task = 'puretone'; nMouse = length(mice);

allMouse = fn_getObjPT_bin40(mice);
mouseMega = wheel2AFCmega(allMouse);

%%
biasL_accR = {};
biasL_accR_stimR = {};
biasU_accR = {};

for i = 1:13
    biasL_accR{i} =[];
    biasL_accR_stimR{i} =[]; 
    biasU_accR{i} = [];
    for j = 1:length(allMouse{i}.biasBlock.blockL.start)
        tempIdx = allMouse{i}.biasBlock.blockL.start(j):allMouse{i}.biasBlock.blockL.end(j);
        tempTable = allMouse{i}.behav(tempIdx,:);
        tempRflag = tempTable.action == 2;
        tempRflag2 = tempTable.stimulus == 2;
        tempR_acc = sum(tempTable.stimulus(tempRflag)==tempTable.action(tempRflag)) / sum(tempRflag);
        biasL_accR{i} = [biasL_accR{i} tempR_acc];
        tempR_acc2 = sum(tempTable.stimulus(tempRflag2)==tempTable.action(tempRflag2)) / sum(tempRflag2);
        biasL_accR_stimR{i} = [biasL_accR_stimR{i} tempR_acc2];        
        tempIdx = [ allMouse{i}.biasBlock.blockL.start(j)-20:allMouse{i}.biasBlock.blockL.start(j)-1 ...
            allMouse{i}.biasBlock.blockL.end(j)+1:allMouse{i}.biasBlock.blockL.end(j)+20];
        try
            tempTable = allMouse{i}.behav(tempIdx,:);
            tempUflag = tempTable.action == 2;
            tempU_acc = sum(tempTable.stimulus(tempUflag)==tempTable.action(tempUflag)) / sum(tempUflag);
            biasU_accR{i} = [biasU_accR{i} tempU_acc];
        catch
            biasU_accR{i} = [biasU_accR{i} nan];
        end
    end
end


biasR_accL = {};
biasL_accL_stimL = {};
biasU_accL = {};
for i = 1:13
    biasR_accL{i} =[];
    biasU_accL{i} =[]; 
    biasL_accL_stimL{i} = [];
    for j = 1:length(allMouse{i}.biasBlock.blockR.start)
        tempIdx = allMouse{i}.biasBlock.blockR.start(j):allMouse{i}.biasBlock.blockR.end(j);
        tempTable = allMouse{i}.behav(tempIdx,:);
        tempRflag = tempTable.action == 1;
        tempRflag2 = tempTable.stimulus == 1;
        tempL_acc = sum(tempTable.stimulus(tempRflag)==tempTable.action(tempRflag)) / sum(tempRflag);
        biasR_accL{i} = [biasR_accL{i} tempL_acc];
        tempL_acc2 = sum(tempTable.stimulus(tempRflag2)==tempTable.action(tempRflag2)) / sum(tempRflag2);
        biasL_accL_stimL{i} = [biasL_accL_stimL{i} tempL_acc2];        
        tempIdx = [ allMouse{i}.biasBlock.blockR.start(j)-20:allMouse{i}.biasBlock.blockR.start(j)-1 ...
            allMouse{i}.biasBlock.blockR.end(j)+1:allMouse{i}.biasBlock.blockR.end(j)+20];
        try
            tempTable = allMouse{i}.behav(tempIdx,:);
            tempUflag = tempTable.action == 1;
            tempU_acc = sum(tempTable.stimulus(tempUflag)==tempTable.action(tempUflag)) / sum(tempUflag);
            biasU_accL{i} = [biasU_accL{i} tempU_acc];
        catch
            biasU_accL{i} = [biasU_accL{i} nan];
        end
    end
end

biasL_accR = cellfun(@(x)(nanmean(x)),biasL_accR);
biasL_accR_stimR = cellfun(@(x)(nanmean(x)),biasL_accR_stimR);
biasU_accR = cellfun(@(x)(nanmean(x)),biasU_accR);

biasR_accL = cellfun(@(x)(nanmean(x)),biasR_accL);
biasL_accL_stimL = cellfun(@(x)(nanmean(x)),biasL_accL_stimL);
biasU_accL = cellfun(@(x)(nanmean(x)),biasU_accL);
%%

fn_plotComparison({nanmean([biasL_accR ; biasR_accL],1), nanmean([biasU_accR ; biasU_accL],1)},'paired',true);
ylim([0.5 1]);

%% plot bias curve of all animals

figure;
for i = 1:13
    subplot(7,2,i); hold on;
    plot([1 length(allMouse{i}.behav.bias)],[0 0],'Color',[0.6 0.6 0.6],'LineWidth',1.5)
    plot(allMouse{i}.behav.bias, 'LineWidth',2); ylim([-1 1]); xlim([1 length(allMouse{i}.behav.bias)])
    xlabel('Trial in training'); ylabel('Choice bias')
end

