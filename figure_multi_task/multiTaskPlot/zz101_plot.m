%clear;

%load('zz101_behav.mat');
behav = mouseObj.behav;

%create a table for only two-tasks days
two_task_table = behav(behav.trainingType=="Two_task", :);
toDelete = two_task_table.day > 75;
two_task_table(toDelete,:) = [];

dayCell = {};
b1_b2_rt={};
b1_b2_acc={};
b1_b2_bias={};
b3_rt={};
b3_acc={};
b3_bias={};
first10_acc={};
before10_acc={};
first10_bias={};
before10_bias={};
first10_rt={};
before10_rt={};


for d = 66:75 % max is 53-75
    if d ~=61
        %inside each day
        date_table=two_task_table(two_task_table.day==d,:);
        %add a column for block 
        date_table.block = zeros(height(date_table),1);
        
        %blocks
        bin=1;
        while (bin+9)<height(date_table)
            u = unique(date_table.stimulus(bin:bin+9));
            if size(u)==size([1;2])
                if u== [1;2]
                date_table{bin:bin+9,'block'}=1;
            
                elseif u== [3;4]
                date_table{bin:bin+9,'block'}=2;
                end 
            else 
                date_table{bin:bin+9,'block'}=3;
            end
            bin=bin+10;
        end
        date_table{bin:end,'block'}=date_table.block(bin-1);

%         %accuracy
%         start_b2_acc = sum(date_table{41:50,'reward'})/10;
%         late_b2_acc=sum(date_table{71:80,'reward'})/10;
%         start_b1_acc = sum(date_table{161:170,'reward'})/10;
%         late_b1_acc=sum(date_table{191:200,'reward'})/10;
%         T(i,:)= {start_b2_acc,late_b2_acc, start_b1_acc, late_b1_acc};
% 
%         i=i+1;

        %fix transitions between block 1 and 2
        [onIdx, offIdx,~,~,~] = fn_getBlockOnOff(date_table.block==3);
        blockLen = (offIdx-onIdx)+1; 
        for nBlock = 1:length(blockLen)
            tempPreBlock = date_table.block(onIdx(nBlock)-1); correction = true; 
            tempCount = 0;
            while correction
                tempIdx = onIdx(nBlock)+tempCount;
                if (tempPreBlock == 1 && date_table.stimulus(tempIdx) <=2) || (tempPreBlock == 2 && date_table.stimulus(tempIdx) >=3 )
                    date_table.block(tempIdx) = tempPreBlock; tempCount = tempCount + 1;
                else
                    correction = false;
                end
            end

            if offIdx(nBlock)+1<length(date_table.block)
                tempPostBlock = date_table.block(offIdx(nBlock)+1); correction = true; 
            else
                correction = false; 
            end
            tempCount = 0;
            while correction
                tempIdx = offIdx(nBlock)-tempCount;
                if (tempPostBlock == 1 && date_table.stimulus(tempIdx) <=2) || (tempPostBlock == 2 && date_table.stimulus(tempIdx) >=3 )
                    date_table.block(tempIdx) = tempPostBlock; tempCount = tempCount + 1;
                else
                    correction = false;
                end
            end
        end
       dayCell = [dayCell;date_table];

    [onIdx, offIdx,blockIdx,blockArea,blockMax] = fn_getBlockOnOff(date_table.block==3);
    disp(offIdx - onIdx+1);disp(onIdx)
     
     %accuracy and bias for block 1 and 2 day by day, add to a cell aray 
     [bias, acc, ~, ~, ~, ~] = fn_getAccBias(date_table.stimulus(date_table.block == 1|2), date_table.responseType(date_table.block == 1|2)==1, date_table.responseType(date_table.block ==1|2)==0);
     b1_b2_acc=[b1_b2_acc; acc ];
     b1_b2_bias=[b1_b2_bias; bias ];
    
     %accuracy and bias for block 3 day by day, add to a cell aray 
     [bias, acc, ~, ~, ~, ~] = fn_getAccBias(date_table.stimulus(date_table.block == 3), date_table.responseType(date_table.block == 3)==1, date_table.responseType(date_table.block ==3)==0);
     b3_acc=[b3_acc; acc ];
     b3_bias=[b3_bias; bias ];

     %reaction time 
     b1_b2_rt =[b1_b2_rt; mean(date_table.reactionTime(date_table.block == 1|2))];
     b3_rt =[b3_rt; mean(date_table.reactionTime(date_table.block == 3))];

     %accuracy and bias for first 10 trials of interleaved
     acc_total=0;
     bias_total=0;
     for i=1:length(onIdx)
         [bias, acc, ~, ~, ~, ~] = fn_getAccBias(date_table.stimulus(onIdx(i):onIdx(i)+9), date_table.responseType(onIdx(i):onIdx(i)+9)==1, date_table.responseType(onIdx(i):onIdx(i)+9)==0);
         acc_total=acc_total+acc;
         bias_total=bias_total+bias;
     end       
     first10_acc=[first10_acc;acc_total/length(onIdx)];
     first10_bias=[first10_bias;bias_total/length(onIdx)];
     first10_rt = [first10_rt;mean(date_table.reactionTime(onIdx(i):onIdx(i)+9))];

     %accuracy and bias for 10 trials before interleaved
     acc_total=0;
     bias_total=0;
     for i=1:length(onIdx)
        [bias, acc, ~, ~, ~, ~] = fn_getAccBias(date_table.stimulus(onIdx(i)-10:onIdx(i)-1), date_table.responseType(onIdx(i)-9:onIdx(i))==1, date_table.responseType(onIdx(i)-9:onIdx(i))==0);
        acc_total=acc_total+acc;
        bias_total=bias_total+bias;
     end
     before10_acc=[before10_acc;acc_total/length(onIdx)];
     before10_bias=[before10_bias;bias_total/length(onIdx)];
     before10_rt = [before10_rt;mean(date_table.reactionTime(onIdx(i)-9:onIdx(i)))];
    end
   
end

%% SECTION Accuracy comparison beween b1+b2 and b1/b2
figure;
p = fn_plotComparison({cell2mat(b1_b2_acc),cell2mat(b3_acc)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks([1 2]); xticklabels({'b1+b2','b1/b2'});title(['p=' num2str(p)])

%% SECTION Bias comparison beween b1+b2 and b1/b2
figure;
p = fn_plotComparison({cell2mat(b1_b2_bias),cell2mat(b3_bias)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks([1 2]); xticklabels({'b1+b2','b1/b2'});title(['p=' num2str(p)])

%% SECTION Reaction Time comparison beween b1+b2 and b1/b2
figure;
p = fn_plotComparison({cell2mat(b1_b2_rt),cell2mat(b3_rt)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks([1 2]); xticklabels({'b1+b2','b1/b2'});title(['p=' num2str(p)])

%% SECTION Accuracy comparison beween before and first 10 trials of interleave
figure;
p = fn_plotComparison({cell2mat(before10_acc),cell2mat(first10_acc)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks([1 2]); xticklabels({'block','interleave'}); title(['p=' num2str(p)])

%% SECTION Bias comparison beween before and first 10 trials of interleave
figure;
p = fn_plotComparison({cell2mat(before10_bias),cell2mat(first10_bias)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks([1 2]); xticklabels({'before','first'});title(['p=' num2str(p)])
%% SECTION Reaction Time comparison beween before and first 10 trials of interleave
figure;
p = fn_plotComparison({cell2mat(before10_rt),cell2mat(first10_rt)},'paired',true,'compType','errorbarWithDot');
xlim([0.6 2.4]); xticks([1 2]); xticklabels({'before','first'});title(['p=' num2str(p)])




