clear;

load('zz101_behav.mat');

%create a table for only two-tasks days
two_task_table = behav(behav.trainingType=="Two_task", :);
toDelete = two_task_table.day > 75;
two_task_table(toDelete,:) = [];

dayCell = {};
for d = 53:75
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
        date_table{bin:end,'block'}=date_table.block(bin-9);

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
%              if blockLen(nBlock) < 20
%                  for j=onIdx(nBlock):offIdx(nBlock)
%                      if date_table.stimulus(j)==1 || date_table.stimulus(j)==2
%                          date_table.block(j)=1;
%                      else
%                          date_table.block(j)=2;
%                      end
% 
%                  end
%                  
%              end 

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
    [onIdx, offIdx,blockIdx,blockArea,blockMax] = fn_getBlockOnOff(date_table.block==3);
    disp(offIdx - onIdx+1);disp(onIdx)

    dayCell = [dayCell;date_table];
    end
   
end








