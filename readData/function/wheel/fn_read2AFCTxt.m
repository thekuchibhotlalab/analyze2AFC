function [wheelSoundOn,wheelPreSound,hit] = fn_read2AFCTxt(txtFilename)

nTrial = 0; soundOnFlag = false; hit = [];
wheelSoundOn = {}; wheelPreSound = {};  
tempSound = {};  tempPreSound = {}; 
if isfile(txtFilename)
    txtData = splitlines(fileread(txtFilename));
    for i = 1:length(txtData)                  
        if ~isempty(txtData{i}) && ~isletter(txtData{i}(1)) % Wheel position if the first char is not letter
            txtDataSplit = strsplit(txtData{i});
            txtDataSplitMat = cell2mat(cellfun(@str2double,txtDataSplit,'UniformOutput',false));
            if sum(isnan(txtDataSplitMat)) ~=0 
                % Get rid of all entries after the first nan (i.e. letter)
                % Due to the exception of 'Sound 1 on!'
                firstNan = fn_findFirst(isnan(txtDataSplitMat));
                txtDataSplitMat(firstNan:end) = [];
            end
            if soundOnFlag
                tempSound{end+1} = txtDataSplitMat;
            else
                tempPreSound{end+1} = txtDataSplitMat;
            end      
        end
        
        
        % Since 'Sonud on' and 'Sound off' always at the end of a line,
        % record wheel first, then check if song is on or off
        if contains(txtData{i},'Sound')
            if contains(txtData{i},'on')
                soundOnFlag = true; nTrial = nTrial+1;
                wheelPreSound{nTrial} = tempPreSound; 
                tempPreSound = {}; tempSound = {};
            elseif contains(txtData{i},'off')
                if nTrial == 0; nTrial = nTrial + 1; end % In case of sound-on not recorde in trial 1
                soundOnFlag = false; hit(nTrial) = 0;
                if ~isempty(tempSound); wheelSoundOn{nTrial} = tempSound{1}; 
                else; wheelSoundOn{nTrial} = []; end
                if length(tempSound) >1
                    disp('WARNING -- multiple movement in sound on'); 
                end
                tempSound = {}; tempPreSound = {};
            end
        elseif contains(txtData{i},'choice reached')
            % Choice reached always follow sound off; only record action
            hit(nTrial) = 1;
        end

    end

    
else
    error(['Error - txtFilename ' txtFilename ' does not exist!']);
end